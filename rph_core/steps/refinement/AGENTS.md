# rph_core/steps/refinement/AGENTS.md

## OVERVIEW

V4 S3/S4 unified refinement engine. A single `RefinementEngine` class driven by
`FidelityProfile` implements both S3 (low-fidelity B97-3c OPT/OptTS + r2SCAN-3c
SP) and S4 (high-fidelity M062X OPT/OptTS/FREQ + wB97M-V SP) stages. Replaces
the legacy `StageCalculator` + `LowLevelEngine` + `HighLevelEngine` triple with
a 3-pass DAG: preflight → primary → rescue → canonical.

All stage differences come from `FidelityProfile` parameters, not `if stage == "S3"/"S4"` branches.

## WHERE TO LOOK

| File | Role |
|------|------|
| `engine.py` | `RefinementEngine` — 3-pass DAG orchestration (~2400 lines) |
| `models.py` | `StructureRequest`, `PreflightOutcome`, `Pass1Outcome` dataclasses |
| `manifest_io.py` | `write_refinement_manifest()` / `read_refinement_manifest()` with legacy `s3_low_level_v3` / `s4_high_level_v4` adapters |
| `__init__.py` | Public exports: `RefinementEngine`, `FidelityProfile`, `StructureRequest`, `PreflightOutcome`, `Pass1Outcome` |
| `../fidelity_profile.py` | `FidelityProfile` dataclass — all S3/S4 differences as frozen parameters |

## 3-PASS DAG

```
Pass 0 (Preflight):
  validate XYZ → charge/multiplicity → forming-bonds range check →
  structure directory layout → AttemptRecorder per structure →
  provenance.json per structure
  ↓
Pass 1 (Primary):
  warmup (role-based: INT/TS always, precursor/product never) →
  primary Opt/OptTS (per-role initial_hessian: precursor/product →
  "model"; INT/TS → "calculate") → independent Freq on opt geometry →
  initial TS/INT/minimum classification → mode alignment score
  ↓  GLOBAL BARRIER (queue-based parallel scheduling)
Pass 2 (Rescue DAG):
  TS rescue:
    L1: ReadHess + target mode (requires validated freq Hessian)
    L2: RecalcHess + trust radius
  INT rescue:
    L1: TightOpt + CalcHess
    L2: Mode displacement (± sign along imaginary mode)
    L3: RecalcHess=5 + TightOpt
  Minimum rescue:
    L1: TightOpt (precursor/product with imaginary modes)
    L2: Mode displacement (± sign)
  TS endpoint discovery:
    IRC (only when TS validated as valid_target_ts)
  ↓
Pass 3 (Canonical):
  canonical selection (sort candidates by stationary_rank → mode_rank →
  hessian_index → converged → alignment_score → gradient_norm) →
  canonical.xyz materialized → Freq on canonical geometry (reuses primary
  Freq if geometry unchanged) → SP → composite thermochemistry →
  property-level ml_usability (9-key dict)
```

## MANIFEST SCHEMA

`refinement_manifest_v1` schema written by `write_refinement_manifest()`:

```python
{
    "schema_version": "refinement_manifest_v1",
    "stage": "S3" | "S4",
    "fidelity": "low" | "high",
    "profile_id": "b97_3c_r2scan_3c_v1" | "m062x_wb97mv_v1",
    "run_id": str | None,
    "structures": [
        {
            "id": str,               # e.g. "product_major_ts"
            "role": str,             # "precursor" | "product" | "intermediate" | "ts"
            "kind": str,             # "minimum" | "ts"
            "opt_status": str,       # "complete" | "failed" | "not_run"
            "opt_xyz": str | None,
            "opt_energy_hartree": float | None,
            "frequency_status": str,
            "frequencies_cm1": list[float],
            "ts_classification": dict | None,
            "int_classification": dict | None,
            "minimum_classification": dict | None,
            "pass2_rescue_attempts": list[dict],
            "canonical_xyz": str | None,
            "canonical_frequency_status": str,
            "canonical_frequencies_cm1": list[float],
            "sp_status": str,
            "sp_energy_hartree": float | None,
            "thermochemistry": dict | None,
            "ml_usability": dict[str, bool] | None,
            "usable_for_ml": bool | None,     # legacy aggregate
            "irc_status": str,
            "irc_endpoints": dict | None,
            "resolved_kind": str | None,
            "resolved_identity": str | None,
            "identity_status": str,
            "status": str,             # "complete" | "degraded" | "failed"
            "error": str | None,
            "attempt_history": list[dict],
            "forming_bonds": list[list[int]],
            # ... provenance / mapping / S1 thermochemistry fields
        },
    ],
    "stale_outputs_archived_to": str | None,
    "scheduling": dict,
    "summary": dict[str, int],
    "provenance": dict,
}
```

Read adapter in `manifest_io.py` accepts legacy `s3_low_level_v3` and
`s4_high_level_v4` manifests and normalizes them to the v1 shape.

## FidelityProfile

All S3/S4 differences are captured as frozen dataclass fields on
`FidelityProfile`, constructed via `FidelityProfile.from_config(config, stage)`.

Key parameter groups:

| Group | Examples |
|-------|----------|
| Methods / bases | `geometry_method: "B97-3c" / "M062X"`, `sp_method: "r2SCAN-3c" / "wB97M-V"` |
| Routes | `route_minimum: "Opt"`, `route_ts: "OptTS"` |
| Grid / SCF | `geometry_grid: None / "DefGrid3"`, `geometry_scf: None / "TightSCF"` |
| Max cycles | `max_cycles_minimum: 100 / 200`, `max_cycles_ts: 150 / 200` |
| Initial Hessian | `initial_hessian_precursor: "model"`, `initial_hessian_ts: "calculate"` |
| Rescue flags | `ts_rescue_enabled: true`, `int_rescue_enabled: true`, `irc_enabled: true` |
| Resources | `cores_per_worker`, `memory_gb_per_worker`, `max_workers`, `timeout_seconds` |
| Thermo | `temperature_k`, `standard_state`, `qrrho`, `ensemble_correction_source` |

`initial_hessian_for_role(role)` returns the policy for a given chemical role.
This is the single source of truth for the per-role initial Hessian strategy.

## KEY DESIGN DECISIONS

- **No `if stage == "S3"/"S4"` in engine code** — all differences come from
  `FidelityProfile` parameters.
- **Per-role initial Hessian**: precursor/product → `model`; intermediate/TS →
  `calculate` (Calc_Hess true per plan 7.2).
- **Warmup is role-based, not seed-state-based**: INT and TS always warm up
  with LooseOpt + bond constraints. Precursor and product never warm up.
- **Frequency is always independent**: never `UseHess`; fresh Freq on final
  geometry per plan 7.2.4.
- **INT detection uses joint criteria**: topology + mapped RMSD + normalized
  reaction progress — the legacy 1.7 A threshold is abolished per plan 8.2.
- **IRC is for TS endpoint discovery, not INT rescue**: only runs when TS
  is validated as `valid_target_ts` per plan 7.3.3.
- **Property-level ml_usability**: 9-key dict (geometry, electronic_energy,
  enthalpy, gibbs_free_energy, frequency_descriptors, ts_descriptors,
  intermediate_descriptors, mechanism_label, multifidelity_pair) replaces the
  legacy single `usable_for_ml: bool`.
- **Canonical selection uses multi-criteria sort**: stationary_rank →
  mode_rank → hessian_index → converged → alignment_score → gradient_norm.
- **Heartbeat + subprocess tracking**: optional progress reporter integration
  for live status dashboards (`rph_watch`).

## STRUCTURE REQUEST FIELDS

`StructureRequest` drives a single refinement job:

- `id` — unique structure identifier (e.g. `product_major_ts`)
- `role` — chemical role: `precursor`, `product`, `intermediate`, `ts`
- `kind` — optimization kind: `minimum` or `ts`
- `input_xyz` — path to starting geometry
- `forming_bonds` — list of `(atom_i, atom_j)` 0-based pairs
- `fallback_xyz` / `original_seed_xyz` — fallback geometry chain
- `source_stage` — `"S2"` for S3, `"S3"` for S4
- `charge` / `multiplicity` — molecular electronic state
- `atom_mapping` — optional path to atom mapping JSON
- `mapping_required` — if true, preflight validates atom_mapping exists
- `S1 thermochemistry fields` — `s1_manifest`, `s1_ensemble_thermodynamics`,
  `s1_thermochemistry_status`, `ensemble_thermochemistry_correction_hartree`
- `Provenance fields` — `structure_id`, `variant_id`, `branch_id`,
  `pathway_id`, `parent_structure_id`, `parent_structure`, `mapping_audit`

## BACKWARD COMPAT

`step3_lowlevel/__init__.py` and `step4_highlevel/__init__.py` export empty
subclasses of `RefinementEngine` named `LowLevelEngine` and `HighLevelEngine`
respectively. They preserve the legacy constructor signature
`(config, run_id, parent_manifest_paths)` and instantiate the appropriate
`FidelityProfile` internally. The orchestrator uses `RefinementEngine`
directly; the aliases exist for downstream consumers (tests, RPH_Postprocess).

Manifest read adapter in `manifest_io.py` handles two legacy schemas:

- `s3_low_level_v3` — adapts to `stage="S3"`, `fidelity="low"`,
  `profile_id="b97_3c_r2scan_3c_v1_legacy"`, converts dict-keyed structures to
  list
- `s4_high_level_v4` — adapts to `stage="S4"`, `fidelity="high"`,
  `profile_id="m062x_wb97mv_v1_legacy"`

## ANTI-PATTERNS

- Do NOT reintroduce `stage_calculator.py` or the
  `_INT_SLID_TO_PRODUCT_MAX_BOND_A = 1.7` constant.
- Do NOT call `qc_jobs.run_*` from outside `RefinementEngine` — all QC flows
  through it.
- Do NOT add `if stage == "S3"` / `if stage == "S4"` branches inside
  `engine.py` — use `FidelityProfile` parameters instead.
- Do NOT use `UseHess` to skip the independent Freq per plan 8.3.
- Do NOT run IRC when TS is not validated per plan 7.3.3.
- Do NOT hardcode theory methods, bases, or routes — put them in
  `FidelityProfile` (which reads from `config/defaults.yaml`).
