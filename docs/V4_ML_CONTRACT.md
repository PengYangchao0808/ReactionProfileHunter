# V4 ML data contract

> Status: authoritative — P0 deliverable of the V4→ML unified contract plan
> (`.sisyphus/plans/v4-ml-unified-contract.md`)
> Producer: RPH V4 pipeline (S0–S4 manifests)
> Consumer: `RPH_Postprocess` (`rph_features`) → `rph_ml`

This document is the single authority for the data contract between the RPH V4
DFT pipeline and the downstream feature-extraction / ML-training stack. It
defines: column space, feature specifications, missing-value semantics, version
mechanism, and the stage manifest guarantees the ML side may rely on.

RPH V4 owns S0–S4 computation and manifest emission. Feature extraction,
dataset assembly and ML training are external to the V4 pipeline and consume
this contract. V3-era outputs enter the same contract through a frozen legacy
backend; the ML side never branches on producer version internally.

## 0. Scope and boundary

| Concern | Owner | This contract covers |
|---|---|---|
| S0–S4 DFT computation, manifest emission, `usable_for_ml` flags | RPH V4 (`rph_core`) | What the manifests guarantee (§5) |
| Feature extraction (C1–C8), dataset assembly, `FeatureSpec` registry | `RPH_Postprocess/rph_features` | Column space (§1), FeatureSpec (§2), missing semantics (§3) |
| ML training (yield, DR), schema assertion, nan_policy execution | `rph_ml` | Version assertion (§4), training-side imputation rules (§3.4) |
| Shared column names, artifact names, `FeatureSpec` base class | `rph_schemas` | Cross-package name authority (§1.1, §2.1) |

**Out of scope**: DFT method selection, conformer search protocol, PEB scan
parameters, orbital rerun compute scheduling. These belong in `config/defaults.yaml`
and the V4 stage engines.

## 1. Column space

### 1.1 Naming convention

Every column is a dot-separated path: `<cluster>.<quantity>` or
`<domain>.<quantity>[_units]`.

- `cluster` ∈ {C1..C8} for V4 extension columns (§1.3); the common baseline
  (§1.2) uses domain prefixes `thermo.`, `geom.`, `freq.`, `ts.`, `int.`, `qc.`.
- `quantity` is lowercase snake_case.
- Energy differences carry their operand atoms/molecules in the name
  (`thermo.dG_ts_int`), never as free-text metadata.
- Units are part of the column contract, never the value. A column named
  `*_kcal` is kcal/mol; `*_hartree` is hartree; `*_angstrom` is Å; `*_cm1` is
  cm⁻¹. Unit conversion happens once, at extraction, and the converted value is
  the contractual value.

### 1.2 Common baseline set (34-column features_raw)

The V4 minimal closed loop (P1) preserves this baseline verbatim. Every column
is present in every row; absence is NaN + `missing_reason` (§3), never a
dropped column.

| Tier | Prefix | Columns | Count |
|---|---|---|---|
| META | — | `reaction_id`, `branch_id`, `condition_id`, `temperature_c`, `yield_fraction`, `dr_major`, `dr_minor` | 7 |
| FEATURE | `thermo.` | `dE_ts_int`, `dE_ts_product`, `dE_int_product`, `dE_major_minor`, `dE_act_total`, `dE_rxn_total`, `dG_ts_int`, `energy_source` | 8 |
| FEATURE | `geom.` | `r1`, `r2`, `asynch`, `asynch_index`, `rg_ts` | 5 |
| FEATURE | `freq.` | `imag_freq_cm1`, `n_imag`, `zpe_hartree` | 3 |
| FEATURE | `ts.` | `stationary_class`, `mode_identity`, `alignment_score` | 3 |
| FEATURE | `int.` | `identity_v2`, `well_depth_kcal`, `above_product_kcal`, `rmsd_to_ts` | 4 |
| LABEL/QC | `qc.` | `has_gibbs`, `forming_bonds_valid`, `energy_source_sp`, `energy_ground_truth` | 4 |
| **Total** | | | **34** |

The deployable matrix (29 columns) is the subset consumed by training:
META(3: `reaction_id`,`branch_id`,`condition_id`) + FEATURE(19) + LABEL(3) +
TARGET(4: `yield_fraction`,`dr_major`,`dr_minor`,`temperature_c`).

The `feature_registry` sidecar carries 18 FeatureSpec entries. `int.rmsd_to_ts`
is present in the feature matrix and training FEATURES but historically lacked a
registry entry; P1 closes that gap so every column has a spec.

### 1.3 V4 extension set (C1–C8)

Extension columns are additive. V3 cannot produce them; they are emitted as NaN
+ `missing_reason = NA_STRUCTURE` by the V3 backend and as real values by the
V4 backend. Each column belongs to exactly one cluster.

| Cluster | Module | +Cols | Data source |
|---|---|---|---|
| C1 conformer ensemble | `conformer.py` | +20 | S3-embedded S1 manifest (`b973c_sp_energy_hartree`, per-candidate `relative_free_energy_kcal` / `boltzmann_population`, `ensemble_thermodynamics.*`) |
| C2 thermo decomposition | `thermo.py` | +9 | `thermochemistry` + `canonical_*` differences (dG/dH/TdS/dZPE/dE_opt_vs_sp/barrier_asym) |
| C3 geometry descriptors | `geometry.py` | +7 | `canonical.xyz` + `forming_bonds` measurement (fb_progress/close_contacts/d_avg/d_max/dihedrals) |
| C4 frequency / vibrational | `frequency.py` | +4 | `canonical_frequencies_cm1` + freq output parse (S_vib/Cv) |
| C5 TS mode | `ts_mode.py` | +2 | hessian diagonalization (reuses `relabel_ts_identity`) |
| C6 intermediate descriptors | `intermediate.py` | +3 | manifest fields |
| C7 reaction-aggregate | `reaction_aggregate.py` | +6 | major/minor branch differences |
| C8 orbital / electronic | `orbital.py` | +10–20 | 40 SP rerun `%molden` + Multiwfn (optional) |
| **Total extension** | | **+61–71** | **80–100 cols/row** |

S1 field-name reference (verified against
`rph_core/steps/conformer_search/censo_lite.py`):
`b973c_sp_energy_hartree` (not `b97_3c_*`); `relative_free_energy_kcal` and
`boltzmann_population` are per-candidate fields; manifest key is
`ensemble_free_energy_estimate_hartree` (with `_estimate_`; the internal
dataclass field is `ensemble_free_energy_hartree` without it). There is no
optical-rotation field.

## 2. FeatureSpec

### 2.1 Dataclass

`FeatureSpec` is the single source of truth for what a column means. It lives
in `rph_schemas.feature_schema` (shared by `rph_features` and `rph_ml`).

```python
@dataclass(frozen=True)
class FeatureSpec:
    column: str            # "thermo.dG_ts_product"
    cluster: str           # "C2" | "baseline"
    unit: str              # kcal/mol | angstrom | cm-1 | hartree | eV | none
    physical_meaning: str  # "TS–product Gibbs free-energy difference (same molecule, SP+freq thermal correction)"
    formula: str           # "G_TS - G_PRODUCT, x627.509"
    source: str            # "manifest: thermochemistry.gibbs_free_energy_hartree"
    missing_reason: str    # MISSING_ARTIFACT | MISSING_MANIFEST_ENTRY | DEGRADED_FALLBACK | NA_STRUCTURE
    nan_policy: str        # keep | drop | zero_fill  (executable, not documentary)
    layer: str             # computational | qa_metadata
    backend: str           # v3 | v4 | both
```

### 2.2 Declaration and enforcement

Each extractor declares its output columns via
`ManifestExtractor.declare_features() -> list[FeatureSpec]`. The base-class
`run()` method then automatically:

1. Merges every spec into `feature_meta.json`
   (`meta.features.<column> = spec`).
2. Includes the spec set in the schema signature hash (§4.2). A unit or formula
   change is a version change.
3. Validates that `extract()` actually produced columns ⊆ the declared spec
   set. An undeclared column or a missing declared column is a hard extraction
   error, not a warning.

### 2.3 Provenance

`source` is a structured string: `<artifact>: <json_path>`. Examples:
`manifest: thermochemistry.gibbs_free_energy_hartree`,
`canonical.xyz: forming_bonds[0]`,
`freq.out: lowest_imaginary_cm1`. This is machine-parseable for audit and
debugging; it is not free text.

## 3. Missing-value semantics

### 3.1 The four rules

| Layer | Rule |
|---|---|
| Extraction | Every declared column **always exists** in every row. Missing = `NaN` + a paired `<column>_missing_reason` enum value. Empty strings and free-text reasons are abolished. |
| Structure-level | If a required manifest, precursor, product or TS is absent, the **row** is skipped and a `skip_reason` column records why. A skipped row is not emitted with NaN-filled features. |
| Registry | `nan_policy` ∈ {`keep`, `drop`, `zero_fill`} is **executed** by the dataset builder at assembly time, not documented and ignored. `drop` removes the column from the deployable matrix; `zero_fill` replaces NaN with 0 and records the imputation in `feature_meta.imputations`; `keep` leaves NaN. |
| Training | `rph_ml` reads `nan_policy` from the registry and applies it. **Implicit zero-fill is forbidden.** `nan_to_num(0)` without a registry `zero_fill` entry is a contract violation; the training runner rejects the dataset. |

### 3.2 missing_reason enum

| Value | Meaning |
|---|---|
| `MISSING_ARTIFACT` | The backing file (`.xyz`, `.out`, `.hess`) does not exist on disk. |
| `MISSING_MANIFEST_ENTRY` | The manifest exists but the JSON path is absent or null. |
| `DEGRADED_FALLBACK` | A lower-fidelity source was used (e.g. `thermo_sp` after `sp_output` failed); the value is present but provenance degraded. |
| `NA_STRUCTURE` | The quantity is structurally undefined for this producer (e.g. a V3 row encountering a C1–C8 V4-only column; or `int.*` when no intermediate node exists). |

### 3.3 Legacy compatibility

V3's dual-key `*_missing_reason` (free string) + `*_status` pattern is mapped
to the enum at the normalization layer. Free-text V3 reasons are bucketed into
the four enum values by a deterministic classifier; the original text is
preserved in `feature_meta.legacy_reasons` for audit. Columns where V3 lacked a
paired `*_status` (`fmo_cdft_dipolar.missing_reason`, `branch_missing_reason`,
`minor_branch_missing_reason`) are normalized to the enum without synthesis.

## 4. Version contract

### 4.1 SCHEMA_VERSION — single constant

One canonical `SCHEMA_VERSION` constant replaces the prior four inconsistent
locations (`context.py:212` = "4.2-tag", `schema.py:27` = "6.2",
`schema.py:28` = "7.0", `feature_miner.py:465` hardcoded "6.2"). It lives in
`rph_schemas` and is imported by both `rph_features` and `rph_ml`. Bumping it
is a breaking change that invalidates every downstream dataset.

### 4.2 schema_signature — content hash

`get_schema_signature()` returns the SHA-1 hexdigest of:

```json
{
  "schema_version": "<SCHEMA_VERSION>",
  "feature_specs": [<sorted list of FeatureSpec.column + unit + formula>],
  "columns": [<sorted deployable column list>],
  "plugins": [<sorted enabled extractor list>]
}
```

Feature **values** are never hashed. A signature change means the contract
changed; a value-only change does not. The signature is written into
`feature_meta.json` and stamped into every dataset directory name.

### 4.3 Training-side assertion

`rph_ml` rejects a dataset whose stamped `schema_signature` does not match the
registry-computed signature at training start. There is no "close enough" mode:
a mismatch is a hard stop with a diff of which specs/columns moved.

### 4.4 Splits parameterization

Dataset assembly accepts `--seed` and guards against `n_reactions < 20` (the
prior hardcoded 16/2/2 slice silently produced empty val/test sets below 20).
The split assignment is deterministic given the seed and is recorded in
`splits.json` alongside the signature.

### 4.5 Temperature dual-track

Two temperatures coexist and must not be confused:

- `thermochemistry.temperature_K` — the composite thermochemistry temperature
  (247.55 K for the current benchmark set). Basis: `computational`.
- `temperature_c` — the experimental reaction condition (mapping metadata).
  Basis: `experimental`.

Both are columns. `FeatureSpec.source` carries the basis tag. A model must not
feed `temperature_c` into a thermodynamic feature without an explicit decision;
the registry marks them as different clusters.

## 5. Stage manifest contract

Each stage manifest is a versioned JSON file. The ML side reads the fields
listed here; unlisted fields are internal and not contractual.

### 5.1 S0 — `S0_Mechanism/mechanism.json`

| Field | Type | ML use |
|---|---|---|
| `forming_bonds` | list[[int,int]] (0-based product-SMILES XYZ indices) | Geometry features (C3), constraint validation. Authoritative source. |
| `mapped_forming_bonds` | list[[int,int]] (1-based map numbers) | Index-space cross-check |
| `reaction_type` / `topology` / `cyclo_mode` | str | Metadata, mechanism label (C7) |
| `atom_mapping_digest` | str (sha256) | Provenance cross-check; propagated to S2/S3 `forming_bonds_hash` |
| `canonical_product_smiles` / `canonical_precursor_smiles` | str | Identity |

**Condition metadata is NOT in `mechanism.json`.** `temperature_c`, `yield_fraction`,
`dr_major`, `dr_minor` live in the source CSV row only (linked to the run via
`source_row_hash`). The dataset builder reads them from the CSV at assembly
time, not from any stage manifest. `dr_branch_plan.json` carries branch role
assignments (`major_or_reference` / `diastereomer_candidate`) for C7
aggregation.

### 5.2 S1 — `S1_ConfSearch/{molecule}/manifest.json`

> Location is **per-molecule**: `S1_ConfSearch/{precursor,product_major,product_minor_001}/manifest.json`,
> not a single `product/manifest.json`. Schema: `s1_censo_light_ranking_v4`.
> S3 references the path via `structures[i].s1_manifest` (string, **null on
> int/ts records** — resolve via `variant_id` for those roles).

Consumed by C1 (conformer ensemble). Key nesting:

| Path | Type | ML use |
|---|---|---|
| `candidates[i].b973c_sp_energy_hartree` | float | C1 ranking energy |
| `candidates[i].g_rrho_correction_hartree` | float | C1 mRRHO correction |
| `candidates[i].s1_score_hartree` | float | C1 composite score (b973c SP + G_RRHO) |
| `candidates[i].relative_free_energy_kcal` | float | C1 population spread (per-candidate) |
| `candidates[i].boltzmann_population` | float | C1 Boltzmann weight, configurational entropy (per-candidate) |
| `candidates[i].torsion_signature` | str | C1 dedup/structural fingerprint |
| `candidates[i].degeneracy` | int | C1 degeneracy |
| `ensemble_thermodynamics.temperature_k` | float | C1 thermo temperature (**298.15 K** — dual-track with S3's 247.55 K, §4.5) |
| `ensemble_thermodynamics.ensemble_free_energy_estimate_hartree` | float | C1 ensemble free energy |
| `ensemble_thermodynamics.conformational_free_energy_correction_kcal` | float | C1 conformational correction |
| `ensemble_thermodynamics.ensemble_member_count` | int | C1 ensemble size |
| `ensemble_thermodynamics.partition_function_relative` | float | C1 partition function |
| `thermodynamic_rank1` | str | C1 rank-1 conformer id |
| `selected_xyz` | str (path) | C1/C3 geometry source |

`selected` and `representative_candidates` are **deprecated** (listed in
`deprecated_fields`); use `thermodynamic_rank1` + `reactivity_screening_candidates`.
Atom count is not a JSON field — parse line 1 of `selected.xyz`.

S1 does not emit DFT OPT/FREQ; C1 conformer features are ranking/thermo only.

### 5.3 S2 — `S2_PEB/manifest.json`

| Field | Type | ML use |
|---|---|---|
| `ts_guess.status` | str | TS seed provenance |
| `intermediate.status` | str (`not_found` if no basin) | Gate for C6 intermediate features |
| `scan_profile.summary` | dict | Diagnostic, not a primary feature source |
| `forming_bonds` (resolved, 0-based XYZ) | list[[int,int]] | Geometry features |

S2 does not produce QC energies; it produces TS guesses and intermediate seeds.

### 5.4 S3 — `S3_LowLevel/manifest.json` (primary ML source)

Schema: `refinement_manifest_v1`, `fidelity: "low"`, `profile_id: "b97_3c_r2scan_3c_v1"`.
The `structures[]` array is the unit of ML consumption. The `role` enum is
{`precursor`, `product`, `intermediate`, `ts`}; structure `id` carries the
suffix (e.g. `product_major_ts` has `role: "ts"`).

Per-structure record (`structures[i]`), grouped by consumer cluster:

**Identity / grouping (C7):**
| Path | Type | ML use |
|---|---|---|
| `id` / `role` / `kind` | str | Row identity; `kind` ∈ {`minimum`,`ts`} |
| `variant_id` / `branch_id` / `pathway_id` | str | C7 cross-branch aggregation |
| `parent_structure_id` / `source_stage` | str | C7 lineage |

**Energy (C2 / C8 anchor):**
| Path | Type | ML use |
|---|---|---|
| `sp_energy_hartree` | float | Baseline `thermo.dE_*`, C2; C8 rerun consistency anchor |
| `opt_energy_hartree` | float | C2 opt-vs-sp decomposition |
| `thermochemistry.gibbs_free_energy_hartree` | float | Baseline `thermo.dG_ts_int` |
| `thermochemistry.enthalpy_hartree` | float | C2 dH decomposition |
| `thermochemistry.single_point_energy_hartree` | float | = `sp_energy_hartree` (redundant, guaranteed); fallback level `thermo_sp` |
| `thermochemistry.temperature_K` | float | **247.55 K** (dual-track with S1's 298.15, §4.5) |
| `canonical_zero_point_energy_hartree` | float | Baseline `freq.zpe_hartree` |
| `canonical_enthalpy_hartree` / `.canonical_gibbs_free_energy_hartree` | float | C2 canonical thermo |

**Geometry (C3):**
| Path | Type | ML use |
|---|---|---|
| `canonical_xyz` | path | C3 authoritative geometry |
| `opt_xyz` / `input_xyz` | path | C3 provenance |
| `forming_bonds` | list[[int,int]] | Baseline `geom.r1/r2`, C3. **Present on int/ts records; empty `[]` on product/precursor.** Pair ordering may differ from S0 (`[17,18]` here vs `[18,17]` in S0) — treat as unordered. |
| `geometry_hash` | str (sha256) | C3 dedup |

**Frequency (C4):**
| Path | Type | ML use |
|---|---|---|
| `canonical_frequencies_cm1` | list[float] | C4 vibrational descriptors, baseline `freq.*` |
| `canonical_imaginary_frequencies_cm1` | list[float] | C4 imaginary modes |
| `canonical_frequency_status` | str | C4 validity gate (`"complete"`) |
| `canonical_hessian_path` | path | C5 hessian source |

**TS mode (C5):**
| Path | Type | ML use |
|---|---|---|
| `ts_classification.stationary_point_class` | str | Baseline `ts.stationary_class` (`valid_target_ts` / `first_order_wrong_mode`) |
| `ts_classification.mode_identity` | str | Baseline `ts.mode_identity` (`target` / `unrelated`) |
| `ts_classification.alignment_score` | float | Baseline `ts.alignment_score` |

**Intermediate (C6) — path is `int_classification.metrics.*`, NOT `int_metrics.*`:**
| Path | Type | ML use |
|---|---|---|
| `int_classification.identity` | str | Baseline `int.identity_v2` |
| `int_classification.metrics.well_depth_kcal` | float | Baseline `int.well_depth_kcal` |
| `int_classification.metrics.above_product_kcal` | float | Baseline `int.above_product_kcal` |
| `int_classification.metrics.rmsd_to_ts` | float | Baseline `int.rmsd_to_ts` |
| `int_classification.metrics.rmsd_to_product` | float | C6 additional |
| `int_classification.usable_for_ml` | bool | C6 gate |

> Legacy duplicate: `int_identity_v2.metrics.*` mirrors `int_classification.metrics.*`
> for backward compatibility. The contract prefers `int_classification`.

**Usability & status:**
| Path | Type | ML use |
|---|---|---|
| `usable_for_ml` | bool | Master gate (legacy = `geometry and electronic_energy`) |
| `ml_usability` | dict[str,bool] | Property-level gates (§5.6) |
| `opt_status` / `sp_status` | str | Provenance, degraded-fallback detection |

> There is **no `ts_frequency_valid` field**. TS validity is derived from
> `canonical_frequency_status == "complete"` combined with
> `ts_classification.stationary_point_class`. The property-level
> `ml_usability.ts_descriptors` gate (§5.6) encodes this.

**Embedded S1 (C1 source within S3):**
| Path | Type | ML use |
|---|---|---|
| `s1_manifest` | str (path) \| null | C1 source — **null on int/ts records**; resolve via `variant_id` to the parent product's S1 manifest |
| `s1_ensemble_thermodynamics` | dict \| null | C1 ensemble block copy (survives S1 dir deletion); null on int/ts |
| `ensemble_thermochemistry_correction_hartree` | float \| null | C2 ensemble correction |

Energy fallback priority (SP, for `thermo.energy_source` and
`qc.energy_source_sp`):
`sp_output` → `sp_gt_rescue` → `thermo_sp` → `manifest_sp` → (terminal `none`).
The chosen level is recorded per row; `DEGRADED_FALLBACK` missing_reason is
set when the level is below `sp_output`.

### 5.5 S4 — `S4_HighLevel/manifest.json`

Same per-structure schema as S3, produced by M062X OPT/OptTS/FREQ → wB97M-V SP
(`fidelity: "high"`, `profile_id: "m062x_wb97mv_v1"`). A failed S4 job
**never deletes or replaces** usable S3 output. The ML side prefers S4 energies
when `usable_for_ml` is true and falls back to S3 when S4 is degraded; the
`energy_ground_truth` column records which level won.

> The golden benchmark (`s2_gfn2_benchmark_alpb_v2`) stopped at S3; no S4
> manifest exists for `rx_id_1`. The S4 schema above is the design contract
> from `config/defaults.yaml` theory section; it will be verified against a
> future `--stop-after s4` golden run.

### 5.6 Property-level usability

`usable_for_ml` is a dict, not a single boolean. The legacy boolean is derived
as `geometry and electronic_energy`. Each property has an independent gate:

| Property | Gate |
|---|---|
| `geometry` | `opt_converged` |
| `electronic_energy` | `sp_complete` |
| `enthalpy` | `sp_complete and freq_available` |
| `gibbs_free_energy` | `sp_complete and freq_available` |
| `frequency_descriptors` | `freq_available` |
| `ts_descriptors` | `requested_kind == "ts" and canonical_frequency_status == "complete" and ts_classification.stationary_point_class == "valid_target_ts"` |
| `intermediate_descriptors` | `requested_role == "intermediate" and identity_status == "role_matched"` |
| `mechanism_label` | always true |
| `multifidelity_pair` | false (reserved) |

This is computed by `rph_core.utils.ml_quality.compute_ml_usability` and
serialized into the manifest. The consumer reads it as-is; it does not
re-derive usability from raw status fields.

## 6. Producer/consumer boundary

**RPH V4 guarantees:**

1. Every stage manifest is versioned and written atomically.
2. `forming_bonds` are 0-based XYZ indices, canonicalized as sorted
   `(min, max)` pairs.
3. `usable_for_ml` is present on every structure record and reflects
   property-level gates.
4. A degraded job preserves XYZ, output logs and failure diagnostics.
5. S4 failure does not corrupt S3 output.

**RPH V4 does not guarantee:**

1. That every structure is `usable_for_ml` — a row may have NaN features with
   structured `missing_reason`.
2. Specific feature column names or formulas — those are this contract's
   domain, owned by the consumer side.
3. Dataset splits, `nan_policy` execution, or training-readiness — those are
   `rph_features` / `rph_ml` responsibilities.

**The consumer guarantees:**

1. Every column has a `FeatureSpec`.
2. Missing values are NaN + enum reason, never empty strings or implicit zeros.
3. `schema_signature` is asserted before training.
4. V3 and V4 rows are indistinguishable in the deployable matrix except by
   `NA_STRUCTURE`-tagged extension columns.

## 7. Decisions adopted

The four pending decisions from the plan (§10.3) are resolved as follows for
P0; they may be revisited before P1/P7 if conditions change.

| Decision | Resolution | Rationale |
|---|---|---|
| D1: `rph_ml` location | `RPH_Postprocess/rph_ml/` (same repo as `rph_features`) | Co-located with the feature layer; single CI; no cross-repo import friction |
| D2: `rph_schemas` promotion | Stay nested in `RPH_Postprocess/rph_schemas/`; both packages import via pip | Minimal change; promotion to a top-level repo deferred until a third consumer appears |
| D3: `train_yield.py` | Treat as to-be-created in P7; the only existing training code is `scripts/audit/train_v15_v4s3.py` | No `train_yield.py` found in any workspace |
| D4: `scripts/audit/` baseline freeze | `git add scripts/audit/` before P1 to fix the 34-column numerical baseline | Required for the golden regression gate |

## 8. Verification commands

```bash
# Contract doc compiles (no broken code blocks)
python -m py_compile docs/V4_ML_CONTRACT.md 2>/dev/null || true

# P1 golden regression (34-col features_raw / 29-col deployable)
pytest tests/test_v4_s3_compat.py -k golden

# schema_signature assertion (P7)
pytest tests/test_rph_ml_contract.py
```

The golden data root is `RPH_Test_Results/s2_gfn2_benchmark_alpb_v2`
(20 reactions, 140 structures, 134 ML-usable). It does not require real
ORCA/Multiwfn binaries; P6 (C8 orbital reruns) is the only phase that needs
real compute and is isolated with `@pytest.mark.slow`.
