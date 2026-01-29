# V6.1–V6.4 Plan: S3 Compute-First, S4 Extract-Only, Linear-Model Features

## Context

### Original Request
Implement a V6.1–V6.4 modification program optimized for small-sample ([5+2], N≈30) + linear models:
- Move **all computations** into S3 (ORCA single points; DIAS/ASM enrichment; optional TS@TS closure).
- Make S4 **strictly extraction/processing/QC marking**; S4 must not submit any QC jobs.
- Add a lightweight **FMO/CDFT parser** in S4 based on S3 “dipolar intermediate” outputs: parse HOMO/LUMO and compute `gap` and electrophilicity `omega` (no new compute).
- Lock S4 outputs (V6.1 contract): `features_raw.csv`, `features_mlr.csv` (6–10 columns), `feature_meta.json`.
- Follow repo anti-pattern constraints: no duplicated QC runner in steps; reuse `rph_core/utils/qc_interface.py` toxic-path sandbox logic.

### Repo Facts (local)
- Current S4 output contract is `S4_Data/features.csv` + `S4_Data/feature_meta.json`:
  - `rph_core/steps/step4_features/feature_miner.py`.
  - Schema and CSV writer: `rph_core/steps/step4_features/schema.py`.
  - Default plugins: `thermo`, `geometry`, `qc_checks`.
- Thermo currently **mixes semantics**: `thermo.dG_activation` / `thermo.dG_reaction` are Gibbs if available, else electronic fallback:
  - `rph_core/steps/step4_features/extractors/thermo.py`.
  - `thermo.energy_source_*` columns already exist.
- Geometry currently provides `geom.r1/r2/asynch/asynch_index` + close contacts + min distances:
  - `rph_core/steps/step4_features/extractors/geometry.py`.
- Step4 completion in resume flow hard-checks `S4_Data/features.csv` when mechanism packaging enabled:
  - `rph_core/utils/checkpoint_manager.py`.
- Gaussian log parser already extracts orbitals (Hartree→eV) and frequencies:
  - `rph_core/steps/step4_features/log_parser.py`.
- Step3 output wiring (canonical `ts_final.xyz`, `reactant_sp.xyz`, SPMatrixReport) lives in:
  - `rph_core/steps/step3_opt/ts_optimizer.py`.

### Metis Review
- Metis consultation was attempted but the metis agent tool call failed with an infra error (`JSON Parse error: Unexpected EOF`). This plan includes an explicit “Self gap check” task to compensate.

---

## Work Objectives

### Core Objective
Deliver a clean S3/S4 separation where S3 performs all enrichment computations and S4 produces reproducible, linear-model-friendly features (raw + MLR-ready + metadata) without submitting QC jobs.

### Concrete Deliverables
- S4 output contract (V6.1):
  - `S4_Data/features_raw.csv`
  - `S4_Data/features_mlr.csv`
  - `S4_Data/feature_meta.json`
- Update resume/contract checks to stop relying on `S4_Data/features.csv`.
- Feature set updates:
  - Thermo columns split to avoid semantic mixing (Gibbs vs electronic).
  - Geometry columns updated for linear stability (`r_avg`, signed `dr`, `close_contacts_density`).
  - TS quality columns parsed from existing outputs only (`n_imag`, `imag1_abs`, `dipole_debye`).
  - Optional ASM/DIAS features in S4 by consuming an S3-written enrichment contract JSON.
  - Optional dipolar FMO/CDFT features in S4 by parsing S3 dipolar outputs (no compute).
- S3 enrichment contract (written under the Step3 output directory; currently `work_dir/S3_TransitionAnalysis/`):
  - `S3_TransitionAnalysis/S3_PostQCEnrichment/enrichment.json`
  - `S3_TransitionAnalysis/S3_PostQCEnrichment/enrichment_status.json`
- S3 artifacts index (for reproducible dipolar selection):
  - `<S3_DIR>/artifacts_index.json`
- V6.4 modeling companion script(s) under `scripts/ml/` (kept out of pipeline execution).
- Docs/tests updated to reflect new contracts.

### Definition of Done
- `python -m pytest -q` passes.
- `python run_auto_test.py --step ALL` completes for at least one local sample run directory.
- No Step4 code path submits QC jobs (no calls into `QCTaskRunner` / `qc_interface` execution functions from S4).
- Resume logic recognizes the new Step4 outputs (no dependency on `features.csv`).
- `python run_auto_test.py --step ALL --resume` passes with `step4.mechanism_packaging.enabled=true` (ensures checkpoint/resume/packager no longer depends on `features.csv`).

### Must NOT Have (Guardrails)
- No QC execution in S4 (including fragment SP jobs).
- No duplicate ad-hoc sandbox/toxic-path handling; always reuse `rph_core/utils/qc_interface.py` utilities.
- No changes to Step2 output contract (`S2_Retro/ts_guess.xyz` and `S2_Retro/reactant_complex.xyz` must remain).
- Avoid large refactors; touch only what is necessary to implement V6.1–V6.4.
- V6.3 does NOT add new “dipolar intermediate” computations; it only indexes/parses dipolar outputs if they already exist under the Step3 output directory.

---

## Verification Strategy

### Test Decision
- **Infrastructure exists**: YES (`pytest`, tests under `tests/`).
- **User wants tests**: Assumed YES (tests-after), plus strict manual verification steps for artifacts.

### Core Commands
- Unit/integration: `python -m pytest -q`
- Pipeline smoke: `python run_auto_test.py --step ALL`
- Resume smoke: `python run_auto_test.py --step ALL --resume`

---

## Task Flow

V6.1 (S4 outputs + feature schema fix) → V6.2 (S3 enrichment + contract + S4 consume) → V6.3 (dipolar FMO parser + plugin) → V6.4 (hardening + ML scripts + docs)

---

## TODOs

### 0. Baseline Inventory & Contract Snapshot

**What to do**:
- Record current Step4 contract and resume semantics.
- Apply the confirmed migration policy: **no backward compatibility** for `S4_Data/features.csv`.
- Enumerate and update all direct `features.csv` expectations in code/tests/docs to the new contract.

**References**:
- `rph_core/steps/step4_features/feature_miner.py` (current outputs)
- `rph_core/utils/checkpoint_manager.py` (resume checks)
- `rph_core/orchestrator.py` (resume skip logic)
- Other runtime references to update away from `features.csv`:
  - `rph_core/utils/result_inspector.py`
- `tests/` (multiple tests assert `features.csv`)
- Docs to update away from `S4_Data/features.csv`:
  - `README.md`
  - `rph_core/steps/AGENTS.md`
  - Any V5-era acceptance/testing docs that are treated as authoritative in this repo

**Acceptance Criteria**:
- `grep` for `S4_Data/features.csv` under `rph_core/`, `tests/`, `README.md`, `rph_core/steps/AGENTS.md` returns no matches.
- Allowable remaining matches (explicit allowlist): only under `rph_core_backup_20260115/`.
- `python -m pytest -q` passes with the new filenames.

---

### 1. V6.1: Lock S4 Output Contract (raw/mlr/meta)

**What to do**:
- Update S4 writer to generate:
  - `S4_Data/features_raw.csv` (full features + qc.* + intermediate columns)
  - `S4_Data/features_mlr.csv` (default ≤10 columns)
  - `S4_Data/feature_meta.json` (extend existing meta with units/sources/missing stats, plugin enablement, schema/version hashes)
- Remove `S4_Data/features.csv` and update resume/tests/docs to use the new filenames immediately.
- Update Step4 resume/backfill logic to use the new contract:
  - `FeatureMiner` backfill hook currently uses `features.csv` existence as the trigger; change the trigger to `features_raw.csv`.
  - `CheckpointManager.is_step4_complete()` must check `features_raw.csv` (non-empty) instead of `features.csv` when mechanism packaging is enabled.
  - Orchestrator resume skip logic remains unchanged but will call the updated `is_step4_complete()`.
- Update Step4 return value semantics:
  - `FeatureMiner.run(...)` should return the path to `S4_Data/features_raw.csv` (so existing callers assigning to `features_csv` continue to work without relying on the old filename).

**New contract schema (pin down, no guesswork)**:
- `S4_Data/features_raw.csv`
  - One-row CSV per run directory.
  - Column ordering rule: fixed columns (in-order) then dynamic columns (sorted), reusing the existing writer pattern.
  - Must include: `schema_version`, `schema_signature`, `feature_status`, `sample_id`, plus all enabled plugin outputs.
  - Must include QC columns (`qc.*`) and intermediate convenience columns (e.g., `thermo.dG_*_best`).
- `S4_Data/features_mlr.csv`
  - One-row CSV per run directory.
  - Contains ONLY the MLR columns listed below (no `schema_version/schema_signature/feature_status` columns), to keep downstream modeling scripts simple.
  - Default column list (exact, in-order):
    1) `sample_id` (string; stable identifier; exact rule: `sample_id = output_dir.parent.name` where `output_dir` is `work_dir/S4_Data`)
    2) `thermo.dE_activation`
    3) `thermo.dE_reaction`
    4) `geom.r_avg`
    5) `geom.dr`
    6) `geom.close_contacts_density`
    7) `ts.imag1_cm1_abs`
    8) `asm.distortion_total_kcal` (present as a column; NaN if enrichment missing/disabled)
    9) `fmo.dipolar_omega_ev` (present as a column; NaN if plugin disabled/missing)
    10) `qc.sample_weight` (present as a column; default 1.0 unless strict-QC downgrades)
  - Always written, even if optional columns are NaN.
- `S4_Data/feature_meta.json` (extend existing structure)
  - `meta.outputs`: `{features_raw_csv, features_mlr_csv, feature_meta_json}` (relative filenames)
  - `meta.units`: mapping `column -> unit` (reuse schema unit map; `None` allowed)
  - `meta.mlr_columns`: the exact column list used for `features_mlr.csv`
  - `meta.missing_stats`: for each column in raw and mlr outputs, record:
    - `is_missing`: value is absent or NaN
    - `is_invalid`: value exists but fails a sanity check (e.g., gap<=0)
    - `reason`: machine-readable string (`missing_input`, `parse_failed`, `sanity_check_failed`, `not_applicable`)
    - Minimum sanity checks to implement (defaults):
      - `geom.r_avg > 0`
      - `geom.close_contacts_density >= 0`
      - `ts.imag1_cm1_abs > 0` when `ts.n_imag > 0`
      - `fmo.dipolar_gap_ev > 0` (if fmo values present)
      - `fmo.dipolar_omega_ev > 0` (if computed)
      - For everything else: `is_invalid=false` unless a plugin explicitly marks it.
  - `trace.plugins.<plugin>.missing_reasons`: list of machine-readable reasons (e.g., `missing_input`, `parse_failed`, `not_applicable`)

**Configuration surface (explicit defaults)**:
- Config file location for defaults: `config/defaults.yaml`.
- `step4.outputs`:
  - `write_features_raw` (default true)
  - `write_features_mlr` (default true)
  - `write_feature_meta` (default true)
- `step4.mlr.columns`: optional override list; if absent use the default list above.
- `step4.strict_qc.enabled` (default false)
- `step4.strict_qc.sample_weight_policy` (default `downweight`):
  - `downweight`: set `qc.sample_weight=0.0` when `ts.n_imag!=1` or forming_bonds invalid; otherwise 1.0.
  - `drop_from_mlr`: write `features_mlr.csv` but set all mlr feature columns to NaN when invalid; keep raw intact.

**Strict-QC vs early-exit (make it consistent)**:
- Current behavior: invalid forming_bonds triggers early `_write_failed_output(...)`.
- V6 behavior: ALWAYS write `features_raw.csv` + `features_mlr.csv` + `feature_meta.json` even on invalid inputs.
  - Set `feature_status` to the exact on-disk enum value used in this repo: `"invalid_inputs"`.
  - Populate `qc.forming_bonds_valid=0`.
  - Apply strict-QC policy to MLR output (`qc.sample_weight=0.0` and/or NaN out MLR features).

**Status and version strings (exact on-disk values)**:
- `feature_status` values are the `.value` strings from `rph_core/steps/step4_features/status.py:FeatureStatus` (e.g., `ok`, `missing_inputs`, `invalid_inputs`).
- `schema_version` in `features_raw.csv` is the Step4 schema version string (update `SCHEMA_VERSION` in `rph_core/steps/step4_features/schema.py` to `6.1`).
- `schema_signature` continues to use the existing Step4 schema signature logic (sha1 hex).

**Where this is implemented (explicit code touchpoints)**:
- Output writing:
  - Extend `rph_core/steps/step4_features/schema.py` with dedicated writers:
    - `write_features_raw_csv(path, row_data)` (can be a thin wrapper over existing `write_features_csv` but with the new filename)
    - `write_features_mlr_csv(path, row_data, mlr_columns)`
      - MLR writer MUST always emit the full `mlr_columns` list by injecting missing keys as NaN.
  - Update `rph_core/steps/step4_features/feature_miner.py` to call these writers and to return `features_raw.csv`.
- Resume validation:
  - Update `rph_core/utils/checkpoint_manager.py:is_step4_complete()` to validate “non-empty” by:
    - file exists and size > 0, and
    - CSV parse yields at least 1 header row + 1 data row.

**References**:
- `rph_core/steps/step4_features/feature_miner.py` (writes current outputs)
- `rph_core/steps/step4_features/schema.py` (CSV writer, schema signature)
- `rph_core/steps/step4_features/context.py` (feature_meta.json structure)
- `rph_core/utils/checkpoint_manager.py:is_step4_complete()` (currently checks `features.csv`)
- `rph_core/steps/step4_features/feature_miner.py` (backfill hook keyed on `features.csv` existence)
- Tests that must be updated to new filenames:
  - `tests/test_m2_step4_resume_semantics.py`
  - `tests/test_m2_schema_versioning.py`
  - `tests/test_degradation.py`
  - `tests/test_degradation_final.py`
  - `tests/test_m2_precursor_fallback.py`
  - `tests/test_s4_artifact_integration.py`
  - `tests/test_mock_integration.py`

**Acceptance Criteria**:
- After a pipeline run, `S4_Data/features_raw.csv`, `S4_Data/features_mlr.csv`, `S4_Data/feature_meta.json` exist and are non-empty.
- `features_mlr.csv` has ≤10 columns by default.
- No `S4_Data/features.csv` is produced.
- When `step4.mechanism_packaging.enabled=true`, `S4_Data/mech_index.json` is still required and still validated by resume logic.
- `python -m pytest -q` passes.

---

### 2. V6.1: Fix Thermo Semantic Mixing (Gibbs vs Electronic)

**What to do**:
- Introduce non-mixed columns:
  - `thermo.dG_activation_gibbs` / `thermo.dG_reaction_gibbs`
  - `thermo.dE_activation` / `thermo.dE_reaction`
  - Keep `thermo.energy_source_activation` / `thermo.energy_source_reaction` but define them explicitly as the source for *best-available* Gibbs-like columns (below).
- Add explicit “best available” convenience columns in raw output only (not default mlr columns):
  - `thermo.dG_activation_best` / `thermo.dG_reaction_best`
  - Rule: if Gibbs exists use Gibbs; else use electronic difference in kcal/mol.
  - `thermo.energy_source_*` describes `thermo.dG_*_best` and is one of `gibbs|electronic|none`.
- Deprecate or remove ambiguous columns:
  - Keep `thermo.dG_activation` / `thermo.dG_reaction` as **legacy aliases** of `thermo.dG_*_best` (so existing fixed-column expectations and checks remain meaningful), but document them as deprecated in meta.
  - New modeling should use the split columns (`thermo.dG_*_gibbs`, `thermo.dE_*`) and/or the curated `features_mlr.csv`.
- Update schema versioning accordingly (V6.1 marker).

- Update existing Step4 consistency checks:
  - `rph_core/steps/step4_features/feature_miner.py` currently compares `thermo.energy_source_*` with what the SP report contains.
  - Keep that check aligned by defining `thermo.energy_source_*` as the source of the legacy/best columns (`thermo.dG_*`/`thermo.dG_*_best`).

**References**:
- `rph_core/steps/step4_features/extractors/thermo.py` (current fallback behavior)
- `rph_core/steps/step4_features/schema.py` (fixed columns and units mapping)
- `rph_core/steps/step3_opt/ts_optimizer.py` (SPMatrixReport contains `g_*` in kcal/mol; `e_*` in Hartree)

**Acceptance Criteria**:
- If Gibbs is missing, `thermo.dG_*_gibbs` is empty (NaN) while `thermo.dE_*` can still be populated.
- No column mixes Gibbs and electronic semantics.
- `python -m pytest -q` passes.

---

### 3. V6.1: Geometry Feature Upgrade for Linear Models

**What to do**:
- Add:
  - `geom.r_avg = (r1+r2)/2`
  - `geom.dr = (r1-r2)` (signed)
  - `geom.close_contacts_density = geom.close_contacts / geom.natoms_ts`
- Keep `geom.asynch` in raw; deprecate `geom.asynch_index` from default mlr columns.

**References**:
- `rph_core/steps/step4_features/extractors/geometry.py` (current r1/r2/asynch/asynch_index)
- `rph_core/steps/step4_features/schema.py` (fixed columns and units)

**Acceptance Criteria**:
- New columns appear in `features_raw.csv`.
- Default `features_mlr.csv` prefers `geom.r_avg`, `geom.dr`, `geom.close_contacts_density`.

---

### 4. V6.1: TS Quality Parsing (No New QC)

**What to do**:
- Add new extractor plugin (e.g., `ts_quality`) that parses *existing* TS output to compute:
  - `ts.n_imag`
  - `ts.imag1_cm1_abs` (absolute value of the most negative frequency)
  - `ts.dipole_debye` (dipole magnitude)
- Parsing sources (prefer in this order):
  1) `context.ts_qm_output` (Gaussian .log or ORCA .out)
  2) `context.ts_log` / `context.ts_orca_out`
- Ensure S4 never triggers QC runs; parsing is text-only.
- Parsing rules (explicit):
  - Gaussian:
    - Frequencies: collect all values from lines matching `Frequencies --` in the final frequency block; `n_imag = count(freq < 0)` and `imag1_cm1_abs = abs(min(freq))` if any.
    - Dipole: parse the last occurrence of the `Dipole moment (field-independent basis, Debye):` block; use the `Tot=` value if present, else compute magnitude from X/Y/Z.
  - ORCA:
    - Frequencies: parse the `VIBRATIONAL FREQUENCIES` section (cm-1); treat negative values as imaginary.
    - Dipole: parse the last `Total Dipole Moment` (Debye) value or the `DIPOLE MOMENT` section if present.

**References**:
- `rph_core/steps/step4_features/context.py` (available paths)
- `rph_core/steps/step4_features/log_parser.py` (Gaussian orbitals/frequencies parser; can be reused/extended)
- `rph_core/steps/step4_features/feature_miner.py` (plugin orchestration)

**Acceptance Criteria**:
- Columns exist in `features_raw.csv`.
- If TS output is missing, columns are empty and a meta missing_reason is recorded (no crash).

**Default enablement**:
- Add `ts_quality` to Step4 `DEFAULT_ENABLED_PLUGINS` so the default MLR column `ts.imag1_cm1_abs` is populated when TS output exists.

---

### 4b. V6.1: Enforce “S4 Never Runs QC” (Tests + Guardrails)

**What to do**:
- Keep Step4 default execution in “extract-only” mode:
  - Ensure `DEFAULT_ENABLED_PLUGINS` stays limited to extract-only plugins.
  - Ensure any plugin that can generate `job_specs` never executes them when `job_run_policy=disallow`.
- Add an additional static safety guard (cheap and strong):
  - Introduce a plugin attribute like `can_submit_jobs` (default false).
  - If `job_run_policy=disallow` and a plugin has `can_submit_jobs=true`, the plugin may still emit `job_specs` but MUST NOT execute them.
  - In this case, the plugin should be forced into “generate-only” behavior:
    - record job_specs
    - set plugin status to `SKIPPED` (or a dedicated `GENERATE_ONLY` status if you add it) and add `missing_reasons=["job_run_policy_disallow"]`
    - emit no computed features
- Add a hard safety test that fails if S4 triggers QC execution:
  - Monkeypatch `subprocess.run` and QC execution entrypoints in `rph_core/utils/qc_interface.py` / `rph_core/utils/qc_task_runner.py` / `rph_core/utils/orca_interface.py` to raise.
  - Run `FeatureMiner.run(...)` with minimal valid inputs and assert it still writes S4 outputs successfully.
  - Concrete entrypoints to patch (repo-real):
    - `rph_core.utils.qc_task_runner.QCTaskRunner.run_ts_opt_cycle`
    - `rph_core.utils.qc_task_runner.QCTaskRunner.run_opt_sp_cycle`
    - `rph_core.utils.qc_interface.GaussianRunner.run`
    - `rph_core.utils.qc_interface.GaussianInterface.optimize`
    - `rph_core.utils.orca_interface.ORCAInterface._run_orca`
    - `subprocess.run` (global patch as last resort)
- Optional additional guardrail (cheap and effective): add a unit test that AST-greps for forbidden calls from within `rph_core/steps/step4_features/` (e.g., `QCTaskRunner(`, `run_in_sandbox(`) and fails if found.

**References**:
- `rph_core/steps/step4_features/feature_miner.py` (default plugins + job_run_policy)
- `rph_core/steps/step4_features/extractors/interaction_analysis.py` (job_specs generation)
- `rph_core/steps/step4_features/fragment_extractor.py` (known anti-pattern: runs Gaussian/xTB; must remain unused)

**Acceptance Criteria**:
- `python -m pytest -q` passes with the new “no QC execution” test.

---

### 5. V6.2: Add S3 Post-QC Enrichment (DIAS/ASM) + Artifact Contract

**What to do**:
- Implement `S3_PostQCEnrichment/` stage that runs after TS success and writes:
- Implement `S3_PostQCEnrichment/` stage that runs after TS success and writes under the Step3 output directory (`work_dir/S3_TransitionAnalysis/` today):
  - `<S3_DIR>/S3_PostQCEnrichment/enrichment.json`
  - `<S3_DIR>/S3_PostQCEnrichment/enrichment_status.json`
- Enrichment must be executed from S3 and must reuse `qc_interface`/`orca_interface` or `QCTaskRunner` (no ad-hoc subprocess or cwd hacks).
- Integration point and failure policy (explicit):
  - Invoke enrichment from `rph_core/steps/step3_opt/ts_optimizer.py:TSOptimizer.run_with_qctaskrunner()` after:
    - TS optimization is converged,
    - `ts_final.xyz` has been written to `<S3_DIR>/ts_final.xyz`, and
    - `SPMatrixReport` is constructed (so `e_ts_final` is available for `ts_TS`).
  - If enrichment fails:
    - Step3 must still succeed overall.
    - Write `<S3_DIR>/S3_PostQCEnrichment/enrichment_status.json` with status=failed and an error message.
    - Step4’s asm plugin must treat this as `missing_reason=enrichment_failed` and output NaNs.

- Implementation placement (explicit, avoid bloating ts_optimizer.py):
  - Add new module: `rph_core/steps/step3_opt/post_qc_enrichment.py`
  - Expose a single entrypoint function, e.g. `run_post_qc_enrichment(s3_dir: Path, config: dict, sp_report: SPMatrixReport, reactant_complex_xyz: Path) -> None`.
  - `ts_optimizer.py` should only call this function and handle exceptions.
- Configuration surface (explicit defaults):
  - Config file location for defaults: `config/defaults.yaml`.
  - `step3.enrichment.enabled` (default false)
  - `step3.enrichment.force_rerun` (default false)
  - `step3.enrichment.connectivity_scale` (default 1.25)
  - `step3.enrichment.bond_min_dist_angstrom` (default 0.6)
  - `step3.enrichment.write_dirname` (default `S3_PostQCEnrichment`)
  - `step3.enrichment.fragment_indices_override` (default null; when set, bypass connectivity split)
  - `step3.enrichment.compute_ts_ts_orca` (default false; if true, run TS@TS ORCA SP for interaction closure)
- Minimum JSON contract fields (suggested):
  - `schema_version`, `orca_template_id`, `orca_template_hash`
    - Definition in this repo (as implemented today): `rph_core/utils/orca_interface.py` renders input from method/basis/aux_basis/solvent/nprocs/maxcore.
    - `orca_template_id`: a stable string (e.g., `orca:{method}/{basis}:{aux_basis}:{solvent}:{nprocs}:{maxcore}`)
    - `orca_template_hash`: sha256 of the exact rendered ORCA input file content excluding the coordinate block (if excluding coordinates is too invasive, hash the full `.inp` content and also record `geometry_hash` separately).
  - `inputs`: hashes of geometries and fragment definitions
  - `energies`: at least `frag1_R`, `frag2_R`, `frag1_TS`, `frag2_TS` (Hartree)
  - Optional: `ts_TS` for interaction closure (Hartree)
  - `units`: explicit per-field (`hartree`, `ev`, etc.)
  - `status`: per-subtask status, error reasons

**Enrichment contract schema (REQUIRED, because caching depends on it)**:
- `<S3_DIR>/S3_PostQCEnrichment/enrichment.json` must contain (top-level required keys):
  - `schema_version`: string (e.g., `"enrichment_v1"`)
  - `created_at`: ISO timestamp
  - `orca_template_id`: string
  - `orca_template_hash`: sha256 hex
  - `geometry_hashes`:
    - `reactant_complex_xyz_sha256`: sha256 of raw bytes of `S2_Retro/reactant_complex.xyz`
    - `ts_final_xyz_sha256`: sha256 of raw bytes of `<S3_DIR>/ts_final.xyz`
  - `fragments`:
    - `frag1_indices`: sorted int list
    - `frag2_indices`: sorted int list
    - `fragment_source`: `override|connectivity`
    - `fragment_hash`: sha256 of JSON({frag1_indices, frag2_indices}) with sort_keys
  - `energies_hartree`:
    - `frag1_R`, `frag2_R`, `frag1_TS`, `frag2_TS`
    - optional `ts_TS` (ONLY when `step3.enrichment.compute_ts_ts_orca=true`; must be computed with the same ORCA template as fragments)
  - `units`: fixed string `{"energies_hartree": "hartree"}`

- `<S3_DIR>/S3_PostQCEnrichment/enrichment_status.json` must contain:
  - `schema_version`: string (e.g., `"enrichment_status_v1"`)
  - `created_at`: ISO timestamp
  - `cache_key`: sha256 of JSON({orca_template_hash, geometry_hashes, fragment_hash}) with sort_keys
  - `cache_hit`: boolean
  - `status`: `ok|skipped|failed`
  - `error`: optional string

**Fragment definition strategy (confirmed)**:
- Use **connectivity split**: infer two fragments by building a bond graph from the reactant complex geometry and extracting connected components.
- Guardrails:
  - Only proceed if exactly 2 components are found; otherwise write `enrichment_status.json` with a clear error and skip enrichment.
  - Record the inferred atom indices for each fragment into the contract JSON.
 - Override escape hatch (MUST HAVE for coverage):
   - If `step3.enrichment.fragment_indices_override` is provided, bypass connectivity split entirely.
   - Record `fragment_source=override` in the contract.

**Connectivity split algorithm (MUST be deterministic)**:
- Input geometry for inference: `S2_Retro/reactant_complex.xyz` (same atom ordering used downstream).
- Build an undirected graph over atoms using a distance-based bond heuristic:
  - For each pair (i, j), compute distance `d_ij`.
  - Lookup covalent radii `r_cov[element]` for both atoms.
  - Add an edge if `d_ij <= SCALE * (r_cov[i] + r_cov[j])`.
  - Use `SCALE = 1.25` (tunable via config) and ignore pairs with `d_ij < 0.6 Å` (pathological overlaps).
  - Element radii source: embed an explicit covalent radii table (Å) for the elements expected in this project’s organic reactions; if an element is unknown, abort inference (status=error, missing_reason=`unknown_element_radius`).
    - H 0.31
    - B 0.85
    - C 0.76
    - N 0.71
    - O 0.66
    - F 0.57
    - Si 1.11
    - P 1.07
    - S 1.05
    - Cl 1.02
    - Br 1.20
    - I 1.39
- Extract connected components.
- Deterministic fragment labeling:
  - `frag1` is the component that contains the smallest atom index; `frag2` is the other.
  - Store both the sorted atom indices list and a stable fingerprint (e.g., sha1 of indices).

**ASM/DIAS energy semantics (no ambiguity)**:
- Define geometries for fragment single points:
  - `frag*_R`: isolated fragment cut from `S2_Retro/reactant_complex.xyz` using the inferred atom indices (no relaxation).
- `frag*_TS`: isolated fragment cut from `<S3_DIR>/ts_final.xyz` (Step3 canonical TS) using the *same* atom indices (no relaxation).
- Define full-system TS energy used for interaction closure:
  - `ts_TS` is ONLY defined when `step3.enrichment.compute_ts_ts_orca=true` and is computed with the same ORCA template as the fragment SPs.
  - Otherwise omit `ts_TS` and omit `asm.interaction_kcal` downstream.
- Derived quantities (S4 consumption):
  - `asm.distortion_total_kcal = Σ[(E_frag_TS − E_frag_R) * 627.509]`
  - `asm.interaction_kcal = (E_ts_TS − Σ(E_frag_TS)) * 627.509` (only if `ts_TS` is present and explicitly marked as Hartree).

**Toxic-path-safe ORCA execution (explicit)**:
- Implementation choice (pin down): extend `rph_core/utils/orca_interface.py` so ORCA execution uses `rph_core/utils/qc_interface.py:is_path_toxic()` (spaces + `[](){} `), not just whitespace.
  - In `_run_orca(...)`, if `is_path_toxic(output_dir)` is true:
    - create a qc_interface `LinuxSandbox` directory
    - copy the `.inp` into the sandbox
    - run ORCA with `cwd=sandbox_dir`
    - copy the `.out` and any auxiliary files back to the real `output_dir`
    - cleanup sandbox
  - Otherwise run ORCA in-place in `output_dir`.
- This ensures enrichment inherits the repo-wide toxic-path policy without adding step-local hacks.

**References**:
- `rph_core/steps/step3_opt/ts_optimizer.py` (S3 lifecycle + where to hook)
- `rph_core/utils/orca_interface.py` (ORCA execution wrapper)
- `rph_core/utils/qc_interface.py` (toxic path + sandbox rules)
- `rph_core/steps/step4_features/fragment_extractor.py` (legacy fragment splitting + QC-in-step4 anti-pattern; reuse only non-QC splitting ideas)

**Acceptance Criteria**:
- With enrichment enabled, at least one sample produces `<S3_DIR>/S3_PostQCEnrichment/enrichment.json` and `<S3_DIR>/S3_PostQCEnrichment/enrichment_status.json`.
- With enrichment disabled, pipeline behavior matches previous version (no new files required).
- `python -m pytest -q` passes.

---

### 6. V6.2: S4 Consumes Enrichment Contract (No Compute)

**What to do**:
- Add S4 extractor plugin (e.g., `asm_enrichment`) that reads `S3_PostQCEnrichment/enrichment.json` and derives:
  - `asm.distortion_total_kcal = Σ(E_frag_TS − E_frag_R) * 627.509`
  - Optional `asm.interaction_kcal = E_TS − Σ(E_frag_TS)` ONLY if `ts_TS` exists AND was computed via `step3.enrichment.compute_ts_ts_orca=true` (same ORCA template).
- If enrichment contract missing or invalid: leave features empty and record missing_reason in meta.
- S4 must locate Step3 output directory robustly:
  - Prefer `work_dir/S3_TransitionAnalysis/` (current orchestrator).
  - Fallback to `work_dir/S3_TS/` for legacy runs.
  - Do not do global tree scans; only check these known step directories.
  - Wiring (explicit):
    - Add `work_dir: Path` and `s3_dir: Optional[Path]` fields to `rph_core/steps/step4_features/context.py:FeatureContext`.
    - In `rph_core/steps/step4_features/feature_miner.py`, set `work_dir = output_dir.parent` and resolve `s3_dir` before plugin execution.
    - Pass `s3_dir` in the context so plugins can read `<S3_DIR>/S3_PostQCEnrichment/...` and dipolar artifacts without guessing.

**References**:
- `rph_core/steps/step4_features/feature_miner.py` (plugin execution)
- `rph_core/steps/step4_features/context.py` (where S3 dir can be resolved; may need path injection)
- `rph_core/steps/step4_features/mech_packager.py` (pattern for locating S3 directories)

**Acceptance Criteria**:
- When contract exists, S4 outputs asm columns without running any QC.
- Add fixture-based unit tests (no QC required):
  - `tests/fixtures/enrichment/enrichment.json` (minimal valid contract)
  - Test asserts `asm.distortion_total_kcal` is computed exactly.
  - Separate test asserts `asm.interaction_kcal` is computed exactly when `ts_TS` is present.

---

### 7. V6.3: Add Dipolar FMO/CDFT Parser (S4 Text Parse Only)

**What to do**:
- Define the dipolar artifact location contract (provenance):
  - If a dipolar intermediate output exists for a reaction, it MUST be placed under `<S3_DIR>/dipolar/`.
  - Supported filenames: `*.log` (Gaussian) and `*.out` (ORCA) containing `dipolar` in the filename.
  - This plan does not add new computations to generate dipolar; it only indexes/parses if present.

- Add Step3 artifact index writer:
  - Always write `<S3_DIR>/artifacts_index.json`.
  - If dipolar artifacts exist under `<S3_DIR>/dipolar/`, select one deterministically (same fallback rules as S4) and record `dipolar.path_rel` + `dipolar.sha256`.

- Add S4 plugin `fmo_cdft_dipolar` that:
  - Locates dipolar intermediate output produced by S3.
  - Parses `EHOMO`/`ELUMO` (eV) and emits:
    - `fmo.dipolar_homo_ev`
    - `fmo.dipolar_lumo_ev`
    - `fmo.dipolar_gap_ev = LUMO - HOMO`
    - `fmo.dipolar_omega_ev` where:
      - `mu = (EHOMO + ELUMO)/2`
      - `eta = ELUMO - EHOMO` (must be >0)
      - `omega = mu^2/(2*eta)`
- File location strategy:
  - REQUIRED: S3 writes `artifacts_index.json` including dipolar relative path + sha256.
  - S4 must read `artifacts_index.json` first for stable, reproducible selection.
  - Fallback (ONLY when index missing): scan within Step3 output directory only for `*dipolar*.log` / `*dipolar*.out`.
  - Config override (explicit): `step4.fmo_dipolar.path_override` (if set, use exactly this file; no scanning).
  - Default for this override should be present (empty) in `config/defaults.yaml`.
- Multi-file strategy (confirmed):
  - Index path is authoritative.
  - Fallback must be deterministic and filesystem-independent:
    - Filter candidates to those with “normal termination” markers (Gaussian: `Normal termination`; ORCA: `ORCA TERMINATED NORMALLY`).
    - Sort remaining candidate relative paths lexicographically and choose the first.
    - Record `selected_relpath` + `selected_sha256` in meta.

**S3 artifacts index contract (REQUIRED)**:
- `<S3_DIR>/artifacts_index.json`:
  - `schema_version`: `artifacts_index_v1`
  - `created_at`: ISO timestamp
  - `dipolar` (optional block, present only if dipolar artifact exists):
    - `path_rel`: relative path from `<S3_DIR>` to the chosen dipolar output
    - `sha256`: sha256 of raw bytes of that file

**Where `artifacts_index.json` is written (explicit)**:
- In Step3, after Step3 completes (and any external dipolar outputs have been placed under `<S3_DIR>/dipolar/`), write/update `<S3_DIR>/artifacts_index.json`.
- If dipolar does not exist, still write the index file (without `dipolar` block) so downstream behavior is deterministic.

**Parsing notes**:
- Gaussian: parse `Alpha  occ. eigenvalues --` and `Alpha virt. eigenvalues --` (Hartree), then convert to eV using 27.2114.
  - Deterministic selection rule:
    - HOMO: last numeric value on the **last** matched `Alpha  occ. eigenvalues --` line in the file.
    - LUMO: first numeric value on the **first** `Alpha virt. eigenvalues --` line that appears *after* that last occ line; if positional matching is too complex, fall back to using the **last** matched virt line.
  - Existing implementation exists in `rph_core/steps/step4_features/log_parser.py` but must be updated to match this rule.
- ORCA: parse the **last** `ORBITAL ENERGIES` table in the file.
  - Identify HOMO as the last orbital row with `OCC > 0.0`.
  - Identify LUMO as the first orbital row after HOMO with `OCC == 0.0`.
  - Prefer the `E(eV)` column if present; otherwise use `E(Eh)` and convert to eV.
  - Sanity checks:
    - Require both HOMO and LUMO parseable.
    - Require `gap = LUMO - HOMO > 0`; if not, mark invalid and leave omega empty.

**References**:
- `rph_core/steps/step4_features/log_parser.py` (Gaussian orbitals parsing; conversion)
- `rph_core/steps/step4_features/electronic_extractor.py` (CDFT indices formula; may be reused conceptually)

**Acceptance Criteria**:
- When dipolar output exists, plugin produces consistent HOMO/LUMO/gap/omega values across repeated runs.
- When missing, plugin does not error; it records missing_reason and leaves columns empty.
- Add fixture-based unit tests (so this is verifiable without requiring S3 to generate dipolar outputs):
  - `tests/fixtures/dipolar/gaussian_dipolar.log` (minimal snippet with occ/virt eigenvalues)
  - `tests/fixtures/dipolar/orca_dipolar.out` (minimal snippet with ORBITAL ENERGIES table)
  - Tests assert HOMO/LUMO/gap/omega parse AND index-priority selection logic.
  - Tests cover deterministic fallback selection (no mtime reliance).

---

### 8. V6.4: Hardening (Caching/Resume/Strict QC Mode)

**What to do**:
- S3 enrichment caching:
  - Derive cache key from `orca_template_hash + geometry_hash + fragment_hash`.
  - Implement cache key computation in `rph_core/utils/cache_key.py` and reuse it from S3/S4.
  - Skip re-running enrichment if the same key already has a successful result, unless `force_rerun_enrichment=true`.
- S4 strict-QC mode (configurable):
  - Use `ts.n_imag != 1`, invalid forming_bonds, etc. to either:
    - keep in raw but drop from mlr, and/or
    - set `sample_weight` in mlr output.
- Disk bloat guardrails:
  - Keep only required ORCA outputs for enrichment; optionally compress/archive large logs.

- Cache-key helper unification (avoid S3→S4 coupling):
  - Create `rph_core/utils/cache_key.py` as the single source of truth for cache keys.
  - Expose two explicit helpers (different formats by design):
    - `generate_step4_plugin_cache_key(...) -> str`: preserves current Step4 behavior (sha1 hex truncated to 16 chars; includes file fingerprints with `mtime`), to keep `tests/test_step4_cache_key.py` consistent.
    - `generate_enrichment_cache_key(...) -> str`: returns full sha256 hex of the JSON payload `{orca_template_hash, geometry_hashes, fragment_hash}` (content-hash-based; no mtime).
  - Update `rph_core/steps/step4_features/schema.py:generate_cache_key()` to call `generate_step4_plugin_cache_key(...)`.

- Disk bloat guardrails (make measurable):
  - Add config: `step3.enrichment.cleanup_aux_files` (default true).
  - On success, keep ONLY:
    - `<S3_DIR>/S3_PostQCEnrichment/enrichment.json`
    - `<S3_DIR>/S3_PostQCEnrichment/enrichment_status.json`
    - ORCA inputs/outputs needed for audit (`*.inp`, `*.out`)
  - Delete all other auxiliary files in `<S3_DIR>/S3_PostQCEnrichment/` unless `cleanup_aux_files=false`.

**References**:
- Create and use a shared cache key helper:
  - New module: `rph_core/utils/cache_key.py`
  - Both S3 enrichment and any future S4 caching import from this module (avoid S3 depending on S4).
- `rph_core/utils/checkpoint_manager.py` (resume semantics)

**Acceptance Criteria**:
- `--resume` run does not duplicate enrichment work when outputs already exist, verified by at least one of:
- `<S3_DIR>/S3_PostQCEnrichment/enrichment_status.json` contains the same `cache_key` and a `cache_hit: true` flag on the second run, and
- no new ORCA run subdirectories are created under `<S3_DIR>/S3_PostQCEnrichment/` between runs (directory listing stable), and
  - `enrichment_status.json` timestamps do not change on cache hit.

---

### 9. V6.4: Modeling Companion Script (LOOCV/VIF/LFER; not pipeline)

**What to do**:
- Add `scripts/ml/train_mlr_loocv.py`:
  - Load `features_mlr.csv`.
  - Run LOOCV.
  - Compute VIF (default threshold 5) and/or perform Lasso-based feature selection.
  - Optional fold-internal LFER: fit `ΔG ≈ a·ΔE + C` on train fold and impute missing ΔG in test fold without leakage.
- Keep this script out of step execution; it’s a post-analysis utility only.

**References**:
- `S4_Data/features_mlr.csv` (new contract)

**Acceptance Criteria**:
- Script runs on a small sample directory and emits a report (stdout + optional CSV/JSON in `scripts/ml/out/`).

---

### 10. Self Gap Check (Contracts & No Strays)

**What to do**:
- Run a repo-wide check (excluding `rph_core_backup_20260115/`) to ensure no runtime code or tests still require `features.csv`.
- Verify new contracts are referenced consistently:
  - `features_raw.csv` / `features_mlr.csv` / `feature_meta.json` in Step4 code and tests.
  - `<S3_DIR>/S3_PostQCEnrichment/enrichment.json` in Step3/S4 linkage.

**Acceptance Criteria**:
- `python -m pytest -q` passes.
- Grep for `features.csv` outside backup/notes is empty.

---

### 11. Docs & Contract Acceptance

**What to do**:
- Add/Update documentation (suggested): `V6_ACCEPTANCE.md`:
  - New S4 output contract and meaning.
  - Enrichment contract schema.
  - Dipolar plugin schema.
  - Recommended default mlr columns.
  - Verification commands.

**Acceptance Criteria**:
- Doc clearly states IN/OUT, contracts, and reproducibility expectations.

---

## Appendix: `feature_meta.json` Example (V6.1)

```json
{
  "meta": {
    "schema_version": "6.1",
    "schema_signature": "<sha1>",
    "feature_status": "OK",
    "method": "<method>",
    "solvent": "<solvent>",
    "temperature_K": 298.15,
    "enabled_plugins": ["thermo", "geometry", "qc_checks", "ts_quality"],
    "outputs": {
      "features_raw_csv": "features_raw.csv",
      "features_mlr_csv": "features_mlr.csv",
      "feature_meta_json": "feature_meta.json"
    },
    "units": {
      "sample_id": null,
      "thermo.dE_activation": "kcal/mol",
      "geom.r_avg": "Å",
      "ts.imag1_cm1_abs": "cm^-1",
      "qc.sample_weight": null
    },
    "mlr_columns": [
      "sample_id",
      "thermo.dE_activation",
      "thermo.dE_reaction",
      "geom.r_avg",
      "geom.dr",
      "geom.close_contacts_density",
      "ts.imag1_cm1_abs",
      "asm.distortion_total_kcal",
      "fmo.dipolar_omega_ev",
      "qc.sample_weight"
    ],
    "missing_stats": {
      "asm.distortion_total_kcal": {"is_missing": true, "is_invalid": false, "reason": "missing_input"},
      "fmo.dipolar_omega_ev": {"is_missing": true, "is_invalid": false, "reason": "plugin_disabled"}
    },
    "policies": {
      "job_run_policy": "disallow",
      "nics_trigger_policy": "generate_only"
    }
  },
  "trace": {
    "inputs_fingerprint": {
      "ts_xyz": {"size": 123, "mtime": 0, "sha256_prefix": "deadbeef"}
    },
    "plugins": {
      "thermo": {
        "status": "OK",
        "runtime_ms": 12.3,
        "warnings": [],
        "errors": [],
        "missing_reasons": []
      },
      "fmo_cdft_dipolar": {
        "status": "SKIPPED",
        "missing_reasons": ["plugin_disabled"]
      }
    }
  }
}
```

---

## Defaults Applied (override if needed)
- Keep `features_mlr.csv` default columns ≤10, starting from:
  - `thermo.dE_activation`, `thermo.dE_reaction`, `geom.r_avg`, `geom.dr`, `geom.close_contacts_density`, `ts.imag1_cm1_abs`
  - plus at most one of: `asm.distortion_total_kcal` (if V6.2 enabled) or `fmo.dipolar_omega_ev`/`fmo.dipolar_gap_ev` (if V6.3 enabled)
- FMO dipolar plugin default: OFF initially.

## Decisions Needed
- None (blocking decisions resolved).
