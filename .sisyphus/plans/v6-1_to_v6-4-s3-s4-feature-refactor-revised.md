# V6.1–V6.4 Plan: S3 Compute-First, S4 Extract-Only, Linear-Model Features (REVISED - Tether-Cut)

## Context

### Original Request
Implement a V6.1–V6.4 modification program optimized for small-sample ([5+2], N≈30) + linear models:
- Move **all computations** into S3 (ORCA single points; DIAS/ASM enrichment; optional TS@TS closure).
- Make S4 **strictly extraction/processing/QC marking**; S4 must not submit any QC jobs.
- Add a lightweight **FMO/CDFT parser** in S4 based on S3 "dipolar intermediate" outputs: parse HOMO/LUMO and compute `gap` and electrophilicity `omega` (no new compute).
- Lock S4 outputs (V6.1 contract): `features_raw.csv`, `features_mlr.csv` (6–10 columns), `feature_meta.json`.
- Follow repo anti-pattern constraints: no duplicated QC runner in steps; reuse `rph_core/utils/qc_interface.py` toxic-path sandbox logic.

### Critical Revision Required (V6.2)
The original plan used **"connectivity split"** fragmenter (direct connected components) which will **FAIL for intramolecular [5+2] systems**:
- Intramolecular oxidopyrylium cycloadditions form a **single connected component**
- Direct connectivity split would find only 1 component → enrichment always skipped/failed
- ASM features would always be empty

**Revised V6.2 strategy: "intramolecular tether-cut"**
- Use forming_bonds + graph topology to identify dipole core (5-atom path)
- Cut the tether-dipole connecting bond
- Apply H-caps to both sides
- This properly handles intramolecular systems while preserving dipole electronic structure

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
- Forming bonds detection already available in Step2/Step3 via `forming_bonds` tuple:
  - `rph_core/orchestrator.py:49` - type definition
  - `rph_core/steps/step4_features/schema.py:21` - validation
  - `rph_core/steps/step4_features/extractors/geometry.py` - uses for r1/r2

### Metis Review
- Metis consultation was attempted but the metis agent tool call failed with an infra error (`JSON Parse error: Unexpected EOF`). This plan includes an explicit "Self gap check" task to compensate.
- **Critical gap identified (post-Metis)**: Original V6.2 "connectivity split" fails for intramolecular systems. This revision addresses it with tether-cut strategy.

---

## Work Objectives

### Core Objective
Deliver a clean S3/S4 separation where S3 performs all enrichment computations (using intramolecular tether-cut fragmenter for ASM) and S4 produces reproducible, linear-model-friendly features (raw + MLR-ready + metadata) without submitting QC jobs.

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
  - Optional ASM/DIAS features in S4 by consuming an S3-written enrichment contract (with tether-cut fragmenter).
  - Optional dipolar FMO/CDFT features in S4 by parsing S3 dipolar outputs (no compute).
- **NEW**: S3 intramolecular tether-cut fragmenter (`rph_core/steps/step3_opt/intramolecular_fragmenter.py`)
- **NEW**: Molecular graph utilities (`rph_core/utils/molecular_graph.py`)
- **NEW**: Fragment manipulation utilities (`rph_core/utils/fragment_manipulation.py`)
- S3 enrichment contract (written under the Step3 output directory; currently `work_dir/S3_TransitionAnalysis/`):
  - `<S3_DIR>/S3_PostQCEnrichment/enrichment.json` (with fragmenter block, charges, multiplicities)
  - `<S3_DIR>/S3_PostQCEnrichment/enrichment_status.json`
- S3 artifacts index (for reproducible dipolar selection):
  - `<S3_DIR>/artifacts_index.json`
- V6.4 modeling companion script(s) under `scripts/ml/` (kept out of pipeline execution).
- Docs/tests updated to reflect new contracts and tether-cut fragmenter.

### Definition of Done
- `python -m pytest -q` passes.
- `python run_auto_test.py --step ALL` completes for at least one local sample run directory.
- No Step4 code path submits QC jobs (no calls into `QCTaskRunner` / `qc_interface` execution functions from S4).
- Resume logic recognizes the new Step4 outputs (no dependency on `features.csv`).
- `python run_auto_test.py --step ALL --resume` passes with `step4.mechanism_packaging.enabled=true` (ensures checkpoint/resume/packager no longer depends on `features.csv`).
- **NEW**: Intramolecular [5+2] samples produce valid ASM features (not always empty/skipped).

### Must NOT Have (Guardrails)
- No QC execution in S4 (including fragment SP jobs).
- No duplicate ad-hoc sandbox/toxic-path handling; always reuse `rph_core/utils/qc_interface.py` utilities.
- No changes to Step2 output contract (`S2_Retro/ts_guess.xyz` and `S2_Retro/reactant_complex.xyz` must remain).
- Avoid large refactors; touch only what is necessary to implement V6.1–V6.4.
- V6.3 does NOT add new "dipolar intermediate" computations; it only indexes/parses dipolar outputs if they already exist under the Step3 output directory.
- **NEW CRITICAL**: Do NOT use "connectivity split" (direct components) - always use "intramolecular tether-cut" for ASM fragmenter.

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

V6.1 (S4 outputs + feature schema fix) → V6.2 (S3 enrichment: tether-cut fragmenter + contract + S4 consume) → V6.3 (dipolar FMO parser + plugin) → V6.4 (hardening + ML scripts + docs)

---

## TODOs

### 0. Baseline Inventory & Contract Snapshot

**What to do**:
- Record current Step4 contract and resume semantics.
- Apply the confirmed migration policy: **no backward compatibility** for `S4_Data/features.csv`.
- Enumerate and update all direct `features.csv` expectations in code/tests/docs to the new contract.
- **Explicit allowlist for `features.csv` grep check**: Since `features.csv` appears in many non-backup locations (docs, historical files, backup directories), define an allowlist:
  - Allow: `rph_core_backup_20260115/` (historical backup, not source code)
  - Exclude: `*.backup_*` files in runtime directories
  - Update grep check to use allowlist; if match is in allowlist, count as pass.

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
- `feature_status` values are the `.value` strings from `rph_core/steps/step4_features/status.py:FeatureStatus` (e.g., `ok`, `missing_inputs`, `invalid_inputs`). **NOTE**: These are lowercase strings, e.g., `"ok"`, not `"OK"`.
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
  - Update `rph_core/utils/checkpoint_manager.py:is_step4_complete()` to validate "non-empty" by:
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
- Add explicit "best available" convenience columns in raw output only (not default mlr columns):
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

### 4b. V6.1: Enforce "S4 Never Runs QC" (Tests + Guardrails)

**What to do**:
- Keep Step4 default execution in "extract-only" mode:
  - Ensure `DEFAULT_ENABLED_PLUGINS` stays limited to extract-only plugins.
  - Ensure any plugin that can generate `job_specs` never executes them when `job_run_policy=disallow`.
- Add an additional static safety guard (cheap and strong):
  - Introduce a plugin attribute like `can_submit_jobs` (default false).
  - If `job_run_policy=disallow` and a plugin has `can_submit_jobs=true`, the plugin may still emit `job_specs` but MUST NOT execute them.
  - In this case, the plugin should be forced into "generate-only" behavior:
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
- **Explicit enforcement location (pin down)**:
  - Document that `can_submit_jobs` + `job_run_policy=disallow` logic is implemented in:
    - `rph_core/steps/step4_features/extractors/base.py` (base plugin class), OR
    - `rph_core/steps/step4_features/feature_miner.py` (orchestrator).
  - This ensures implementers know exactly where the guardrail lives.

**References**:
- `rph_core/steps/step4_features/feature_miner.py` (default plugins + job_run_policy)
- `rph_core/steps/step4_features/extractors/interaction_analysis.py` (job_specs generation)
- `rph_core/steps/step4_features/fragment_extractor.py` (known anti-pattern: runs Gaussian/xTB; must remain unused)

**Acceptance Criteria**:
- `python -m pytest -q` passes with the new "no QC execution" test.

---

### 5. V6.2: Implement Molecular Graph Utilities (NEW MODULE)

**What to do**:
- Create new module: `rph_core/utils/molecular_graph.py`
- Implement graph utilities for tether-cut fragmenter:

```python
def build_bond_graph(coords: np.ndarray, symbols: List[str],
                   scale: float = 1.25, min_dist: float = 0.6) -> Dict[int, List[int]]:
    """
    Build adjacency list graph using covalent radius-based bond heuristic.
    
    Args:
        coords: (N, 3) coordinates in Å
        symbols: List of element symbols
        scale: Distance threshold multiplier
        min_dist: Minimum distance to ignore (pathological overlaps)
    
    Returns:
        Adjacency list: {atom_idx: [neighbor_indices]}
    
    Raises:
        ValueError: If unknown element encountered
    """
    # Covalent radii table (Å)
    COVALENT_RADII = {
        'H': 0.31, 'B': 0.85, 'C': 0.76, 'N': 0.71, 'O': 0.66,
        'F': 0.57, 'Si': 1.11, 'P': 1.07, 'S': 1.05,
        'Cl': 1.02, 'Br': 1.20, 'I': 1.39
    }
    
    n_atoms = len(coords)
    graph = {i: [] for i in range(n_atoms)}
    
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            d_ij = GeometryUtils.calculate_distance(coords, i, j)
            
            # Skip too close (pathological)
            if d_ij < min_dist:
                continue
            
            # Covalent radius sum threshold
            r_i = COVALENT_RADII.get(symbols[i])
            r_j = COVALENT_RADII.get(symbols[j])
            
            if r_i is None or r_j is None:
                raise ValueError(f"Unknown element radius: {symbols[i]}, {symbols[j]}")
            
            threshold = scale * (r_i + r_j)
            
            if d_ij <= threshold:
                graph[i].append(j)
                graph[j].append(i)
    
    return graph

def get_connected_components(graph: Dict[int, List[int]]) -> List[List[int]]:
    """
    Extract all connected components from adjacency list graph using BFS.
    
    Args:
        graph: Adjacency list
    
    Returns:
        List of components, each is a list of atom indices
    """
    visited = set()
    components = []
    
    for atom in graph:
        if atom not in visited:
            # BFS from this atom
            component = []
            queue = [atom]
            visited.add(atom)
            
            while queue:
                node = queue.pop(0)
                component.append(node)
                
                for neighbor in graph[node]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)
            
            components.append(component)
    
    return components

def find_shortest_path(graph: Dict[int, List[int]], start: int, end: int) -> List[int]:
    """
    Find shortest path between two atoms using BFS.
    
    Args:
        graph: Adjacency list
        start: Start atom index
        end: End atom index
    
    Returns:
        List of atom indices representing path (inclusive)
        Returns empty list if no path exists
    """
    if start == end:
        return [start]
    
    # BFS with parent tracking
    queue = [start]
    parent = {start: None}
    visited = {start}
    
    while queue:
        node = queue.pop(0)  # BFS uses queue
        
        if node == end:
            # Reconstruct path
            path = []
            while node is not None:
                path.append(node)
                node = parent[node]
            return path[::-1]
        
        for neighbor in graph[node]:
            if neighbor not in visited:
                visited.add(neighbor)
                parent[neighbor] = node
                queue.append(neighbor)
    
    return []  # No path found
```

**Acceptance Criteria**:
- Unit tests for each function in `tests/test_molecular_graph.py`.
- Graph building correctly identifies bonds for simple molecules.
- Shortest path finds correct routes.
- Connected components extracted correctly.

---

### 6. V6.2: Implement Fragment Manipulation Utilities (NEW MODULE)

**What to do**:
- Create new module: `rph_core/utils/fragment_manipulation.py`
- Implement H-capping utility:

```python
def h_cap_fragment(
    coords: np.ndarray,
    symbols: List[str],
    cap_positions: List[Tuple[int, np.ndarray]],  # [(atom_idx, cap_direction_vec), ...]
    cap_bond_lengths: Dict[str, float] = None
) -> Tuple[np.ndarray, List[str]]:
    """
    Add hydrogen atoms to cap open valences.
    
    Args:
        coords: (N, 3) coordinates
        symbols: List of element symbols
        cap_positions: List of (atom_idx, cap_direction_vec)
            - atom_idx: Index of atom to cap
            - cap_direction_vec: Vector pointing from atom toward cap position
              (should be colinear with bond being cut, pointing away from cut bond)
        cap_bond_lengths: Bond length by element (Å)
            Default: {'C': 1.09, 'N': 1.01, 'O': 0.96, 'H': 0.76}
    
    Returns:
        (capped_coords, capped_symbols) - Extended coordinates and symbols
    """
    if cap_bond_lengths is None:
        cap_bond_lengths = {
            'C': 1.09, 'N': 1.01, 'O': 0.96, 'H': 0.76,
            'Si': 1.48, 'P': 1.42, 'S': 1.35, 'Cl': 1.27,
            'Br': 1.42, 'I': 1.54
        }
    
    # Add cap atoms
    capped_coords = coords.copy()
    capped_symbols = symbols.copy()
    
    for atom_idx, direction_vec in cap_positions:
        # Normalize direction
        direction = direction_vec / np.linalg.norm(direction_vec)
        
        # Element and bond length
        element = symbols[atom_idx]
        bond_length = cap_bond_lengths.get(element, 1.09)  # Default to C-H
        
        # Cap atom position (colinear, opposite direction of cut bond)
        # direction_vec points from cap atom -> bonded atom, so we go opposite
        cap_pos = coords[atom_idx] - direction * bond_length
        
        # Add cap atom
        capped_coords = np.vstack([capped_coords, cap_pos])
        capped_symbols.append('H')
    
    return capped_coords, capped_symbols
```

**Acceptance Criteria**:
- Unit test: H-cap creates H atom at correct position (colinear, correct bond length).
- Verify cap for different elements (C, N, O).

---

### 7. V6.2: Implement Intramolecular Tether-Cut Fragmenter (NEW MODULE)

**What to do**:
- Create new module: `rph_core/steps/step3_opt/intramolecular_fragmenter.py`
- Implement tether-cut fragmenter class:

```python
from dataclasses import dataclass
from typing import Tuple, List, Optional
import numpy as np

@dataclass
class FragmenterResult:
    """Result of tether-cut fragmentation."""
    fragA_indices: List[int]
    fragB_indices: List[int]
    cut_bond_indices: Tuple[int, int]
    dipole_end_indices: Tuple[int, int]
    alkene_end_indices: Tuple[int, int]
    dipole_core_indices: List[int]
    fragA_coords_R: np.ndarray  # Capped
    fragB_coords_R: np.ndarray  # Capped
    fragA_coords_TS: np.ndarray  # Capped
    fragB_coords_TS: np.ndarray  # Capped
    fragA_symbols_R: List[str]  # Capped
    fragB_symbols_R: List[str]  # Capped
    fragA_symbols_TS: List[str]  # Capped
    fragB_symbols_TS: List[str]  # Capped
    status: str  # "ok" | "failed"
    reason: str = ""  # Error message if failed

class IntramolecularFragmenter(LoggerMixin):
    """
    Intramolecular fragmenter using tether-cut strategy.
    
    For oxidopyrylium [5+2] systems:
    - Fragment A: Dipole (5-atom core + substituents, no tether)
    - Fragment B: Dipolarophile + tether
    
    Uses forming_bonds to identify reaction atoms, then:
    1. Classify alkene vs dipole ends
    2. Find 5-atom dipole core path
    3. Locate and cut tether-dipole connecting bond
    4. H-cap both fragments
    """
    
    def fragment(
        self,
        reactant_coords: np.ndarray,
        reactant_symbols: List[str],
        ts_coords: np.ndarray,
        ts_symbols: List[str],
        forming_bonds: Tuple[Tuple[int, int], Tuple[int, int]],
        config: dict = None
    ) -> FragmenterResult:
        """
        Execute tether-cut fragmentation algorithm (7 steps).
        
        Args:
            reactant_coords: Reactant complex coordinates (N, 3) in Å
            reactant_symbols: Element symbols
            ts_coords: Transition state coordinates (N, 3) in Å
            ts_symbols: Element symbols
            forming_bonds: ((u1, v1), (u2, v2)) forming bond atom pairs
            config: Optional config dict (scale, min_dist, etc.)
        
        Returns:
            FragmenterResult with capped fragment geometries and metadata
        """
        if config is None:
            config = {}
        
        # Step A: Build topology graph
        graph = self._build_graph(
            reactant_coords, reactant_symbols,
            scale=config.get('connectivity_scale', 1.25),
            min_dist=config.get('bond_min_dist_angstrom', 0.6)
        )
        
        # Step B: Classify reaction atoms (alkene vs dipole)
        alkene_pair, dipole_pair = self._classify_reaction_ends(
            graph, reactant_coords, forming_bonds
        )
        if alkene_pair is None or dipole_pair is None:
            return FragmenterResult(
                fragA_indices=[], fragB_indices=[],
                cut_bond_indices=(0, 0),
                dipole_end_indices=(0, 0),
                alkene_end_indices=(0, 0),
                dipole_core_indices=[],
                fragA_coords_R=np.array([]), fragB_coords_R=np.array([]),
                fragA_coords_TS=np.array([]), fragB_coords_TS=np.array([]),
                fragA_symbols_R=[], fragB_symbols_R=[],
                fragA_symbols_TS=[], fragB_symbols_TS=[],
                status="failed",
                reason="failed_to_classify_reaction_ends"
            )
        
        # Step C: Define dipole core (5-atom path)
        dipole_core = self._find_dipole_core_path(graph, dipole_pair)
        if dipole_core is None:
            return FragmenterResult(
                fragA_indices=[], fragB_indices=[],
                cut_bond_indices=(0, 0),
                dipole_end_indices=dipole_pair,
                alkene_end_indices=alkene_pair,
                dipole_core_indices=[],
                fragA_coords_R=np.array([]), fragB_coords_R=np.array([]),
                fragA_coords_TS=np.array([]), fragB_coords_TS=np.array([]),
                fragA_symbols_R=[], fragB_symbols_R=[],
                fragA_symbols_TS=[], fragB_symbols_TS=[],
                status="failed",
                reason="dipole_path_not_5_atoms"
            )
        
        # Step D: Find tether + cut bond
        cut_bond = self._find_cut_bond(
            graph, dipole_core, alkene_pair, reactant_symbols
        )
        if cut_bond is None:
            return FragmenterResult(
                fragA_indices=[], fragB_indices=[],
                cut_bond_indices=(0, 0),
                dipole_end_indices=dipole_pair,
                alkene_end_indices=alkene_pair,
                dipole_core_indices=dipole_core,
                fragA_coords_R=np.array([]), fragB_coords_R=np.array([]),
                fragA_coords_TS=np.array([]), fragB_coords_TS=np.array([]),
                fragA_symbols_R=[], fragB_symbols_R=[],
                fragA_symbols_TS=[], fragB_symbols_TS=[],
                status="failed",
                reason="failed_to_find_cut_bond"
            )
        
        # Step E: Cut bond → connected components
        fragA_indices, fragB_indices = self._cut_and_get_components(
            graph, cut_bond, dipole_pair, alkene_pair
        )
        
        if len(fragA_indices) == 0 or len(fragB_indices) == 0:
            return FragmenterResult(
                fragA_indices=fragA_indices, fragB_indices=fragB_indices,
                cut_bond_indices=cut_bond,
                dipole_end_indices=dipole_pair,
                alkene_end_indices=alkene_pair,
                dipole_core_indices=dipole_core,
                fragA_coords_R=np.array([]), fragB_coords_R=np.array([]),
                fragA_coords_TS=np.array([]), fragB_coords_TS=np.array([]),
                fragA_symbols_R=[], fragB_symbols_R=[],
                fragA_symbols_TS=[], fragB_symbols_TS=[],
                status="failed",
                reason="component_validation_failed"
            )
        
        # Step F: H-cap fragments (R and TS)
        # NOTE: cap_positions are computed in GLOBAL coordinate space and must be
        # transformed to FRAGMENT-LOCAL indices for each fragment.
        
        # R geometries
        fragA_coords_R, fragA_symbols_R = self._h_cap_geometry(
            reactant_coords, reactant_symbols, fragA_indices, cut_bond
        )
        fragB_coords_R, fragB_symbols_R = self._h_cap_geometry(
            reactant_coords, reactant_symbols, fragB_indices, cut_bond
        )
        
        # TS geometries
        fragA_coords_TS, fragA_symbols_TS = self._h_cap_geometry(
            ts_coords, ts_symbols, fragA_indices, cut_bond
        )
        fragB_coords_TS, fragB_symbols_TS = self._h_cap_geometry(
            ts_coords, ts_symbols, fragB_indices, cut_bond
        )
        
        return FragmenterResult(
            fragA_indices=fragA_indices,
            fragB_indices=fragB_indices,
            cut_bond_indices=cut_bond,
            dipole_end_indices=dipole_pair,
            alkene_end_indices=alkene_pair,
            dipole_core_indices=dipole_core,
            fragA_coords_R=fragA_coords_R,
            fragB_coords_R=fragB_coords_R,
            fragA_coords_TS=fragA_coords_TS,
            fragB_coords_TS=fragB_coords_TS,
            fragA_symbols_R=fragA_symbols_R,
            fragB_symbols_R=fragB_symbols_R,
            fragA_symbols_TS=fragA_symbols_TS,
            fragB_symbols_TS=fragB_symbols_TS,
            status="ok",
            reason=""
        )
    
    def _build_graph(self, coords, symbols, scale, min_dist):
        """Step A: Build topology graph using bond graph utilities."""
        return build_bond_graph(coords, symbols, scale=scale, min_dist=min_dist)
    
    def _classify_reaction_ends(self, graph, coords, forming_bonds):
        """
        Step B: Identify 4 reaction atoms and classify alkene vs dipole ends.
        
        Returns:
            (alkene_pair, dipole_pair) where each is (i, j) tuple
            Returns (None, None) if classification fails
        """
        # Flatten forming bonds
        reaction_atoms = set()
        for bond in forming_bonds:
            reaction_atoms.add(bond[0])
            reaction_atoms.add(bond[1])
        reaction_atoms = sorted(reaction_atoms)
        
        # Find atom pairs with dist=1 in graph (likely C=C)
        alkene_candidates = []
        for i in range(len(reaction_atoms)):
            for j in range(i + 1, len(reaction_atoms)):
                a, b = reaction_atoms[i], reaction_atoms[j]
                path = find_shortest_path(graph, a, b)
                if len(path) == 2:  # Direct bond (dist=1)
                    # Check geometric distance (typical C=C ~1.34 Å)
                    dist = GeometryUtils.calculate_distance(coords, a, b)
                    if dist < 1.55:  # Reasonable C=C upper bound
                        alkene_candidates.append(((a, b), dist))
        
        if not alkene_candidates:
            return None, None
        
        # Choose shortest bond
        alkene_pair = min(alkene_candidates, key=lambda x: x[1])[0]
        
        # Remaining atoms are dipole ends
        all_atoms = set(reaction_atoms)
        dipole_atoms = all_atoms - set(alkene_pair)
        dipole_pair = tuple(sorted(dipole_atoms))
        
        return alkene_pair, dipole_pair
    
    def _find_dipole_core_path(self, graph, dipole_pair):
        """
        Step C: Find 5-atom dipole core path.
        
        Returns:
            List of 5 atom indices (dipole_core_indices)
            Returns None if path length != 4 edges (5 atoms)
        """
        d1, d2 = dipole_pair
        path = find_shortest_path(graph, d1, d2)
        
        if len(path) != 5:
            return None  # Not 5 atoms
        
        return path
    
    def _find_cut_bond(self, graph, dipole_core, alkene_pair, symbols):
        """
        Step D: Find tether-dipole cut bond.
        
        Returns:
            (i, j) cut bond indices
            Returns None if cut would break a ring
        """
        a, b = alkene_pair
        
        # Find shortest paths from dipole core to each alkene end
        path_to_a = self._shortest_path_to_set(graph, dipole_core, a)
        path_to_b = self._shortest_path_to_set(graph, dipole_core, b)
        
        # Choose shorter path
        if len(path_to_a) <= len(path_to_b):
            tether_path = path_to_a
        else:
            tether_path = path_to_b
        
        if len(tether_path) < 2:
            return None  # Should not happen
        
        # Cut bond is first edge leaving dipole core
        p0, p1 = tether_path[0], tether_path[1]
        
        # Check if cut bond is within a ring (weak check: both ends have degree >= 2)
        # If this is too weak, implement SSSR detection (more complex)
        degree_p0 = len(graph[p0])
        degree_p1 = len(graph[p1])
        
        if degree_p0 >= 2 and degree_p1 >= 2:
            # Likely ring bond - try next edge
            if len(tether_path) >= 3:
                p1, p2 = tether_path[1], tether_path[2]
                if not (len(graph[p1]) >= 2 and len(graph[p2]) >= 2):
                    return (p1, p2)
        
        return (p0, p1)
    
    def _shortest_path_to_set(self, graph, source_set, target):
        """Find shortest path from any atom in source_set to target."""
        best_path = None
        
        for source in source_set:
            path = find_shortest_path(graph, source, target)
            if path:
                if best_path is None or len(path) < len(best_path):
                    best_path = path
        
        return best_path if best_path else []
    
    def _cut_and_get_components(self, graph, cut_bond, dipole_pair, alkene_pair):
        """
        Step E: Remove cut bond and extract components.
        
        Returns:
            (fragA_indices, fragB_indices)
            fragA must contain dipole_pair, fragB must contain alkene_pair
        """
        # Create copy of graph
        cut_graph = {k: v.copy() for k, v in graph.items()}
        
        # Remove cut bond
        i, j = cut_bond
        cut_graph[i].remove(j)
        cut_graph[j].remove(i)
        
        # Get components
        components = get_connected_components(cut_graph)
        
        if len(components) != 2:
            return [], []
        
        # Label components by which reaction atoms they contain
        comp0_set = set(components[0])
        comp1_set = set(components[1])
        
        dipole_set = set(dipole_pair)
        alkene_set = set(alkene_pair)
        
        # Component containing dipole is fragA
        if dipole_set.issubset(comp0_set):
            fragA_indices, fragB_indices = components[0], components[1]
        else:
            fragA_indices, fragB_indices = components[1], components[0]
        
        # Validation
        if not dipole_set.issubset(set(fragA_indices)):
            return [], []  # Validation failed
        if not alkene_set.issubset(set(fragB_indices)):
            return [], []
        
        return fragA_indices, fragB_indices
    
    def _compute_cap_positions(self, coords, cut_bond):
        """
        Compute H-cap positions for cut bond.
        
        Returns:
            List of [(atom_idx, direction_vec)] for each cut endpoint
        """
        i, j = cut_bond
        vec_ij = coords[j] - coords[i]
        vec_ji = coords[i] - coords[j]
        
        return [
            (i, vec_ij),  # Cap at i, pointing toward j (so we go opposite)
            (j, vec_ji)   # Cap at j, pointing toward i (so we go opposite)
        ]
```

    def _h_cap_geometry(self, full_coords, full_symbols, frag_indices, cut_bond, cap_rules=None):
        """
        H-cap a fragment geometry with correct index transformation.
        
        Args:
            full_coords: Full system coordinates (N, 3)
            full_symbols: Full system element symbols
            frag_indices: Fragment atom indices (global indices)
            cut_bond: (i, j) global indices of cut bond
            cap_rules: Optional bond length by element
        
        Returns:
            (capped_coords, capped_symbols) - Capped fragment geometry
            
        Algorithm:
            1. Slice full coordinates to fragment-local indices (0 to len(frag)-1)
            2. Map cut_bond endpoints from global to fragment-local space
            3. Compute cap positions in fragment-local space
            4. Call h_cap_fragment with correct arguments
        """
        if cap_rules is None:
            cap_rules = {"C": 1.09, "N": 1.01, "O": 0.96, "H": 0.76}
        
        # Slice to fragment-local coordinates
        frag_coords = full_coords[frag_indices]
        frag_symbols = [full_symbols[i] for i in frag_indices]
        
        # Create global-to-local index map
        global_to_local = {g: l for l, g in enumerate(frag_indices)}
        
        # Find cut bond endpoints in fragment-local space
        i_global, j_global = cut_bond
        
        # Only cap the endpoint that's actually IN this fragment
        cap_positions = []
        if i_global in global_to_local:
            i_local = global_to_local[i_global]
            # Direction: from i toward j in global space
            cap_dir = full_coords[j_global] - full_coords[i_global]
            cap_positions.append((i_local, cap_dir))
        if j_global in global_to_local:
            j_local = global_to_local[j_global]
            # Direction: from j toward i in global space
            cap_dir = full_coords[i_global] - full_coords[j_global]
            cap_positions.append((j_local, cap_dir))
        
        # Call existing h_cap_fragment with fragment-local data
        capped_coords, capped_symbols = h_cap_fragment(
            frag_coords, frag_symbols, cap_positions, cap_bond_lengths=cap_rules
        )
        
        return capped_coords, capped_symbols


**Acceptance Criteria**:
- Unit tests for each step of fragmenter.
- Test with known intramolecular [5+2] oxidopyrylium example produces valid fragments.
- Validation tests verify dipole contains dipole ends, fragment B contains alkene ends.

---

### 8. V6.2: Add S3 Post-QC Enrichment (DIAS/ASM with Tether-Cut)

**What to do**:
- Implement `S3_PostQCEnrichment/` stage that runs after TS success and writes:
   - `<S3_DIR>/S3_PostQCEnrichment/enrichment.json` (with fragmenter block, charges, multiplicities)
   - `<S3_DIR>/S3_PostQCEnrichment/enrichment_status.json`
- Enrichment must be executed from S3 and must reuse `qc_interface`/`orca_interface` or `QCTaskRunner` (no ad-hoc subprocess or cwd hacks).

**S3 directory resolution (for S4 context injection)**:
- Preferred: `work_dir/S3_TransitionAnalysis/` (current orchestrator).
- Fallback: `work_dir/S3_TS/` (legacy runs).
- Do NOT do global tree scans; only check these known step directories.
- Implementation: S4 context must resolve `s3_dir` consistently using this preference order (matching `mech_packager` resolution pattern if needed).
- Integration point and failure policy (explicit):
  - Invoke enrichment from `rph_core/steps/step3_opt/ts_optimizer.py:TSOptimizer.run_with_qctaskrunner()` after:
    - TS optimization is converged,
    - `ts_final.xyz` has been written to `<S3_DIR>/ts_final.xyz`, and
    - `SPMatrixReport` is constructed (so `e_ts_final` is available for `ts_TS`).
  - If enrichment fails:
    - Step3 must still succeed overall.
    - Write `<S3_DIR>/S3_PostQCEnrichment/enrichment_status.json` with status=failed and an error message.
    - Step4's asm plugin must treat this as `missing_reason=enrichment_failed` and output NaNs.
- Implementation placement (explicit, avoid bloating ts_optimizer.py):
  - Add new module: `rph_core/steps/step3_opt/post_qc_enrichment.py`
  - Expose a single entrypoint function, e.g. `run_post_qc_enrichment(s3_dir: Path, config: dict, sp_report: SPMatrixReport, reactant_complex_xyz: Path, forming_bonds: Tuple[Tuple[int, int], Tuple[int, int]]) -> None`.
  - `ts_optimizer.py` should only call this function and handle exceptions.

**Fragment definition strategy (REPLACED with tether-cut)**:
- Use **intramolecular tether-cut** (see TODO 7 for implementation details).
- Do NOT use "connectivity split" (direct components) - this fails for intramolecular systems.
- Guardrails:
  - If fragmenter returns `status=failed`, write `enrichment_status.json` with the reason and skip enrichment (no ORCA runs).
  - Record all fragmenter metadata (cut_bond, dipole_core, etc.) into the contract JSON.
- Override escape hatch (MUST HAVE for edge cases):
  - If `step3.enrichment.fragment_indices_override` is provided, bypass tether-cut entirely.
  - Record `fragment_source=override` in the contract.

**Configuration surface (explicit defaults)**:
- Config file location for defaults: `config/defaults.yaml`.
- `step3.enrichment.enabled` (default false)
- `step3.enrichment.force_rerun` (default false)
- `step3.enrichment.connectivity_scale` (default 1.25)
- `step3.enrichment.bond_min_dist_angstrom` (default 0.6)
- `step3.enrichment.write_dirname` (default `S3_PostQCEnrichment`)
- `step3.enrichment.fragment_indices_override` (default null; when set, bypass tether-cut)
- `step3.enrichment.compute_ts_ts_orca` (default false; if true, run TS@TS ORCA SP for interaction closure)
- `step3.enrichment.fragmenter.type` (default `intramolecular_tether_cut`)

**ORCA template identity definition (REQUIRED for caching/provenance)**:
- `orca_template_id`: Canonical string uniquely identifying the ORCA settings used for enrichment.
  - Format: `"orca:{method}/{basis}:{aux_basis}:{solvent}:{nprocs}:{maxcore}"`
  - Example: `"orca:M062X/def2-TZVPP:def2/J:acetone:16:8000"`
  - This string allows two developers to reproduce the exact ORCA configuration without guessing.
- `orca_template_hash`: SHA256 hash of the rendered ORCA `.inp` file content.
  - Hash the ENTIRE input file (including charge/multiplicity, coordinates, keywords).
  - This ensures any change to the ORCA method/basis/solvent/keywords invalidates the cache.

/tmp/insert_system_charge.txt
- **For `.out` files**: Line 138-141 shows it attempts to parse charge/spin via `CoordinateExtractor.get_charge_spin_from_orca_out()`.
- **For `.xyz` files**: No explicit charge/spin parsing logic exists. Default is charge=0, spin=1.
  - In generated `.inp` file (line 62), these defaults are embedded: `* xyz {charge} {spin}`.
  - **IMPORTANT**: The `spin` value in `.inp` corresponds to multiplicity `S` (singlet → spin=0, doublet → spin=1, triplet → spin=2).
  - Formula: `spin = 2*S - 1` where S is multiplicity.
- **Implementation approach for enrichment**:
  - Option A (Recommended): Extend `ORCAInterface.single_point()` to accept explicit `charge` and `multiplicity` parameters.
    - Modify `ORCAInterface._generate_input()` to accept `charge=None, multiplicity=None` parameters.
    - If charge/multiplicity are provided, use them explicitly in the generated `.inp` file: `* xyz {charge} {spin}` where `spin = 2*mult - 1`.
    - If not provided, fall back to current behavior (defaults=charge=0, spin=1 or parsed from `.out`).
    - Add new wrapper function in enrichment: `single_point_with_charge(xyz, output_dir, charge, multiplicity)`.
  - Option B (Alternative): Have `IntramolecularFragmenter` write XYZ files with ORCA-compatible charge/spin headers.
    - Append comment line to fragment XYZ files: e.g., `charge=2\nspin=1\n` on line after atom count.
    - ORCAInterface already reads this via `CoordinateExtractor.get_charge_spin_from_orca_out()` (line 141-145).
    - This requires modifying `h_cap_fragment()` to write headers.
    - **Decision**: Option A is preferred. Option B is acceptable fallback but requires documenting exact header format.
- **Fragment charge distribution**:
  - Full system charge from `SPMatrixReport` is used to assign fragment charges.
  - Default policy (configurable): assign full system charge to Fragment A (dipole), Fragment B gets charge 0.
  - Oxidopyrylium systems typically have charge=+1, so fragA=+1, fragB=0.
  - Override: `step3.enrichment.fragment_charges_override` can provide explicit `{fragA: int, fragB: int}`.
- **Fragment multiplicities**:
  - Default: 1 (singlet) for both fragments.
  - Override: `step3.enrichment.fragment_multiplicities_override` can provide explicit values.
  - Spin parameter in ORCA: `spin = 2*S - 1` where S is multiplicity (S=1 → spin=0).

**Enrichment contract schema (REQUIRED, because caching depends on it)**:
- `<S3_DIR>/S3_PostQCEnrichment/enrichment.json` must contain (top-level required keys):
  - `schema_version`: string (e.g., `"enrichment_v2"`)
  - `created_at`: ISO timestamp
  - `orca_template_id`: string
  - `orca_template_hash`: sha256 hex
  - `geometry_hashes`:
    - `reactant_complex_xyz_sha256`: sha256 of raw bytes of `S2_Retro/reactant_complex.xyz`
    - `ts_final_xyz_sha256`: sha256 of raw bytes of `<S3_DIR>/ts_final.xyz`
  - `fragmenter`: (NEW BLOCK)
    - `type`: `intramolecular_tether_cut|override`
    - `version`: `1.0`
    - `cut_bond_indices`: [i, j]
    - `fragA_indices`: sorted int list
    - `fragB_indices`: sorted int list
    - `fragment_source`: `override|tether_cut`
    - `fragment_hash`: sha256 of JSON({fragA_indices, fragB_indices, cut_bond_indices}) with sort_keys
    - `dipole_end_indices`: [d1, d2]
    - `alkene_end_indices`: [a, b]
    - `dipole_core_indices`: [5 atom indices]
    - `cap_rules`: {"C":1.09,"N":1.01,"O":0.96,...}
    - `cap_method`: `colinear_back_projection`
  - `fragment_charges`: (NEW BLOCK)
    - `fragA`: charge for fragment A (default +1 for oxidopyrylium)
    - `fragB`: charge for fragment B (default 0)
  - `fragment_multiplicities`: (NEW BLOCK)
    - `fragA`: multiplicity for fragment A (default 1)
    - `fragB`: multiplicity for fragment B (default 1)
  - `energies_hartree`:
    - `fragA_R`, `fragB_R`, `fragA_TS`, `fragB_TS`
    - optional `ts_TS` (ONLY when `step3.enrichment.compute_ts_ts_orca=true`; must be computed with the same ORCA template as fragments)
  - `units`: fixed string `{"energies_hartree": "hartree"}`

- `<S3_DIR>/S3_PostQCEnrichment/enrichment_status.json` must contain:
  - `schema_version`: string (e.g., `"enrichment_status_v1"`)
  - `created_at`: ISO timestamp
  - `cache_key`: sha256 of JSON({orca_template_hash, geometry_hashes, fragment_hash}) with sort_keys (fragment_hash now includes cut_bond and cap_rules)
  - `cache_hit`: boolean
  - `status`: `ok|skipped|failed`
  - `error`: optional string

**Enrichment failure/skip semantics (explicit)**:
- **When fragmenter fails** (dipole_path_not_5_atoms, failed_to_find_cut_bond, component_validation_failed):
  - Write `enrichment_status.json` with `status=failed`.
  - Set `error` to the exact fragmenter `reason` (e.g., `"dipole_path_not_5_atoms"`).
  - Do NOT run any ORCA calculations.
  - S4 reads this and outputs `asm.distortion_total_kcal=NaN` with `missing_reason="enrichment_failed"`.
- **When forming_bonds is missing/invalid** (detected before enrichment):
  - Write `enrichment_status.json` with `status=skipped`.
  - Set `error="forming_bonds_invalid"` or `"forming_bonds_missing"`.
  - Do NOT run fragmenter or ORCA.
  - S4 asm_plugin treats this as `missing_reason="enrichment_skipped"` and outputs NaN.
- **When enrichment succeeds** (fragmenter ok + all ORCA SPs converged):
  - Write `enrichment_status.json` with `status=ok`.
  - `error` field is empty/null.
- **When ORCA SP fails** (fragmenter ok but ORCA converges with errors):
  - Write `enrichment_status.json` with `status=failed`.
  - Set `error` to the ORCA error message (e.g., `"ORCA SCF not converged"`).
  - Partial energies may be written (e.g., only 2 of 4 fragments succeeded) but `status=failed`.
- **When cache is hit** (valid existing enrichment.json with same cache_key and `force_rerun=false`):
  - Write `enrichment_status.json` with `status=ok`, `cache_hit=true`.
  - No new ORCA runs.
  - `created_at` timestamp updates.

**ASM/DIAS energy semantics (no ambiguity)**:
- Define geometries for fragment single points:
  - `frag*_R`: capped fragment cut from `S2_Retro/reactant_complex.xyz` using inferred atom indices (no relaxation).
  - `frag*_TS`: capped fragment cut from `<S3_DIR>/ts_final.xyz` (Step3 canonical TS) using the *same* atom indices (no relaxation).
- Define full-system TS energy used for interaction closure:
  - `ts_TS` is ONLY defined when `step3.enrichment.compute_ts_ts_orca=true` and is computed with the same ORCA template as the fragment SPs.
  - Otherwise omit `ts_TS` and omit `asm.interaction_kcal` downstream.
- Derived quantities (S4 consumption):
  - `asm.distortion_total_kcal = Σ[(E_frag_TS − E_frag_R) * 627.509]`
  - `asm.interaction_kcal = (E_ts_TS − Σ(E_frag_TS)) * 627.509` (only if `ts_TS` is present and explicitly marked as Hartree).

**Toxic-path-safe ORCA execution (explicit)**:
- Implementation choice (pin down): extend `rph_core/utils/orca_interface.py` so ORCA execution uses `rph_core/utils/qc_interface.py:is_path_toxic()` (spaces + `[](){} `), not just whitespace.
- **Base directory choice**: Use a known-safe base dir that works across platforms (e.g., `/tmp` or `tempfile.gettempdir()`).
- In `_run_orca(...)`, implement the following logic:
  - Determine `actual_output_dir`: If `is_path_toxic(output_dir)` is true, use a temporary sandbox dir; otherwise use `output_dir` directly.
  - If sandbox is needed:
    1. Create `LinuxSandbox` directory (from `qc_interface`) with a unique name (e.g., `/tmp/orca_sandbox_<random>`).
    2. Copy the `.inp` file into the sandbox.
    3. Run ORCA with `cwd=sandbox_dir` (not toxic `output_dir`).
    4. Copy the `.out` file and any auxiliary files needed for audit (e.g., `.gbw`, `.molden`) from sandbox back to the real `output_dir`.
    5. Cleanup sandbox (delete directory).
  - If no sandbox needed, run ORCA in-place in `output_dir`.
- This ensures enrichment inherits the repo-wide toxic-path policy without adding step-local hacks, and works correctly even when run dir is toxic.

**References**:
- `rph_core/steps/step3_opt/ts_optimizer.py` (S3 lifecycle + where to hook)
- `rph_core/utils/orca_interface.py` (ORCA execution wrapper)
- `rph_core/utils/qc_interface.py` (toxic path + sandbox rules)
- `rph_core/steps/step4_features/fragment_extractor.py` (legacy fragment splitting + QC-in-step4 anti-pattern; reuse only non-QC splitting ideas)

**Acceptance Criteria**:
- With enrichment enabled, at least one intramolecular [5+2] sample produces valid `<S3_DIR>/S3_PostQCEnrichment/enrichment.json` with fragmenter block and ASM energies.
- Fragmenter correctly identifies dipole core and cut bond for oxidopyrylium systems.
- With enrichment disabled, pipeline behavior matches previous version (no new files required).
- `python -m pytest -q` passes.

---

### 9. V6.2: S4 Consumes Enrichment Contract (No Compute)

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
- Intramolecular [5+2] samples with valid enrichment produce non-NaN ASM features.
- Add fixture-based unit tests (no QC required):
  - `tests/fixtures/enrichment/enrichment_tether_cut.json` (minimal valid contract with fragmenter block)
  - Test asserts `asm.distortion_total_kcal` is computed exactly.
  - Separate test asserts `asm.interaction_kcal` is computed exactly when `ts_TS` is present.

---

### 10. V6.3: Add Dipolar FMO/CDFT Parser (S4 Text Parse Only)

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
    - Filter candidates to those with "normal termination" markers (Gaussian: `Normal termination`; ORCA: `ORCA TERMINATED NORMALLY`).
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

**Test baseline repair (REQUIRED before refactor)**:
- Current test suite has structural issues (e.g., syntax/indentation errors in `tests/test_m2_schema_versioning.py`).
- Add explicit task: Restore test suite to runnable baseline.
  - Fix any Python syntax/indentation errors.
  - Ensure `python -m pytest -q` passes on clean repo state.
  - Document that test repairs are part of V6.1 (not blocking V6.2+).
- This prevents confusion between "existing bugs" and "bugs introduced by refactor".

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

### 11. V6.4: Hardening (Caching/Resume/Strict QC Mode)

**What to do**:
- S3 enrichment caching:
  - Derive cache key from `orca_template_hash + geometry_hash + fragment_hash` (fragment_hash now includes cut_bond and cap_rules from enrichment.json).
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
- Cache key correctly includes fragmenter cut_bond and cap_rules (changing these invalidates cache).

---

### 12. V6.4: Modeling Companion Script (LOOCV/VIF/LFER; not pipeline)

**What to do**:
- Add `scripts/ml/train_mlr_loocv.py`:
  - Load `features_mlr.csv`.
  - Run LOOCV.
  - Compute VIF (default threshold 5) and/or perform Lasso-based feature selection.
  - Optional fold-internal LFER: fit `ΔG ≈ a·ΔE + C` on train fold and impute missing ΔG in test fold without leakage.
- Keep this script out of step execution; it's a post-analysis utility only.

**References**:
- `S4_Data/features_mlr.csv` (new contract)

**Acceptance Criteria**:
- Script runs on a small sample directory and emits a report (stdout + optional CSV/JSON in `scripts/ml/out/`).

---

### 13. Self Gap Check (Contracts & No Strays)

**What to do**:
- Run a repo-wide check (excluding `rph_core_backup_20260115/`) to ensure no runtime code or tests still require `features.csv`.
- Verify new contracts are referenced consistently:
  - `features_raw.csv` / `features_mlr.csv` / `feature_meta.json` in Step4 code and tests.
  - `<S3_DIR>/S3_PostQCEnrichment/enrichment.json` in Step3/S4 linkage.
- **NEW**: Verify "connectivity split" is NOT mentioned in runtime code (replaced by "intramolecular_tether_cut").

**Acceptance Criteria**:
- `python -m pytest -q` passes.
- Grep for `features.csv` outside backup/notes is empty.
- Grep for `connectivity split` in runtime code (excluding tests/docs) returns no matches.

---

### 14. Docs & Contract Acceptance

**What to do**:
- Add/Update documentation (suggested): `V6_ACCEPTANCE.md`:
  - New S4 output contract and meaning.
  - Enrichment contract schema (including fragmenter block, charges, multiplicities).
  - Dipolar plugin schema.
  - Recommended default mlr columns.
  - Verification commands.
  - **NEW**: Intramolecular tether-cut fragmenter algorithm explanation.
  - **NEW**: Why "connectivity split" fails for intramolecular systems.

**Acceptance Criteria**:
- Doc clearly states IN/OUT, contracts, and reproducibility expectations.
- Fragmenter algorithm is clearly documented with physical rationale.

---

## Appendix: `feature_meta.json` Example (V6.1)

```json
{
  "meta": {
    "schema_version": "6.1",
    "schema_signature": "<sha1>",
    "feature_status": "ok",
    "method": "<method>",
    "solvent": "<solvent>",
    "temperature_K": 298.15,
    "enabled_plugins": ["thermo", "geometry", "qc_checks", "ts_quality", "asm_enrichment"],
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
      "asm.distortion_total_kcal": {"is_missing": false, "is_invalid": false, "reason": ""},
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
      "asm_enrichment": {
        "status": "OK",
        "runtime_ms": 45.6,
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

## Appendix: `enrichment.json` Example (V6.2 Tether-Cut)

```json
{
  "schema_version": "enrichment_v2",
  "created_at": "2026-01-27T12:34:56Z",
  "orca_template_id": "orca:B3LYP/def2-SVP:::scrf(acetone):16:8000",
  "orca_template_hash": "a1b2c3d4e5f6...",
  "geometry_hashes": {
    "reactant_complex_xyz_sha256": "abc123...",
    "ts_final_xyz_sha256": "def456..."
  },
  "fragmenter": {
    "type": "intramolecular_tether_cut",
    "version": "1.0",
    "cut_bond_indices": [5, 12],
    "fragA_indices": [0, 1, 2, 3, 4, 5],
    "fragB_indices": [6, 7, 8, 9, 10, 11, 12, 13, 14],
    "fragment_source": "tether_cut",
    "fragment_hash": "789abc...",
    "dipole_end_indices": [0, 4],
    "alkene_end_indices": [8, 9],
    "dipole_core_indices": [0, 1, 2, 3, 4],
    "cap_rules": {"C": 1.09, "N": 1.01, "O": 0.96},
    "cap_method": "colinear_back_projection"
  },
  "fragment_charges": {
    "fragA": 1,
    "fragB": 0
  },
  "fragment_multiplicities": {
    "fragA": 1,
    "fragB": 1
  },
  "energies_hartree": {
    "fragA_R": -384.567890,
    "fragB_R": -234.567890,
    "fragA_TS": -384.523456,
    "fragB_TS": -234.534567
  },
  "units": {
    "energies_hartree": "hartree"
  }
}
```

---

## Defaults Applied (override if needed)
- Keep `features_mlr.csv` default columns ≤10, starting from:
  - `thermo.dE_activation`, `thermo.dE_reaction`, `geom.r_avg`, `geom.dr`, `geom.close_contacts_density`, `ts.imag1_cm1_abs`
  - plus at most one of: `asm.distortion_total_kcal` (if V6.2 enabled) or `fmo.dipolar_omega_ev`/`fmo.dipolar_gap_ev` (if V6.3 enabled)
- FMO dipolar plugin default: OFF initially.
- Fragmenter type: `intramolecular_tether_cut` (NOT `connectivity_split`).

## Decisions Needed
- None (blocking decisions resolved).
