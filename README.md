
 ReactionProfileHunter v6.1.0
ReactionProfileHunter (RPH) is a product-first reaction mechanistic pipeline:
S1 (anchor/conformer) -> S2 (retro scan) -> S3 (TS optimize/rescue) -> S4 (feature extraction + packaging).
V6.1 goal: S3 performs all computation; S4 is extract-only and writes fixed 3-file outputs.
 Key Changes In v6.1
- S4 outputs: `S4_Data/features_raw.csv`, `S4_Data/features_mlr.csv`, `S4_Data/feature_meta.json`.
- S4 does not submit QC jobs; only parses existing outputs.
- NBO compute and artifact collection are disabled by default.
 Entrypoints
- CLI wrapper: `bin/rph_run`
- Python CLI: `rph_core/orchestrator.py` (`main()`)
 Install
pip install -e ".[dev]"
Run
Single (override config.run.single.product_smiles):
bin/rph_run --smiles "C=C(C)C(=O)O" --output ./Output/rx_manual
Config-driven run (recommended): edit config/defaults.yaml then:
bin/rph_run
Architecture (Code Map)
- Orchestrator / wiring: rph_core/orchestrator.py
- Step1 anchor/conformer:
  - rph_core/steps/anchor/
  - rph_core/steps/conformer_search/
- Step2 retro scan: rph_core/steps/step2_retro/retro_scanner.py
- Step3 TS optimize/rescue: rph_core/steps/step3_opt/ts_optimizer.py
- Step4 feature extraction + packaging:
  - rph_core/steps/step4_features/feature_miner.py
  - rph_core/steps/step4_features/mech_packager.py
- QC sandbox execution: rph_core/utils/qc_interface.py
- OPT/Freq driver and route assembly: rph_core/utils/qc_task_runner.py
- Resume/checkpoint: rph_core/utils/checkpoint_manager.py
IO Contract (Required Outputs)
RPH is built around a strict step output contract:
S1 output (required):
- S1_ConfGeneration/product/product_min.xyz
- S1_ConfGeneration/precursor/precursor_min.xyz (optional)
S2 output (required; must output BOTH):
- S2_Retro/ts_guess.xyz
- S2_Retro/reactant_complex.xyz
S3 output (required):
- S3_TS/ts_final.xyz
- S3_TS/reactant_sp.xyz (canonical reactant geometry for downstream)
- S3_TS/reactant_opt/standard/ or S3_TS/reactant_opt/rescue/ (reactant OPT+Freq run directory)
S4 output (required):
- S4_Data/features_raw.csv
- S4_Data/mech_index.json
S4 QC artifacts (optional; v5.4 only NBO):
- S4_Data/qc_nbo.37 (only if a candidate NBO file is found in S3)
NBO Integration (v5.4)
- Reactant OPT+Freq is executed via QCTaskRunner.run_opt_sp_cycle(..., enable_nbo=True):
  - route gets `Pop=NBO` appended by default (see rph_core/utils/qc_task_runner.py).
  - if a keylist is provided, the runner uses `Pop=(NBORead)`.
- No separate SP-NBO task is submitted.
- NBO files are collected in Step4 from S3 whitelisted subdirs:
  - S3_TS/nbo_analysis/
  - S3_TS/nbo/
  - S3_TS/reactant_opt/standard/
  - S3_TS/reactant_opt/rescue/
  (see QC_ARTIFACT_SUBDIRS in rph_core/steps/step4_features/mech_packager.py)
- Packager copies the picked file to S4_Data/qc_nbo.37 and records qc_artifacts.nbo_outputs in mech_index.json.
Removed In v5.4
- Hirshfeld charge parsing/collection
- NMR GIAO calculation/collection
- Gaussian templates for NMR/Hirshfeld
- Any optional QC task runner that previously generated NMR/Hirshfeld outputs
Output Directory Layout (Example)
Output/rx_xxx/
  S1_ConfGeneration/
    product/
      product_min.xyz
    precursor/
      precursor_min.xyz
  S2_Retro/
    ts_guess.xyz
    reactant_complex.xyz
  S3_TS/
    ts_final.xyz
    reactant_sp.xyz
    reactant_opt/
      standard/  (reactant OPT+Freq + NBO)
      rescue/    (rescue OPT+Freq + NBO)
  S4_Data/
    features_raw.csv
    mech_index.json
    qc_nbo.37      (optional)
  rph.log
Testing
python -m pytest -q
Note: tests/conftest.py adds repo root to sys.path so tests can import rph_core without requiring editable install.
## `PLAN.md` (v5.4) 建议完整内容
```markdown
# ReactionProfileHunter v5.4 Plan (Architecture + IO + Feature Spec)
Document date: 2026-01-23
Target version: 5.4.0
Goal:
- Integrate NBO into reactant OPT+Freq.
- Remove Hirshfeld and NMR end-to-end (compute + packaging + config + templates + docs/tests).
## Scope
In scope (v5.4):
- NBO is produced during reactant OPT+Freq and packaged into S4.
- Step4 qc_artifacts becomes NBO-only.
Out of scope (v5.4):
- NEB/IRC workflow expansion
- Additional QC engines beyond current gaussian/orca/xtb integration
- ML training/UI
## Current Pipeline
Serial 4-step pipeline (product-first):
- S1 Anchor/Conformer -> S2 Retro scan -> S3 TS optimize/rescue -> S4 Feature extraction + packager
Key wiring:
- `rph_core/orchestrator.py` runs the pipeline and enforces step dependencies.
## Output Contract (Hard Requirements)
- S2 MUST output both `ts_guess.xyz` and `reactant_complex.xyz`.
- S3 MUST output `ts_final.xyz`.
- S4 MUST output `features_raw.csv` and `mech_index.json`.
Locations:
- Contract summary: `rph_core/steps/AGENTS.md`
## V5.4 Implementation Notes
### 1) NBO integrated into reactant OPT+Freq
Where enabled:
- `rph_core/steps/step3_opt/ts_optimizer.py` calls:
  - `QCTaskRunner.run_opt_sp_cycle(..., enable_nbo=True)` for the reactant branch only.
How implemented:
- `rph_core/utils/qc_task_runner.py`
  - `run_opt_sp_cycle(..., enable_nbo: bool = False)`
  - `_try_normal_optimization()` / `_try_normal_rescue()` append `Pop=NBO` by default when enabled.
Expected artifacts:
- In S3, NBO files can appear as `*.37`, `*.nbo`, `*.nbo7` depending on backend.
- S4 package uses fixed name `qc_nbo.37`.
### 2) Remove Hirshfeld and NMR end-to-end
Removed surfaces (v5.4):
- config/templates:
  - remove `config/templates/gaussian_nmr.gjf`
  - remove `config/templates/gaussian_hirshfeld.gjf`
- config:
  - remove `qc_tasks` section from `config/defaults.yaml`
- QC interface:
  - `rph_core/utils/qc_interface.py` no longer exposes TaskKind entries or helpers for NMR/Hirshfeld.
- Packager:
  - `rph_core/steps/step4_features/mech_packager.py` qc_artifacts constants are NBO-only.
### 3) Step4 qc_artifacts (NBO-only)
- Whitelisted search roots (relative to S3 root):
  - `nbo_analysis/`, `nbo/`, `reactant_opt/standard/`, `reactant_opt/rescue/`
- Patterns:
  - `*.37`, `*.nbo`, `*.nbo7`
- Target:
  - `qc_nbo.37`
- Result recorded in `mech_index.json` under `qc_artifacts.nbo_outputs`.
Implementation:
- `rph_core/steps/step4_features/mech_packager.py`
## Tests And Docs Policy (v5.4)
- No test should reference NMR/Hirshfeld templates or outputs.
- Mock tests should only assert NBO artifact behavior.
- `tests/conftest.py` ensures `import rph_core` works without editable install.
## Acceptance Checklist (v5.4)
- `python -m compileall rph_core` passes
- `python -m pytest -q` passes
- Grep residue scan:
  - no `NMR_GIAO`, `NMR=GIAO`, `hirshfeld`, `gaussian_nmr.gjf`, `gaussian_hirshfeld.gjf` in runtime code/tests/docs (except historical docs if explicitly kept)
## Known Docs To Update
- `README.md` should reflect v5.4 changes (NBO-only qc_artifacts; remove NMR/Hirshfeld).
- Legacy docs like `V5_ACCEPTANCE.md` may still mention NMR/Hirshfeld and should be revised if kept as authoritative.
