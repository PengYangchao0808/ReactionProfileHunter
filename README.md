# ReactionProfileHunter

[English](README.md) | [中文](README.zh-CN.md)

<div align="center">

**Product-driven DFT reaction mechanism pipeline (S0-S4)**

[![Version](https://img.shields.io/badge/version-4.0.0-blue.svg)](https://github.com/PengYangchao0808/ReactionProfileHunter)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Code Style](https://img.shields.io/badge/code%20style-agents%20ready-success.svg)](AGENTS.md)

</div>

> **Project Stats:** ~50k lines, 198 Python files | **Agent Guide:** See [AGENTS.md](AGENTS.md) for coding standards

---

## Overview

ReactionProfileHunter (RPH) is a product-driven automated DFT reaction mechanism pipeline (S0-S4) for transition state search and geometry optimization. Feature extraction and ML training are handled externally by the companion `RPH_Postprocess` package.

### Core capabilities

- Product-driven strategy: start from products and search reaction pathways backward
- Five-stage DFT pipeline: mechanism validation -> CENSO-LITE conformer search -> PEB retro scan -> low-level QC -> high-precision QC
- Fixed theory contract: CREST/GFN2 + B97-3c SP (S1), ORCA B97-3c -> r2SCAN-3c (S3), Gaussian M062X -> ORCA wB97M-V (S4)
- Multi-engine support: Gaussian, ORCA, xTB, CREST
- Manifest-based checkpoint/resume with per-structure failure isolation
- WSL live status viewer (`rph_watch`)
- Dataset-only entry: no SMILES mode

### Use cases

- Mechanistic studies of organic reactions (cycloadditions, rearrangements, substitutions)
- Transition state prediction and validation
- Reaction dataset building and high-throughput screening
- Machine-learning dataset generation (feature extraction via RPH_Postprocess)

---

## Key features

### V4 architecture — fixed S0-S4 pipeline

The V4 pipeline replaces the V3 Berny/QST2/IRC rescue chain and SMILES entry with a clean, reproducible protocol:

| Stage | Protocol | Detail |
|-------|----------|--------|
| S0 | Trusted dataset mechanism validation | Reads Reaxys-cleaned CSV; resolves forming bonds, atom mapping, reaction type |
| S1 | CENSO-LITE conformer search | CREST/GFN2 sampling + ORCA B97-3c SP ranking + xTB mRRHO. No DFT OPT/FREQ. |
| S2 | PEB retro scan | xTB potential-energy-bond scan from S1 selected conformer; produces TS guess + intermediate |
| S3 | Low-level QC | ORCA B97-3c OPT/OptTS + Freq; ORCA r2SCAN-3c SP (CPCM acetone) |
| S4 | High-precision QC | Gaussian M062X OPT/OptTS + Freq; ORCA wB97M-V SP (CPCM acetone) |

Pipeline flow:

```
CSV reaction record -> S0 (mechanism) -> S1 (CENSO-LITE) -> S2 (PEB)
                                                              |
                                                              v
                                               S3 (B97-3c OPT + r2SCAN-3c SP)
                                                              |
                                                              v
                                               S4 (M062X OPT + wB97M-V SP)
                                                              |
                                                              v
                                               RPH_Postprocess (external)
                                                 -> features, ML datasets
```

### Removed V3 features

The following V3 capabilities are not part of V4:

- SMILES-only pipeline entry (removed; CSV dataset is the only input)
- Berny / QST2 rescue / IRC validation chain
- S1 DFT OPT or DFT SP (S1 is CENSO-LITE only)
- Two-stage conformer funnel (GFN0 -> GFN2 -> DFT)
- Lewis acid calibration subsystem
- `reaction_profiles` config section
- DR aggregator and condition_feature_merger
- `QCTaskRunner` and legacy `CheckpointManager`
- S4 as feature extraction (now high-precision QC)

### Technical highlights

- Molecular autonomy: each molecule has its own working directory (S1_ConfSearch/[molecule]/)
- Per-structure failure isolation: a failed S4 job never deletes usable S3 outputs
- Manifest-based handoffs: each stage writes a versioned manifest; the next stage reads it
- Sandbox isolation: each QC task runs in an isolated sandbox
- `pipeline.state` driven resumability: re-run the same command to continue

---

## Project Structure

```
ReactionProfileHunter/
├── rph_core/              # Core source — absolute imports only
│   ├── v4_orchestrator.py # V4 pipeline orchestration + CLI (V4Orchestrator)
│   ├── __main__.py        # python -m rph_core entry point
│   ├── v4_watch.py        # WSL live status viewer
│   ├── scheduling/        # Batch scheduling & task orchestration
│   ├── steps/             # S0-S4 implementations
│   │   ├── conformer_search/  # S1 CENSO-LITE engine
│   │   ├── mechanism_classifier/  # S0 mechanism validation
│   │   ├── step2_retro/       # S2 PEB scanner
│   │   ├── step3_lowlevel/    # S3 low-level QC (B97-3c / r2SCAN-3c)
│   │   ├── step4_highlevel/   # S4 high-precision QC (M062X / wB97M-V)
│   │   └── stage_calculator.py  # Shared OPT+SP+FRFQ dispatch for S3/S4
│   └── utils/             # QC/IO/checkpoint infrastructure
│       ├── qc_interface.py     # All QC calls route here
│       ├── v4_checkpoint.py    # V4 manifest-based checkpoint
│       ├── qc_jobs.py          # QC job spec mapping
│       └── qc_models.py        # QC data models
├── bin/                  # CLI wrappers (sys.path + main())
│   ├── rph_run           # V4 primary entry point
│   ├── rph_watch         # WSL live status viewer
│   └── rph               # Short alias
├── config/               # Runtime configuration (single source of truth)
│   ├── defaults.yaml     # All config keys
│   └── templates/        # Gaussian .gjf/.com templates
├── tests/                # pytest test suite (~82 files)
│   ├── conftest.py       # Adds repo root to sys.path
│   └── fixtures/         # Static test data
├── scripts/              # Utility scripts & CI tools
│   ├── ci/check_imports.py   # Import style gate
│   └── *.py / *.sh           # Data build, analysis
├── rph_benchmark/        # Stand-alone evaluation suite
│   ├── confsearch/       # Conformer search protocol benchmarks
│   └── dft_theory/       # DFT method benchmarks
├── docs/                 # Design docs & analysis reports
├── AGENTS.md             # V4 agentic coding guide
├── README.md             # <- You are here
└── README.zh-CN.md       # Chinese version
```

---

## Quick start

### System requirements

| Component | Requirement |
|-----------|-------------|
| Operating system | Linux (Ubuntu 20.04+ or CentOS 7+ recommended) |
| Python | 3.8 or newer |
| Memory | 16 GB minimum (64 GB recommended) |
| CPU | Multi-core CPU (16+ cores recommended) |

### Quantum chemistry dependencies

V4 requires all of the following for full S0-S4 operation:

- Gaussian 16 — S4 geometry optimization
- ORCA — S3/S4 single-point energies, S1 B97-3c ranking
- xTB — S1 mRRHO thermochemistry, S2 PEB scan
- CREST — S1 conformer sampling

### Installation

1. Clone the repo
   ```bash
   git clone https://github.com/PengYangchao0808/ReactionProfileHunter.git
   cd ReactionProfileHunter
   ```

2. Install Python dependencies (if packaging metadata is available)
   ```bash
   pip install -e .
   ```

3. Configure QC executables

   Edit `config/defaults.yaml` and fill in installed paths:
   ```yaml
   executables:
     gaussian:
       path: "/path/to/g16/g16"
       root: "/path/to/g16"
       profile: "/path/to/g16/g16.profile"
     orca:
       path: "/path/to/orca/orca"
       ld_library_path: "/path/to/orca"
     xtb:
       path: "/path/to/xtb/bin/xtb"
     crest:
       path: "/path/to/crest/crest"
   ```

4. Verify the environment
   ```bash
   python -m pytest -q
   ```

### Basic usage

#### CLI entry (dataset-only)

V4 does not accept SMILES input. You must provide a trusted Reaxys-cleaned CSV:

```bash
bin/rph_run --csv data/trusted_reactions.csv --rx-id RXN_000001 --output ./Output/RXN_000001
```

#### Partial pipeline run

Use `--stop-after` to run only selected stages:

```bash
# Run S0 through S1 only (conformer search)
bin/rph_run --csv data/trusted_reactions.csv --rx-id RXN_000001 --output ./Output/RXN_000001 --stop-after s1

# Run through S3 only (skip high-precision S4)
bin/rph_run --csv data/trusted_reactions.csv --rx-id RXN_000001 --output ./Output/RXN_000001 --stop-after s3
```

#### WSL live status viewer

While a pipeline is running, monitor progress in another terminal:

```bash
bin/rph_watch --output ./Output/RXN_000001 --watch
```

#### Custom config

```bash
bin/rph_run --csv data/trusted_reactions.csv --rx-id RXN_000001 --output ./Output/RXN_000001 --config config/custom.yaml
```

#### Python API

```python
from pathlib import Path
from rph_core.v4_orchestrator import V4Orchestrator
from rph_core.steps.mechanism_classifier.s0_record import load_s0_reaction_record

record = load_s0_reaction_record(Path("data/trusted_reactions.csv"), "RXN_000001")
orchestrator = V4Orchestrator()
result = orchestrator.run(record, Path("./Output/RXN_000001"), stop_after="s4")
if result.get("success"):
    print(f"Pipeline completed. Manifests: {result}")
```

#### Entry points

| Method | Command | Notes |
|--------|---------|-------|
| CLI script | `bin/rph_run --csv <file> --rx-id <id> --output <dir>` | Adds repo root to sys.path |
| Short alias | `bin/rph --csv <file> --rx-id <id> --output <dir>` | Identical to rph_run |
| Python module | `python -m rph_core --csv <file> --rx-id <id> --output <dir>` | Same flags, relies on Python path |
| Programmatic | `V4Orchestrator().run()` | Full control, returns dict |

---

## Detailed documentation

### Architecture overview

V4 is a sequential pipeline driven by manifest-based handoffs:

```
  CSV row
     |
     v
S0_Mechanism/mechanism.json -- forming_bonds, reaction_type, mapping
     |
     v
S1_ConfSearch/product/manifest.json -- selected.xyz, candidates, energies
     |
     v
S2_PEB/manifest.json -- peb_peak (TS guess), intermediate_seed, scan_profile
     |
     v
S3_LowLevel/manifest.json -- per-structure opt, freq, sp results
     |
     v
S4_HighLevel/manifest.json -- per-structure high-precision opt, freq, sp results
     |
     v
pipeline.result.json -- final summary
```

Stages are idempotent. Each writes its manifest and a signature. Re-running the same command reuses valid checkpoints and re-executes only incomplete stages.

### Inter-stage handoffs

| Source | Destination | Artifacts |
|--------|-------------|-----------|
| S0 record | S1 | forming_bonds (map-space), reaction_type |
| S1 manifest | S2 | `selected.xyz`, candidate energies |
| S2 manifest | S3 | `peb_peak` (TS guess), `intermediate_seed`, forming_bonds (xyz indices) |
| S3 manifest | S4 | Per-structure opt/sp/freq results, opt_xyz, status, usable_for_ml |
| S4 manifest | RPH_Postprocess | High-precision energies, optimized geometries, frequency data |

### Theory contract

```
S1: CREST/GFN2 conformational sampling
    + ORCA B97-3c single-point ranking
    + xTB GFN1 mRRHO thermochemistry
    No DFT optimization or DFT single-point for geometry refinement.

S3: ORCA B97-3c OPT (minimum) / OptTS (TS) + independent Freq (TS only)
    -> ORCA r2SCAN-3c SP
    Solvent: CPCM(acetone)

S4: Gaussian M062X/def2-SVP OPT / OptTS + independent Freq (TS only)
    -> ORCA wB97M-V/def2-TZVPP SP
    Solvent: CPCM(acetone)
```

All methods, paths, resources and timeouts are defined in `config/defaults.yaml` under `theory.s3_low_level` and `theory.s4_high_precision`.

### Output contract

The pipeline produces these authoritative files:

```
Output/<rx_id>/
├── S0_Mechanism/mechanism.json
├── S1_ConfSearch/product/manifest.json
├── S2_PEB/manifest.json
├── S3_LowLevel/manifest.json
├── S4_HighLevel/manifest.json
├── pipeline.state
├── pipeline.result.json
├── rph_v4.log
```

`forming_bonds` are authoritative in the S0 manifest and use 0-based XYZ indices. Each S3/S4 structure record retains input geometry, optimized geometry (when available), output files, energies, status and `usable_for_ml`.

### Example output tree

```
Output/RXN_000001/
├── .rph_schema.json
├── pipeline.result.json
├── pipeline.state
├── rph_v4.log
├── S0_Mechanism/
│   ├── mechanism.json
├── S1_ConfSearch/
│   ├── product/
│   │   ├── manifest.json
│   │   ├── selected.xyz
│   │   └── raw_censo/
│   └── manifest.json
├── S2_PEB/
│   ├── manifest.json
│   ├── scan_profile.json
│   ├── ts_guess.xyz
│   └── intermediate.xyz
├── S3_LowLevel/
│   ├── manifest.json
│   ├── status.json
│   ├── events.jsonl
│   ├── product/
│   │   ├── opt/
│   │   └── sp/
│   ├── product_int/
│   │   ├── opt/
│   │   └── sp/
│   └── product_ts/
│       ├── opt_ts/
│       ├── freq/
│       └── sp/
├── S4_HighLevel/
│   ├── manifest.json
│   ├── status.json
│   ├── events.jsonl
│   ├── s4.log
│   ├── product/
│   │   ├── opt/
│   │   └── sp/
│   ├── product_int/
│   │   ├── opt/
│   │   └── sp/
│   └── product_ts/
│       ├── opt_ts/
│       ├── freq/
│       └── sp/
```

### Feature extraction (external)

Feature extraction is not part of the RPH V4 pipeline. It is handled by the companion `RPH_Postprocess` package, which consumes the S4 manifests and produces ML-ready feature datasets.

Run it after RPH completes:

```bash
rph-features extract --rph-run ./Output/RXN_000001 --output ./Output/RXN_000001/S4_Data
```

RPH_Postprocess produces:
- `features_raw.csv` - all extracted features
- `features_mlr.csv` - ML-ready feature matrix
- `feature_meta.json` - provenance and config metadata

This separation keeps RPH focused on DFT computation while RPH_Postprocess handles all post-hoc analysis.

---

## Configuration guide

### QC executables

```yaml
executables:
  gaussian:
    path: "/opt/software/gaussian/g16/g16"
    root: "/opt/software/gaussian/g16"
  orca:
    path: "/opt/software/orca/orca"
    ld_library_path: "/opt/software/orca"
  crest:
    path: "/opt/software/crest/crest"
  xtb:
    path: "/opt/software/xtb/bin/xtb"
```

### Resources

```yaml
resources:
  mem: "32GB"
  nproc: 16
  orca_maxcore_safety: 0.65
```

### S3 low-level theory

```yaml
theory:
  s3_low_level:
    optimization:
      engine: orca
      method: B97-3c
      solvent: acetone
      solvent_model: CPCM
      route_minimum: "Opt"
      route_ts: "OptTS"
      frequency:
        enabled_for_ts: true
        imaginary_cutoff_cm1: -50.0
    single_point:
      engine: orca
      method: r2SCAN-3c
      solvent: acetone
      solvent_model: CPCM
```

### S4 high-precision theory

```yaml
theory:
  s4_high_precision:
    optimization:
      engine: gaussian
      method: M062X
      basis: def2-SVP
      solvent: acetone
      solvent_model: CPCM
      route_ts: "Opt=(TS,CalcFC,NoEigenTest)"
      grid: UltraFine
      scf: XQC
      frequency:
        enabled_for_ts: true
        imaginary_cutoff_cm1: -50.0
    single_point:
      engine: orca
      method: wB97M-V
      basis: def2-TZVPP
      aux_basis: def2/J
      solvent: acetone
      solvent_model: CPCM
```

### S1 CENSO-LITE

```yaml
step1:
  protocol: censo_lite
  censo_lite:
    crest:
      gfn_level: 2
      search_mode: imtd_smtd
      energy_window_kcal: 6.0
      solvent: acetone
    ranking:
      engine: orca
      method: B97-3c
      solvent: acetone
    xtb_thermo:
      enabled: true
      gfn_level: 1
      temperature_k: 298.15
    deduplication:
      backend: torsion_signature
      heavy_atom_rmsd_prefilter_A: 0.25
      torsion_bin_deg: 20.0
```

### PEB scan

```yaml
step2:
  scan:
    topology_guard_enabled: true
    scan_start_distance: 4.0
    scan_end_distance: 1.5
    scan_steps: 24
    retry_on_boundary_maximum: true
```

### Config file locations

- Main config: `config/defaults.yaml`
- Templates: `config/templates/`

---

## Tests

### Run Tests
```bash
# All tests
pytest -v tests/

# V4 protocol contract tests
pytest tests/test_v4_protocol_contract.py -v

# V4 checkpoint tests
pytest tests/test_v4_checkpoint.py -v

# V4 stage calculator tests
pytest tests/test_v4_stage_calculator.py -v

# With coverage
pytest --cov=rph_core --cov-report=html
```

### Import Style Check (CI Gate)
```bash
python scripts/ci/check_imports.py rph_core
```

### Notes
- `tests/conftest.py` adds the repo root to `sys.path`, so tests can run without editable install.
- Integration tests use mocked QC calculations and do not require real Gaussian/ORCA binaries.
- V3 tests have been archived to `tests/deprecated_v3/`.
- See [AGENTS.md](AGENTS.md) for verification commands.

---

## Troubleshooting

### 1. Gaussian executable not found

Error:
```
FileNotFoundError: Gaussian executable not found: /opt/software/gaussian/g16/g16
```

Fix:
- Verify `executables.gaussian.path` in `config/defaults.yaml`
- Ensure Gaussian is installed and environment variables are configured
- Try running `g16 < test.gjf` to validate availability

### 2. Out-of-memory errors

Error:
```
Error termination via Lnk1e in /root/g16/g16
```

Fix:
- Lower `resources.mem` and `resources.nproc`
- Reduce `orca_maxcore_safety` in resources
- Use fewer scan steps in `step2.scan`

### 3. CENSO-LITE conformer search fails

Error:
```
RuntimeError: S1 manifest has no selected candidate
```

Fix:
- Increase CREST sampling time or energy window in `step1.censo_lite.crest`
- Check that ORCA and xTB executables are configured correctly
- Verify the input CSV row has a valid product SMILES

### 4. S3/S4 optimization fails for a specific structure

The pipeline continues with remaining structures when `continue_on_structure_failure` is enabled (default). Check the per-structure status in the manifest:

- If `opt_status` is not `complete`, inspect the ORCA/Gaussian output in the opt directory
- If `sp_status` is not `complete`, check the single-point output
- If `ts_frequency_valid` is `false`, inspect the frequency output for the imaginary mode

### 5. Pipeline does not resume correctly

Delete the `pipeline.state` and target stage manifest, then re-run:

```bash
rm -f Output/RXN_000001/pipeline.state
bin/rph_run --csv data/trusted_reactions.csv --rx-id RXN_000001 --output ./Output/RXN_000001
```

### Logs

V4 logs to `rph_v4.log` in the output directory:

```bash
tail -f Output/RXN_000001/rph_v4.log
```

### Debug mode

```bash
bin/rph_run --csv data/trusted_reactions.csv --rx-id RXN_000001 --output ./Output/RXN_000001 --log-level DEBUG
```

### Getting help

1. Check GitHub Issues: https://github.com/PengYangchao0808/ReactionProfileHunter/issues
2. Provide logs and config files
3. Provide a minimal reproducible example

---

## Documentation Index

| Document | Contents |
|----------|----------|
| [`AGENTS.md`](AGENTS.md) | V4 coding guide -- build/test commands, code style, conventions |
| [`rph_core/AGENTS.md`](rph_core/AGENTS.md) | Core package architecture & code map |
| [`rph_core/steps/AGENTS.md`](rph_core/steps/AGENTS.md) | Step architecture & output contracts |
| [`rph_core/steps/conformer_search/AGENTS.md`](rph_core/steps/conformer_search/AGENTS.md) | CENSO-LITE conformer search |
| [`rph_core/steps/step2_retro/AGENTS.md`](rph_core/steps/step2_retro/AGENTS.md) | PEB scanner |
| [`rph_core/steps/step3_lowlevel/AGENTS.md`](rph_core/steps/step3_lowlevel/AGENTS.md) | S3 low-level QC |
| [`rph_core/steps/step4_highlevel/AGENTS.md`](rph_core/steps/step4_highlevel/AGENTS.md) | S4 high-precision QC |
| [`rph_core/utils/AGENTS.md`](rph_core/utils/AGENTS.md) | QC/IO/checkpoint infra reference |
| [`config/AGENTS.md`](config/AGENTS.md) | Config structure notes |
| [`tests/AGENTS.md`](tests/AGENTS.md) | Test organization and conventions |
| [`scripts/AGENTS.md`](scripts/AGENTS.md) | Scripts & CI tools notes |
| [`rph_benchmark/AGENTS.md`](rph_benchmark/AGENTS.md) | Benchmark suite documentation |
| [`docs/RPH_V4_COMPLETION_MASTER_PLAN.md`](docs/RPH_V4_COMPLETION_MASTER_PLAN.md) | V4 architecture master plan |
| [`docs/WSL_TEST_PLAN_V4.md`](docs/WSL_TEST_PLAN_V4.md) | WSL test plan and gates T0-T6 |

---

## Citation

If you use ReactionProfileHunter in research, please cite:

```bibtex
@software{reactionprofilehunter2025,
  title = {ReactionProfileHunter: Automated Reaction Mechanism Exploration},
  author = {Peng Yangchao},
  year = {2026},
  url = {https://github.com/PengYangchao0808/ReactionProfileHunter},
  version = {4.0.0}
}
```

### Resources

- Theory background: Houk group methodology, Grimme GFN-xTB papers, wB97M-V benchmark reports
- Tool docs:
  - https://gaussian.com/man/
  - https://sites.google.com/site/orcainputlibrary/
  - https://xtb-docs.readthedocs.io/

---

## Contributing

1. Fork the repo
2. Create a feature branch: `git checkout -b feature/amazing-feature`
3. Commit changes: `git commit -m "Add amazing feature"`
4. Push: `git push origin feature/amazing-feature`
5. Open a pull request

### Code Standards

- Use `pathlib.Path` for all paths (no string paths)
- Use `logging.getLogger(__name__)` or `LoggerMixin` (no `print()` in library code)
- Keep output directories idempotent: reuse existing outputs to avoid reruns
- All QC calls must go through `utils/qc_interface.py`
- Route every external QC job through `qc_jobs.py` and the existing QC interfaces
- **Absolute imports only** -- no multi-dot relative imports (e.g., `from rph_core.utils...` not `from ...utils`)
- Never add `subprocess.run` to a stage module
- Preserve XYZ, output logs and failure diagnostics when a job degrades
- A failed S4 job must never delete or replace usable S3 outputs
- Every stage writes a versioned manifest and is resumable through `pipeline.state`
- See [AGENTS.md](AGENTS.md) for complete coding guidelines

### Tests

- New features must include unit tests
- All tests must pass: `python -m pytest -q`
- Keep coverage above 80%

### Documentation

- Update relevant `AGENTS.md` files
- Add usage examples for new features
- Update version in `rph_core/version.py` and README badges

---

## License

This project is released under the MIT License. See `LICENSE`.

---

## Contact

- Author: Peng Yangchao
- Project: https://github.com/PengYangchao0808/ReactionProfileHunter
- Issues: https://github.com/PengYangchao0808/ReactionProfileHunter/issues

---

## Acknowledgements

- Prof. Houk: dual-level computation guidance
- Grimme group: xTB and CREST
- Gaussian and ORCA teams: QC tooling
- Community contributors

---

<div align="center">

**If this project is useful, please star it.**

</div>
