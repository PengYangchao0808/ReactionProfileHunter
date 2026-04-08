# ReactionProfileHunter

[English](README.md) | [中文](README.zh-CN.md)

<div align="center">

**Product-driven reaction mechanism exploration and feature extraction**

[![Version](https://img.shields.io/badge/version-2.0.0-blue.svg)](https://github.com/yourusername/ReactionProfileHunter)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Code Style](https://img.shields.io/badge/code%20style-agents%20ready-success.svg)](AGENTS.md)

</div>

> **Project Stats:** ~50k lines, 198 Python files | **Agent Guide:** See [AGENTS.md](AGENTS.md) for coding standards

---

## Overview

ReactionProfileHunter (RPH) is a product-driven automated reaction mechanism pipeline designed for transition state search, geometry optimization, and feature extraction.

### Core capabilities

- Product-driven strategy: start from products and search reaction pathways backward
- Four-step pipeline: anchor/conformer search -> retro scan -> TS optimization -> feature extraction
- Dual-level computation: xTB pre-optimization -> B3LYP geometry optimization -> wB97X-D3BJ high-accuracy single point
- Multi-engine support: Gaussian, ORCA, xTB, CREST, Multiwfn
- Feature engineering: reaction barriers, geometric parameters, electronic descriptors
- Checkpoint/resume support

### Use cases

- Mechanistic studies of organic reactions (cycloadditions, rearrangements, substitutions)
- Transition state prediction and validation
- Reaction dataset building and high-throughput screening
- Machine-learning feature extraction for reaction properties

---

## Key features

### v2.0.0 highlights (forward-scan + major refactoring)

- Step1 activation features: Boltzmann/Gibbs weighted energy, conformer entropy/flexibility, leaving-group geometry, alpha-H gating
- Step2 cyclization features: kinetics/thermochemistry, TS geometry (forming bonds), CDFT metrics (Fukui/Dual Descriptor/QTAIM)
- **Forward-scan TS search (NEW)**: xTB native `$scan` for arbitrary cycloadditions ([4+2], [4+3], [3+2], [5+2])
- **Reaction profiles (NEW)**: configuration-driven scan parameters per reaction type
- GEDT/CDFT enhancements: forming-bond-based charge transfer, eV unit locking, range validation
- Multiwfn integration: tier-1 support (Fukui/Dual Descriptor/QTAIM), non-interactive batch, failure-tolerant, caching
- Strict contract: feature_meta records config snapshot and provenance; missing features degrade with NaN + warning

### v1.x historical (legacy versions from v6.x era)

- Simplified S4 module: feature extraction only (no QC execution)
- Fixed output format: features_raw.csv, features_mlr.csv, feature_meta.json
- NBO disabled by default to reduce runtime; can be enabled when needed
- Improved error handling and logging

### Technical highlights

- Molecular autonomy: each molecule has its own working directory (S1_ConfGeneration/[Molecule]/)
- OPT-SP coupling loop: optimized geometry and single-point energy tightly coupled
- Path compatibility: auto adapts to legacy directory layouts
- Sandbox isolation: each QC task runs in an isolated sandbox

---

## Quick start

### System requirements

| Component | Requirement |
|------|------|
| Operating system | Linux (Ubuntu 20.04+ or CentOS 7+ recommended) |
| Python | 3.8 or newer |
| Memory | 16 GB minimum (64 GB recommended) |
| CPU | Multi-core CPU (16+ cores recommended) |

### Quantum chemistry dependencies

Install at least one stack:

Option 1: Gaussian
- Gaussian 16 - geometry optimization
- xTB - pre-optimization (optional but recommended)
- CREST - conformer search (optional)

Option 2: ORCA
- ORCA - geometry optimization and single point
- xTB - pre-optimization (optional but recommended)
- CREST - conformer search (optional)

Option 3: Mixed (recommended)
- Gaussian 16 - geometry optimization
- ORCA - high-accuracy single point
- xTB + CREST - pre-optimization and conformer search

### Installation

1. Clone the repo
   ```bash
   git clone https://github.com/yourusername/ReactionProfileHunter.git
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

#### Method 1: CLI (recommended)

```bash
# Run a single reaction from a SMILES string
bin/rph_run --smiles "C=C(C)C(=O)O" --output ./Output/rx_manual

# Specify reaction type for forward-scan (NEW)
bin/rph_run --smiles "C=C(C)C(=O)O" --reaction-type "[4+3]_default" --output ./Output/rx_4p3
```

Supported reaction types:
- `[5+2]_default` - 5+2 cycloaddition (retro_scan, default)
- `[4+3]_default` - 4+3 cycloaddition (forward_scan)
- `[4+2]_default` - Diels-Alder (forward_scan)
- `[3+2]_default` - 1,3-dipolar cycloaddition (forward_scan)

`rx_id` comes from the config (`run.single.rx_id`) or from the dataset id column.

#### Method 2: Config-driven

1) Set run mode in `config/defaults.yaml`:

```yaml
run:
  source: single
  single:
    product_smiles: "C=C(C)C(=O)O"
    rx_id: "rx_001"
```

2) Run:

```bash
bin/rph_run
```

#### Method 3: Python API

```python
from pathlib import Path
from rph_core.orchestrator import ReactionProfileHunter

hunter = ReactionProfileHunter(config_path="config/defaults.yaml")

result = hunter.run_pipeline(
    product_smiles="C=C(C)C(=O)O",
    work_dir=Path("./Output/rx_001")
)

if result.success:
    print(f"OK: {result.features_csv}")
else:
    print(f"Failed: {result.error_message}")
```

---

## Detailed documentation

### Architecture overview

ReactionProfileHunter uses a sequential four-step pipeline:

```
S1: Anchor/conformer search
    -> (product_min.xyz, precursor_min.xyz)
S2: Retro scan
    -> (ts_guess.xyz, reactant_complex.xyz)
S3: TS optimization/rescue
    -> (ts_final.xyz, reactant_opt/, NBO artifacts)
S4: Feature extraction and packaging
    -> (features_raw.csv, features_mlr.csv, feature_meta.json)
```

#### Step breakdown

| Step | Purpose | Core modules | Key outputs |
|------|--------|--------------|-------------|
| S1 | Build 3D structure from SMILES; conformer search; DFT optimization | `steps/anchor/`, `steps/conformer_search/` | `product_min.xyz`, `precursor_min.xyz` |
| S2 | Retro/Forward scan from products; TS guess; reactant complex | `steps/step2_retro/` (retro_scan or forward_scan) | `ts_guess.xyz`, `reactant_complex.xyz` |
| S3 | TS optimization and frequency analysis; reactant optimization; rescue strategies | `steps/step3_opt/` | `ts_final.xyz`, `reactant_opt/` |
| S4 | Energy extraction; geometric features; NBO/FMO analysis | `steps/step4_features/` | `features_raw.csv`, `feature_meta.json` |

### Inputs and outputs

#### Required inputs

- Product SMILES string (e.g., "C=C(C)C(=O)O")
- Reaction ID (optional; used for output directory naming)

#### Strict output contract

**S1 outputs (required)**
- `S1_ConfGeneration/product/product_min.xyz`
- `S1_ConfGeneration/precursor/precursor_min.xyz` (optional)

**S2 outputs (required; both)**
- `S2_Retro/ts_guess.xyz` - TS guess structure
- `S2_Retro/reactant_complex.xyz` - reactant complex

**S3 outputs (required)**
- `S3_TS/ts_final.xyz` - optimized TS structure
- `S3_TS/reactant_sp.xyz` - optimized reactant geometry
- `S3_TS/reactant_opt/standard/` or `S3_TS/reactant_opt/rescue/` - OPT+Freq run directory

**S4 outputs (required)**
- `S4_Data/features_raw.csv` - raw features
- `S4_Data/features_mlr.csv` - ML-ready features
- `S4_Data/feature_meta.json` - feature metadata (version, provenance)

**S4 optional outputs (NBO)**
- `S4_Data/qc_nbo.37` - NBO file (if found in S3)

#### Example output tree

```
Output/rx_001/
├── S1_ConfGeneration/
│   ├── product/
│   │   └── product_min.xyz
│   └── precursor/
│       └── precursor_min.xyz
├── S2_Retro/
│   ├── ts_guess.xyz
│   └── reactant_complex.xyz
├── S3_TS/
│   ├── ts_final.xyz
│   ├── reactant_sp.xyz
│   └── reactant_opt/
│       └── standard/
│           ├── input.gjf
│           ├── output.log
│           └── *.fchk
├── S4_Data/
│   ├── features_raw.csv
│   ├── features_mlr.csv
│   ├── feature_meta.json
│   └── qc_nbo.37
└── rph.log
```

### Configuration guide

#### 1) QC executables (`executables`)
```yaml
executables:
  gaussian:
    path: "/root/g16/g16"
    use_wrapper: true
    wrapper_path: "./scripts/run_g16_worker.sh"
```

#### 2) Resources (`resources`)
```yaml
resources:
  mem: "64GB"
  nproc: 16
  orca_maxcore_safety: 0.8
```

#### 3) Theory levels (`theory`)
```yaml
  preoptimization:
    enabled: true
    gfn_level: 2
    overlap_threshold: 1.0
  optimization:
    method: B3LYP
    basis: def2-SVP
    dispersion: GD3BJ
    engine: gaussian
  single_point:
    method: wB97X-D3BJ
    basis: def2-TZVPP
    engine: orca
```

#### 4) Reaction profiles (`reaction_profiles`) - NEW
```yaml
reaction_profiles:
  "[4+3]_default":
    forming_bond_count: 2
    s2_strategy: forward_scan       # Use xTB $scan
    scan:
      scan_start_distance: 1.8      # Initial bond distance (Å)
      scan_end_distance: 3.2       # Final bond distance (Å)
      scan_steps: 20                # Number of scan steps
      scan_mode: concerted          # or "sequential"
      scan_force_constant: 0.5      # Constraint force constant

  "[4+2]_default":
    forming_bond_count: 2
    s2_strategy: forward_scan
    scan:
      scan_start_distance: 2.0
      scan_end_distance: 3.5
      scan_steps: 15

  "[5+2]_default":
    forming_bond_count: 2
    s2_strategy: retro_scan         # Legacy backward scan

  "_universal":
    forming_bond_count: 2
    s2_strategy: forward_scan
    scan:
      scan_start_distance: 2.2
      scan_end_distance: 3.5
```

#### 5) Optimization control (`optimization_control`)
```yaml
optimization_control:
  max_cycles: 100
  convergence:
    level: normal
  hessian:
    initial: calcfc
    recalc_every: 10
```

#### Config file locations

- Main config: `config/defaults.yaml`
- Templates: `config/templates/`

---

## Advanced usage

### Batch mode

Create a TSV file `reactions.tsv`:

```tsv
rx_id	product_smiles	precursor_smiles
rx_001	C=C(C)C(=O)O	C=C(C)C(=O)OC
rx_002	CC(=O)OC	CC(=O)Cl
```

Set `config/defaults.yaml`:

```yaml
run:
  source: batch
  batch:
    input_file: "reactions.tsv"
    smiles_column: "product_smiles"
```

Run:

```bash
bin/rph_run
```

### Dataset mode

```yaml
run:
  source: dataset
  dataset:
    format: "csv"
    path: "data/reaxys_cleaned.csv"
    product_smiles_col: "product_smiles_main"
    id_col: "rx_id"
    delimiter: ","
```

### Resume and rehydrate

RPH resumes from checkpoints by default. Re-run the same command to continue:

```bash
bin/rph_run --smiles "C=C(C)C(=O)O" --output ./Output/rx_001
```

To control resume behavior, edit config:

```yaml
run:
  resume: true
  resume_rehydrate: true
  resume_rehydrate_policy: best_effort
```

### Custom parameters

Create a custom config file and pass it via `--config`:

```bash
bin/rph_run --config config/custom.yaml
```

Example overrides:

```yaml
theory:
  optimization:
    method: M06-2X
resources:
  nproc: 32
```

### Enable NBO analysis

Config option:

```yaml
step3:
  reactant_opt:
    enable_nbo: true
```

NBO file collection rules:
- Search in `S3_TS` subdirectories:
  - `nbo_analysis/`
  - `nbo/`
  - `reactant_opt/standard/`
  - `reactant_opt/rescue/`
- Recognized extensions: `*.37`, `*.nbo`, `*.nbo7`
- Copy to `S4_Data/qc_nbo.37` (normalized name)

---

## Tests

### Run Tests
```bash
# All tests
pytest -v tests/

# Single test file
pytest tests/test_s4_no_qc_execution.py -v

# Single test function
pytest tests/test_s4_no_qc_execution.py::test_extractor_degrades_gracefully -v

# Fast CI gate (import smoke + no-QC tests)
pytest tests/test_imports_step4_features.py tests/test_s4_no_qc_execution.py -v

# S4 contract tests only
pytest tests/test_s4_*.py tests/test_m2_*.py tests/test_m4_*.py -v

# With coverage
pytest --cov=rph_core --cov-report=html
```

### Import Style Check (CI Gate)
```bash
# Block multi-dot relative imports (exit 1 on violation)
python scripts/ci/check_imports.py rph_core
```

### Notes
- `tests/conftest.py` adds the repo root to `sys.path`, so tests can run without editable install.
- Integration tests use mocked QC calculations and do not require real Gaussian/ORCA binaries.
- See [tests/AGENTS.md](tests/AGENTS.md) for detailed test organization.

---

## Troubleshooting

### 1. Gaussian executable not found

Error:
```
FileNotFoundError: Gaussian executable not found: /root/g16/g16
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
- Use a smaller basis set (e.g., def2-SVP -> 6-31G*)
- Enable ORCA `maxcore` limits

### 3. Conformer search fails

Error:
```
ConformerSearchError: No conformers found below energy threshold
```

Fix:
- Increase CREST sampling time
- Loosen the energy window
- Verify SMILES input

### 4. TS optimization fails

Error:
```
TSOptimizationError: TS optimization failed after 5 attempts
```

Fix:
- Inspect `ts_guess.xyz` (imaginary frequencies, geometry)
- Enable rescue strategies (default)
- Adjust TS guess geometry

### 5. Missing NBO data

Error:
```
Warning: NBO file not found in S3 subdirectories
```

Fix:
- Ensure `step3.reactant_opt.enable_nbo: true`
- Check for `*.37` or `*.nbo` under `S3_TS/reactant_opt/standard/`
- Confirm NBO succeeded in Gaussian/ORCA logs

### Logs

Logs are stored in `Output/[rx_id]/rph.log`.

```bash
tail -f Output/rx_001/rph.log
```

### Debug mode

```bash
bin/rph_run --log-level DEBUG
```

### Getting help

1. Check GitHub Issues: https://github.com/yourusername/ReactionProfileHunter/issues
2. Provide logs and config files
3. Provide a minimal reproducible example

---

## Further reading

### Documentation Index

| Document | Contents |
|------|------|
| [`AGENTS.md`](AGENTS.md) | **Agentic coding guide** — build/test commands, code style, conventions |
| [`rph_core/steps/AGENTS.md`](rph_core/steps/AGENTS.md) | Step architecture notes |
| [`config/AGENTS.md`](config/AGENTS.md) | Config structure notes |
| [`tests/AGENTS.md`](tests/AGENTS.md) | Test organization and conventions |
| [`scripts/AGENTS.md`](scripts/AGENTS.md) | Scripts & CI tools notes |
| [`docs/DUAL_LEVEL_STRATEGY_SUMMARY.md`](docs/DUAL_LEVEL_STRATEGY_SUMMARY.md) | Dual-level computation strategy |
| [`docs/QUICK_START_DUAL_LEVEL.md`](docs/QUICK_START_DUAL_LEVEL.md) | Quick start for dual-level strategy |
| [`docs/BUGFIX_STEP3_QCRESULT.md`](docs/BUGFIX_STEP3_QCRESULT.md) | S3 QCResult refactor notes |
| [`docs/S4_FEATURES_SUMMARY.md`](docs/S4_FEATURES_SUMMARY.md) | S4 feature extraction summary |
| [`docs/TESTING_GUIDE.md`](docs/TESTING_GUIDE.md) | Testing guide and environment setup |

### Citation

If you use ReactionProfileHunter in research, please cite:

```bibtex
@software{reactionprofilehunter2025,
  title = {ReactionProfileHunter: Automated Reaction Mechanism Exploration},
  author = {Your Name},
  year = {2025},
  url = {https://github.com/yourusername/ReactionProfileHunter},
  version = {2.0.0}
}
```

### Resources

- Theory background: Houk group methodology, Grimme GFN-xTB papers, wB97X-D3 benchmark reports
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
- **Absolute imports only** — no multi-dot relative imports (e.g., `from rph_core.utils...` not `from ...utils`)
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

- Author: QCcalc Team
- Email: your.email@example.com
- Project: https://github.com/yourusername/ReactionProfileHunter
- Issues: https://github.com/yourusername/ReactionProfileHunter/issues

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
