# ReactionProfileHunter Changelog

## v2.0.0 Release (Current)

### Release Date: March 18, 2025

---

## New Features

### 1. Forward-Scan Transition State Search
- **Description**: Implement xTB native `$scan` for arbitrary cycloaddition transition state search
- **Supported Reaction Types**:
  - `[4+3]` - 4+3 cycloaddition (forward-scan)
  - `[4+2]` - Diels-Alder reaction (forward-scan)
  - `[3+2]` - 1,3-dipolar cycloaddition (forward-scan)
  - `[5+2]` - 5+2 cycloaddition (retro-scan)
- **Related Files**:
  - `rph_core/steps/step2_retro/retro_scanner.py`
  - `rph_core/utils/xtb_runner.py`

### 2. Reaction Profiles Configuration
- **Description**: Configuration-driven scan parameters per reaction type
- **Config Example**:
  ```yaml
  reaction_profiles:
    "[4+3]_default":
      forming_bond_count: 2
      s2_strategy: forward_scan
      scan:
        scan_start_distance: 1.8
        scan_end_distance: 3.2
        scan_steps: 20
  ```
- **Related Files**: `config/defaults.yaml`

### 3. Step1 Activation Features
- Boltzmann/Gibbs weighted energy
- Conformer entropy and flexibility
- Leaving group geometry
- Alpha-H gating mechanism
- **Related Files**:
  - `rph_core/steps/step4_features/extractors/step1_activation.py`

### 4. Step2 Cyclization Features
- Kinetics/thermochemistry
- Transition state geometry (forming bonds)
- CDFT metrics (Fukui/Dual Descriptor/QTAIM)
- **Related Files**:
  - `rph_core/steps/step4_features/extractors/step2_cyclization.py`

### 5. Multiwfn Integration
- Tier-1 support (Fukui/Dual Descriptor/QTAIM)
- Non-interactive batch processing
- Failure tolerance
- Caching mechanism
- **Related Files**:
  - `rph_core/utils/multiwfn_runner.py`
  - `rph_core/steps/step4_features/extractors/multiwfn_features.py`

### 6. Small Molecule Cache System
- Global small molecule reuse
- Cache key generation
- Detection and matching
- **Related Files**:
  - `rph_core/utils/small_molecule_cache.py`
  - `rph_core/utils/small_molecule_catalog.py`

---

## Refactoring & Optimizations

### 1. Directory Structure Refactor
- **Change**: Renamed `S1_Product` to `S1_ConfGeneration`
- **Reason**: Clearer semantic expression
- **Related Commits**:
  - `b37a2d1` - refactor(s1): rename S1_Product to S1_ConfGeneration

### 2. S2 Path Updates
- Update S2-related path references
- Ensure backward compatibility
- **Related Commits**:
  - `b37a2d1` - refactor(s1): rename S1_Product to S1_ConfGeneration and update S2 paths

### 3. Project Structure Cleanup
- Remove deprecated backup directory (`rph_core_backup_20260115/`)
- Update `.gitignore`
- Optimize test file organization
- **Related Commits**:
  - `ca5e2f9` - chore: cleanup project structure and update .gitignore

### 4. Anchor Integration with SmallMoleculeCache
- Integrate SmallMoleculeCache for small molecule reuse
- **Related Commits**:
  - `19ebe2a` - feat(anchor): integrate SmallMoleculeCache for small molecule reuse (Task 4)

---

## Testing Enhancements

### New Test Files
| File | Description |
|------|-------------|
| `tests/test_s1_progress_parser.py` | S1 progress parser tests |
| `tests/test_forward_scan_wiring.py` | Forward-scan wiring tests |
| `tests/test_geometry_guard.py` | Geometry guard tests |
| `tests/test_s2_boundary_degrade.py` | S2 boundary degradation tests |
| `tests/test_s3_checkpoint.py` | S3 checkpoint tests |
| `tests/test_s4_extractor_degrade_behavior.py` | S4 extractor degradation behavior tests |
| `tests/test_s4_meta_warnings_and_weights.py` | S4 metadata warnings and weights tests |
| `tests/test_s4_v62_final_verification.py` | S4 V6.2 final verification tests |
| `tests/test_qc_interface_gaussian_failures.py` | Gaussian failure handling tests |
| `tests/test_xtb_path_integration.py` | xTB path integration tests |
| `tests/test_xtb_scan_input.py` | xTB scan input tests |
| `tests/test_xtb_scan_params.py` | xTB scan parameters tests |
| `tests/test_cleaner_adapter.py` | Cleaner adapter tests |
| `tests/test_dataset_loader.py` | Dataset loader tests |
| `tests/test_fchk_reader_multiline.py` | FCHK multiline reader tests |
| `tests/test_orca_interface.py` | ORCA interface tests |
| `tests/test_gau_xtb_interface.py` | Gaussian-xTB interface tests |

### Test Optimizations
- Archive redundant tests to `tests/deprecated/`
- Add fast CI gate tests
- Improve Mock QC integration

---

## Configuration Updates

### New Config Options

```yaml
# Reaction Profiles
reaction_profiles:
  "[4+3]_default":
    forming_bond_count: 2
    s2_strategy: forward_scan
    scan:
      scan_start_distance: 1.8
      scan_end_distance: 3.2
      scan_steps: 20
      scan_mode: concerted
      scan_force_constant: 0.5

# Step1 Feature Extraction
step1:
  features:
    activation:
      boltzmann_weighted_energy: true
      conformer_entropy: true
      leaving_group_geometry: true
      alpha_h_gating: true

# Multiwfn Configuration
step4:
  multiwfn:
    enabled: true
    tiers:
      - fukui
      - dual_descriptor
      - qtaim
```

---

## Documentation Updates

### New Documents
| Document | Description |
|----------|-------------|
| `README.zh-CN.md` | Chinese version of README |
| `docs/S2_S1_S2_2_Plan_20260315.md` | S2 planning document |
| `docs/S2_XTB_Path_Integration_Report.md` | xTB path integration report |
| `docs/S4_FEATURES_SUMMARY.md` | S4 feature summary |
| `docs/S4_REPORT_COMPLETE.md` | S4 completion report |
| `docs/TESTING_GUIDE.md` | Testing guide |
| `docs/QC_VALIDATION_GUIDE.md` | QC validation guide |
| `docs/DEV_STATUS_REPORT_20250311.md` | Development status report |
| `docs/S2_S3_Architecture_Report.md` | S2/S3 architecture report |
| `docs/MODIFICATION_REPORT_20250311.md` | Modification report |
| `docs/DUAL_LEVEL_STRATEGY_SUMMARY.md` | Dual-level strategy summary |

---

## Dependency Updates

### Python Dependencies
- Maintain Python 3.8+ compatibility
- Optimize package management configuration

### QC Tool Support
- Gaussian 16
- ORCA
- xTB
- CREST
- Multiwfn

---

## Known Issues

1. **NBO Analysis**: Disabled by default to reduce runtime; can be enabled via configuration
2. **Large Molecule Conformer Search**: May require longer computation time
3. **Parallelization**: Current version primarily supports sequential execution

---

## Upgrade Guide

### Upgrading from v6.1.x

1. Update config file to support new `reaction_profiles`:
   ```yaml
   run:
     reaction_type: "[4+3]_default"  # New
   ```

2. To enable Multiwfn analysis:
   ```yaml
   step4:
     multiwfn:
       enabled: true
   ```

3. Run tests to verify:
   ```bash
   pytest -v tests/
   ```

---

## Contributors

- QCcalc Team

---

## Historical Versions

### v6.1.0
- Simplified S4 module: feature extraction only (no QC execution)
- Fixed output format: features_raw.csv, features_mlr.csv, feature_meta.json
- NBO disabled by default to reduce runtime
- Improved error handling and logging

### v6.0.0
- Initial release
- Four-step pipeline architecture

---

**Last Updated**: March 18, 2025

**Version**: 6.2.0
