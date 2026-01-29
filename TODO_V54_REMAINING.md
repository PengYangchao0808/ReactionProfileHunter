# V5.4: Complete removal of NMR/Hirshfeld

## pytest execution (solved)
- Created `tests/conftest.py` to add project root to sys.path for imports
- This resolves `ModuleNotFoundError: No module named 'rph_core'` when running tests

## Next steps

### 1. Clean up core code residues
- [ ] `rph_core/utils/qc_interface.py`: remove all NMR/Hirshfeld parameters and behavior
- [ ] `rph_core/utils/optional_qc_tasks.py`: delete file entirely

### 2. Config/template cleanup
- [ ] `config/defaults.yaml`: remove entire `qc_tasks:` section

### 3. S4 test unification to NBO-only
- [ ] `tests/test_m4_qc_artifacts_mech_index.py`: remove hirshfeld/nmr assertions
- [ ] `tests/test_m4_qc_artifacts_structure.py`: remove hirshfeld tests
- [ ] `tests/test_s4_mech_packager.py`: remove NMR/Hirshfeld artifact creation
- [ ] `tests/test_m4_template_structure.py`: remove NMR/Hirshfeld template assertions

### 4. M3/QC tests unification to NBO-only
- [ ] `tests/test_qc_interface_v52.py`: remove all NMR/Hirshfeld tests
- [ ] `tests/test_m3_gaussian_templates.py`: remove NMR template tests

### 5. Mock e2e fixes
- [ ] `tests/test_mock_qc_e2e.py`: ensure XYZ format is valid, step_dirs meet contract

### 6. Full repo residue scan
- [ ] Run `rg "NMR_GIAO|nmr_giao|NMR_WHITELIST|NMR=GIAO" -S`
- [ ] Run `rg "HIRSHFELD|hirshfeld|Pop=Hirshfeld|parse_hirsh" -S`
- [ ] Run `rg "gaussian_nmr.gjf|gaussian_hirshfeld.gjf" -S`

### 7. pytest staged run (small to large)
- [ ] `python -m pytest tests/test_m3_qc_mock_simple.py`
- [ ] `python -m pytest tests/test_qc_interface_v52.py`
- [ ] `python -m pytest tests/test_m4_qc_artifacts_mech_index.py`
- [ ] `python -m pytest tests/test_s4_mech_packager.py`
- [ ] `python -m pytest` (full suite)

### 8. Final verification
All above steps must pass before V5.4 is considered complete.
