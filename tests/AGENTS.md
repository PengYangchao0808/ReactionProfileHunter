# tests/AGENTS.md

## OVERVIEW
**53 pytest files** covering import smoke tests, unit tests (mocked QC), integration tests, S4 degrade/contract tests, and path-compat/sandbox tests. External QC binaries (Gaussian/ORCA/xTB/CREST) are NOT required for most tests.

> **Refactoring Note (2026-03-11)**: 已从75个文件精简到53个，删除了16个冗余测试文件（主要是重复的formchk测试）。详见 `TEST_REFACTOR_REPORT.md`。

## WHERE TO LOOK
| Category | Files | Notes |
|----------|-------|-------|
| **Import smoke (fast CI)** | `test_imports_step4_features.py` | Import-only; no QC required; run first |
| **Forward scan wiring (NEW)** | `test_forward_scan_wiring.py`, `test_xtb_scan_input.py` | Profile routing, cleaner_data, scan params |
| **S4 extractor behavior** | `test_s4_extractor_degrade_behavior.py`, `test_s4_no_qc_execution.py`, `test_s4_meta_warnings_and_weights.py`, `test_s4_v62_final_verification.py` | No QC binaries needed |
| **S4 contracts** | `test_s4_mech_packager.py`, `test_s4_artifact_integration.py`, `test_s4_gedt_labeling.py`, `test_s4_path_compat.py`, `test_step4_cache_key.py` | Schema/contract checks |
| **S4 schema versioning** | `test_m2_schema_versioning.py`, `test_m2_step4_resume_semantics.py`, `test_m2_precursor_fallback.py` | |
| **S3 optimization** | `test_s3_checkpoint.py`, `test_m3_gaussian_templates.py`, `test_m3_qc_collection_mock.py`, `test_m3_qc_mock_simple.py` | Mock-based |
| **S2 retro scan** | `test_retro_scanner_v52.py`, `test_step2_path_compat.py` | |
| **S1 conformer** | `test_two_stage_conformer.py`, `test_e2e_precursor_leaving_group.py`, `test_small_molecule_cache.py`, `test_small_molecule_catalog.py` | |
| **QC interface** | `test_qc_interface_v52.py`, `test_orca_interface.py`, `test_qctaskrunner_integration.py` | Real-QC tests skip by default |
| **Sandbox/toxic path** | `test_sandbox_toxic_paths.py` | `pytest.skip()` at runtime if no xTB/CREST |
| **Degradation** | `test_degradation_final.py` | |
| **Integration** | `test_integration_final.py`, `test_mock_integration.py`, `test_mock_qc_e2e.py` | `test_mock_*` uses mocked QC |
| **fchk parsing** | `test_fchk_reader_multiline.py` | Formchk相关重复测试已合并移除 |
| **Misc utils** | `test_molecular_graph.py`, `test_molecule_utils.py`, `test_fragment_manipulation.py`, `test_sp_report.py`, `test_thermo_validation_v52.py`, `test_tsv_loader.py` | |
| **M4 mechanism** | `test_m4_mech_*.py`, `test_m4_qc_artifacts*.py`, `test_m4_template_structure.py` | |
| **Manual verify scripts** | `verify_sandbox.py` | Not pytest — manual use only |

## SKIP PATTERNS
```python
# Default skip for tests that need real QC binaries:
@pytest.mark.skipif(True, reason="Requires real ORCA environment")
def test_orca_real_run(): ...

# Runtime skip (used in test body when binary missing):
if not shutil.which("xtb"):
    pytest.skip("需要 XTB 可执行文件")

# Mocking QC calls (preferred for CI):
with unittest.mock.patch("rph_core.utils.qc_interface.subprocess.run") as mock_run:
    mock_run.return_value = ...
```

To enable real-QC integration tests: change `skipif(True, ...)` to `skipif(os.environ.get("ORCA_PATH") is None, ...)` and set the env var.

## FIXTURES
- **`tests/conftest.py`**: only adds repo root to `sys.path` — allows `import rph_core` without `pip install`.
- **Committed S1 fixtures** (`tests/tmp_v2_2_test/da_reaction/S1_Anchor/`): Diels-Alder example with full xTB/CREST + DFT artifacts (`.gjf`, `.log`, `.chk`, `.fchk`, `crest_conformers.xyz`, `cluster.xyz`, `isostat.log`). **Do NOT delete or gitignore.**
- **Per-test fixtures**: defined inside test modules (not in conftest); e.g., `base_config` in `test_qc_interface_v52.py`, `sample_xyz` in `test_orca_interface.py`.
- **No committed S2/S3 fixtures** — tests that need S2/S3 artifacts generate synthetic inputs or use mocks.

## CONVENTIONS
- `test_imports_*.py` + `test_s4_no_qc_execution.py` are the minimal fast-CI gate.
- Real-QC integration tests marked `skipif(True, ...)` by default; enable by changing condition.
- Do not write QC binary paths into tests — use env vars or `is_path_toxic()` checks.
- `verify_*.py` scripts are manual verification tools, not in the pytest suite.

## COMMANDS
```bash
# Fast offline tests only (CI gate)
pytest tests/test_imports_step4_features.py tests/test_s4_no_qc_execution.py -v

# Core S4 tests (recommended for most validation)
pytest tests/test_s4_*.py tests/test_degradation_final.py -v

# All tests
pytest -v tests/

# Specific step
pytest -v tests/test_s3_checkpoint.py tests/test_retro_scanner_v52.py

# S4 contract tests
pytest -v tests/test_s4_*.py tests/test_m2_*.py tests/test_m4_*.py
```

## DEPRECATED TESTS
已移除的冗余测试文件存放在 `tests/deprecated/` 目录：
- 9个重复的formchk/try_formchk测试
- test_degradation.py（保留test_degradation_final.py）
- test_integration.py（保留test_integration_final.py）
- verify_sandbox.py（手动验证脚本）

这些文件可安全删除，如需恢复可从deprecated目录移回。

## ANTI-PATTERNS
- Hardcoding absolute paths in test bodies (especially paths with spaces/brackets — use toxic-path test helpers).
- Treating skip as test failure — many real-QC tests will always skip in CI.
- Modifying `tests/tmp_v2_2_test/` committed fixtures — treat as read-only reference data.
- **Duplicating test logic** — if a feature is already tested, extend existing test rather than creating new file.
- **Creating multiple test files for the same module** — consolidate related tests into single comprehensive test file.
