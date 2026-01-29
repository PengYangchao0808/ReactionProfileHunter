# tests/AGENTS.md

## OVERVIEW
pytest 测试覆盖模块导入、数据流编排、降级策略、以及 toxic path/sandbox 行为。

## WHERE TO LOOK
- Integration: `test_integration.py`, `test_integration_final.py`
- Degradation: `test_degradation.py`, `test_degradation_final.py`
- Toxic path/sandbox: `test_sandbox_toxic_paths.py`
- ORCA: `test_orca_interface.py`
- QC task runner: `test_qctaskrunner_integration.py`
- Step4: `test_imports_step4_features.py`, `test_step4_cache_key.py`, `test_s4_artifact_integration.py`

## RUN
```bash
pytest -v tests/
pytest -v tests/test_sandbox_toxic_paths.py
pytest -v tests/test_integration.py
```

## CONVENTIONS
- 依赖真实外部程序（XTB/CREST/Gaussian/ORCA）的用例常会 skip；不要把 skip 当失败。
- 仓库还保留 verify_*.py 作为手工验证脚本（不等同 pytest）。

## ANTI-PATTERNS
- 测试里写死本机路径（尤其是包含空格/括号的路径；应覆盖 toxic path 逻辑）
