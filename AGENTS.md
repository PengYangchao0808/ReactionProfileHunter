# AGENTS.md (ReactionProfileHunter)

## OVERVIEW
ReactionProfileHunter 是 product-first 反应机理流水线：S1(锚定/构象) → S2(逆向扫描) → S3(TS优化/救援) → S4(特征提取)。

## ENTRYPOINTS
- CLI: `bin/rph_run`
- Python CLI main: `rph_core/orchestrator.py`（含 `main()` / `__main__`）
- 集成跑法: `run_auto_test.py`（按 step / resume）

## STRUCTURE
- `rph_core/`                核心库（orchestrator + steps + utils）
- `config/`                  defaults + Gaussian 模板
- `tests/`                   pytest（含 toxic path / degradation / integration）
- `bin/`                     CLI wrapper（本地运行时把 repo root 加入 sys.path）
- `scripts/`                 外部程序 wrapper（例如 Gaussian g16 worker）
- `ci/`                      CI 辅助脚本与导入风格检查
- `test_results/`            运行产物/样例输出（非源码，但要纳入理解）
- `rph_core_backup_20260115/`历史冻结备份（对照用途）

## WHERE TO LOOK
| Task | Location |
|------|----------|
| Pipeline wiring / data flow | `rph_core/orchestrator.py` |
| Step2 Retro scan | `rph_core/steps/step2_retro/retro_scanner.py` |
| Step3 TS optimize/rescue | `rph_core/steps/step3_opt/` |
| Step4 Feature extraction | `rph_core/steps/step4_features/` |
| QC execution + sandbox/toxic path | `rph_core/utils/qc_interface.py` |
| ORCA / xTB runners | `rph_core/utils/orca_interface.py`, `rph_core/utils/xtb_runner.py` |
| Checkpoint/resume | `rph_core/utils/checkpoint_manager.py` |

## COMMANDS
```bash
pip install -e ".[dev]"
pytest -v tests/
python run_auto_test.py --step ALL
python run_auto_test.py --step ALL --resume
```

## TEST DISCOVERY NOTE
- `pytest.ini` 限制 `testpaths=tests`：repo root 下 `test_*.py` 多为手工验证脚本，不应被 pytest 自动收集。

## PROJECT GOTCHAS (HIGH SIGNAL)
- **Toxic path**：空格与 `[](){}`
  - 相关：`rph_core/utils/qc_interface.py` + `tests/test_sandbox_toxic_paths.py`
  - 任何落盘 QC 输出路径都应复用该逻辑；不要各处自写 subprocess cwd。
- **Step2 输出契约**：必须同时产出 `ts_guess.xyz` 与 `reactant_complex.xyz`
  - 后者被 Step3(QST2 rescue) 与 Step4(畸变/片段)依赖。
- `test_results/` 是"运行产物/样例"，不是主实现；修复应优先改 `rph_core/`。
- `rph_core_backup_20260115/` 为冻结备份：用于对照，不作为开发目标。

## ANTI-PATTERNS
- 在 steps 内重复实现 QC runner（应复用 `rph_core/utils`）
- 修 bug 顺手大重构（repo 内有备份与大体积产物目录，diff 风险高）
