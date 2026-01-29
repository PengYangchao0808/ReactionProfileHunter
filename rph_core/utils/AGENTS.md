# rph_core/utils/AGENTS.md

## OVERVIEW
跨步骤基础设施：QC 执行与失败分类、sandbox/toxic path、安全落盘、日志、checkpoint、外部程序接口。

## WHERE TO LOOK
- sandbox/toxic path + QC facade: `qc_interface.py`
- retry/failure typing: `qc_runner.py`
- task runner / orchestration helper: `qc_task_runner.py`
- ORCA: `orca_interface.py`
- xTB: `xtb_runner.py`
- logging: `log_manager.py`
- checkpoint: `checkpoint_manager.py`

## PROJECT-SPECIFIC RULES
- Toxic path 规则：空格与 `[](){}`
  - 任何会创建/使用 QC 输出目录的逻辑必须复用 `is_path_toxic()`/sandbox 执行。
- 外部可执行文件查找：优先 config → PATH → fallback（xTB 已实现该模式）。

## ANTI-PATTERNS
- 在多个模块各自实现"可执行文件查找/重试/日志格式"
- 在异常路径吞掉 stderr/stdout 导致不可诊断
