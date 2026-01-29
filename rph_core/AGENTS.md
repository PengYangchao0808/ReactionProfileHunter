# rph_core/AGENTS.md

## OVERVIEW
核心实现区：orchestrator 负责串联 S1~S4；steps 实现业务；utils 提供 QC/IO/日志/恢复。

## WHERE TO LOOK
- Orchestrator: `orchestrator.py`
- Steps: `steps/`
- QC + infra: `utils/`
- Package: `__init__.py`

## CONVENTIONS
- 路径统一用 `pathlib.Path` 传递。
- 日志用 `LoggerMixin`（禁止在核心路径 print）。
- 输出目录尽量幂等：存在可复用输出则跳过重算（避免长任务重复）。

## ANTI-PATTERNS
- step 内直接 `subprocess.run` 调外部程序而绕过 `utils/qc_interface.py`
- 通过复制粘贴把 `test_results/` 行为固化为代码逻辑
