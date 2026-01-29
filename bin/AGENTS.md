# bin/AGENTS.md

## OVERVIEW
本地 CLI wrapper：保证“未安装 editable 包”的情况下也能从 repo root 运行 `rph_core.orchestrator:main`。

## WHERE TO LOOK
- `rph_run` / `rph`：插入 `PROJECT_ROOT` 到 `sys.path` 后调用 `rph_core.orchestrator.main()`

## CONVENTIONS
- `bin/` 只做薄封装；业务逻辑留在 `rph_core/`。

## ANTI-PATTERNS
- 在 wrapper 内解析/修改 pipeline config（应由 `rph_core/orchestrator.py` 统一处理）。
