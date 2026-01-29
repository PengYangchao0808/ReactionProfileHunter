# ci/AGENTS.md

## OVERVIEW
CI 辅助脚本：以“导入链健康/导入风格”为主，避免重构后出现 ModuleNotFoundError。

## WHERE TO LOOK
- Import style checker: `check_imports.py`
- Integration guide: `README.md`
- Import rules (repo-wide): `IMPORT_GUIDELINES.md`

## CONVENTIONS
- 禁止多点相对导入（例如 `from ...utils`）；统一改为 `from rph_core.utils ...`。
- import smoke tests：`tests/test_imports_*.py`（用于 CI 早期失败）。

## ANTI-PATTERNS
- 通过增加 `sys.path` hack 来“修复”导入（应改 import 语句本身）。
- 在 `rph_core_backup_20260115/` 里修 import 再整体替换（diff 风险高）。
