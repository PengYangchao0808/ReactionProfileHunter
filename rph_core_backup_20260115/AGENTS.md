# rph_core_backup_20260115/AGENTS.md

## OVERVIEW
历史冻结备份：用于对照/定位回归/参考回滚策略，不作为当前开发入口。

## USAGE
- ✅ 对照差异：从 backup 找"以前怎么做"
- ❌ 不在此目录修 bug；修复应落在 `rph_core/`

## ANTI-PATTERNS
- 在 backup 里修完再整体替换主目录（高风险、diff 不可控）
