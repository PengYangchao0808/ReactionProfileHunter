# config/AGENTS.md

## OVERVIEW
defaults + QC 输入模板（Gaussian .gjf/.com）。

## WHERE TO LOOK
- `defaults.yaml`
- `templates/gaussian_*.gjf` / `templates/gaussian_*.com`

## CONVENTIONS
- defaults.yaml 是配置单一来源；代码应通过统一 config 传递获取，而非散落默认值。
- 模板修改需考虑可追溯性与复用策略（例如 oldchk/Guess=Read 若启用）。

## ANTI-PATTERNS
- 在代码中复制模板字符串（应读 templates 文件）
- 维护多份 defaults 并让行为分叉（仓库存在 .bak：避免逻辑依赖它）
