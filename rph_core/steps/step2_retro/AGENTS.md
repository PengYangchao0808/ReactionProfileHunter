# step2_retro/AGENTS.md

## OVERVIEW
从产物识别 forming bonds，生成：
- TS guess（拉伸到 TS 距离 + 受限优化）
- Reactant complex（拉伸到断裂距离 + 放松优化）

## WHERE TO LOOK
- Engine: `retro_scanner.py`
- SMARTS: `smarts_matcher.py`
- Stretch geometry: `bond_stretcher.py`

## MUST OUTPUTS
- `ts_guess.xyz`（S3 Berny 初始）
- `reactant_complex.xyz`（S3 QST2 + S4 畸变/片段）
- forming bonds indices（S4 fragment split 依据）

## CONVENTIONS
- 内部原子索引多为 0-index；给外部约束文件时转 1-index。
- product input 兼容 v2.1(文件) 与 v3.0/v6.1(S1_ConfGeneration 目录) 两种形态（retro_scanner 内有分支）。

## ANTI-PATTERNS
- 只产出 TS guess 不产出 reactant_complex
- 在这里做批量调度（应由 orchestrator/batch 层负责）
