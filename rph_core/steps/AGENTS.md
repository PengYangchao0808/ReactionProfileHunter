# rph_core/steps/AGENTS.md

## OVERVIEW
步骤化实现：S1 锚定/构象 → S2 逆向扫描 → S3 TS 优化/救援 → S4 特征提取。

## WHERE TO LOOK
- S2: `step2_retro/`
- S3: `step3_opt/`
- S4: `step4_features/`
- Anchor/conformer: `anchor/`, `conformer_search/`

## OUTPUT CONTRACT
- S1 MUST: `product/product_min.xyz` (and optionally `precursor/precursor_min.xyz`) in `S1_ConfGeneration/`
- S2 MUST: `ts_guess.xyz` + `reactant_complex.xyz` (+ forming bonds indices)
- S3 MUST: `ts_final.xyz`（并保留可追溯 logs/chk）
- S4 MUST: `features_raw.csv`, `features_mlr.csv`, `feature_meta.json` 或等价结构化产物

## ANTI-PATTERNS
- 缺失 reactant_complex 仍继续进入 S3/S4（应 fail-fast 或补救）
- 目录结构变更不更新 orchestrator/path accessor
