# step4_features/AGENTS.md

## OVERVIEW
从 S1/S2/S3 的结构与 QC 输出抽取几何/电子/热化学特征；包含片段畸变/相互作用（fragment-based）逻辑。

## WHERE TO LOOK
- Orchestrator: `feature_miner.py`
- Data plumbing: `context.py`, `status.py`, `schema.py`, `path_accessor.py`
- Distortion: `distortion_calculator.py`
- Fragments: `fragment_extractor.py`
- Parsing: `log_parser.py`
- Extractors: `extractors/`（geometry/thermo/nics/nbo_e2/interaction_analysis/qc_checks）

## CONVENTIONS
- extractor 不应硬编码目录/文件名：统一从 context/path accessor 取路径。
- 对缺失 artifact 的行为应"可降级但可追溯"：产出部分 features，同时在 status/log 记录缺失原因。
- forming bonds（来自 S2）是 fragment split 的主线索，避免各处各算一遍。

## ANTI-PATTERNS
- 在多个 extractor 内复制 fragment split 逻辑（应集中在 fragment_extractor）
- 缺失 fchk/log 时无记录地"成功"
