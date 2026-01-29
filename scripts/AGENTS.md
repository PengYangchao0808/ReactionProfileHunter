# scripts/AGENTS.md

## OVERVIEW
外部程序 wrapper（主要是 Gaussian g16）：处理 scratch 隔离、低磁盘保护、环境探测与自清理。

## WHERE TO LOOK
- Gaussian worker: `run_g16_worker.sh`

## INTEGRATION POINTS
- 默认 wrapper 路径：`./scripts/run_g16_worker.sh`
- 引用位置（示例）：
  - `rph_core/utils/qc_interface.py`
  - `rph_core/steps/step3_opt/qst2_rescue.py`
  - `rph_core/steps/step3_opt/irc_driver.py`
  - `rph_core/steps/conformer_search/engine.py`

## ANTI-PATTERNS
- 在 Python 侧自行清理/复用同一 GAUSS_SCRDIR（wrapper 已按 job 隔离并清理）。
- 在 wrapper 内写死机器相关路径（保持 config/PATH 探测优先级）。
