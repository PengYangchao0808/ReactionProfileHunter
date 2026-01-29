# step3_opt/AGENTS.md

## OVERVIEW
TS 优化：Berny 为主；失败触发 QST2 rescue；可选 IRC；validator 负责收敛/虚频一致性校验。

## WHERE TO LOOK
- Coordinator: `ts_optimizer.py`
- Berny: `berny_driver.py`
- QST2: `qst2_rescue.py`
- IRC: `irc_driver.py`
- Validation: `validator.py`
- Rescue levels: `rph_core/utils/oscillation_detector.py`

## CONTRACT
Inputs:
- `ts_guess.xyz`（来自 S2）
- `reactant_complex.xyz`（来自 S2）+ product endpoint（来自 S1，用于 QST2）
Output:
- `ts_final.xyz`（validated）

## ANTI-PATTERNS
- Berny 失败后静默继续（必须 rescue 或明确失败）
- 不保留 log/chk 导致不可复现/不可做下游特征
