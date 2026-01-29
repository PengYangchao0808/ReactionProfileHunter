# rph_core/steps/anchor/AGENTS.md

## OVERVIEW
S1 锚定阶段：对底物池与产物做构象搜索 + OPT/SP，产出每个分子的全局最低能结构与 SP 能量。

## WHERE TO LOOK
- Coordinator: `handler.py`（`AnchorPhase.run()`）
- Conformer search engine: `rph_core/steps/conformer_search/engine.py`

## IO / OUTPUTS
- 输出写入 `S1_Product/<molecule_name>/`（由 `ConformerEngine` 管理 `xtb2/`, `cluster/`, `dft/`）
- `anchored_molecules[name]` 结构：`xyz`(Path) + `e_sp`(float) + 可选 `log/chk/fchk`

## CONVENTIONS
- 目录统一用 `pathlib.Path`；不要在此层自行拼接相对路径。
- Gaussian `.chk` -> `.fchk` 转换复用 `rph_core/utils/qc_interface.py:try_formchk()`。

## ANTI-PATTERNS
- 在 S1 内重复实现 QC runner / sandbox（应复用 `rph_core/utils/qc_interface.py`）。
- 产物/底物输出命名漂移导致 S2/S3/S4 找不到输入（改名需同步 orchestrator/contract）。
