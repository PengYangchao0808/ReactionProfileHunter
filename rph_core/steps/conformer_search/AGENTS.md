# rph_core/steps/conformer_search/AGENTS.md

## OVERVIEW
S1 的 Unified Conformer Engine：RDKit 构象生成 → CREST/xTB 采样与聚类(isostat) → DFT OPT/SP 耦合，选全局最低能构象。

## WHERE TO LOOK
- Engine: `engine.py`（`ConformerEngine`）
- QC facade: `rph_core/utils/qc_interface.py`（XTB/CREST/Gaussian factory）
- OPT/SP loop: `rph_core/utils/qc_task_runner.py`
- Geometry/log parsing: `rph_core/utils/geometry_tools.py`

## CONFIG SURFACES
- `config['step1']`: `crest`, `conformer_search`（window/rmsd/ngeom/isostat_*）
- `config['executables']`: `isostat`, `shermo`, `gaussian.wrapper_path`（默认 `./scripts/run_g16_worker.sh`）

## DIRECTORY LAYOUT
- `S1_ConfGeneration/<molecule_name>/xtb2/`：CREST/xTB 工件
- `S1_ConfGeneration/<molecule_name>/cluster/`：聚类/isostat
- `S1_ConfGeneration/<molecule_name>/dft/`：DFT OPT/SP（Gaussian/ORCA）

## GOTCHAS
- `engine.py` 里存在 `subprocess.run(..., shell=True, cwd=...)`（路径/转义敏感）；优先复用 `utils` 的 sandbox/toxic-path 逻辑。

## ANTI-PATTERNS
- 直接在这里硬编码 g16/orca/xtb 二进制路径（应走 `config['executables']` 或 utils 统一查找）。
- 复制 `test_results/` 行为作为逻辑分支（仅用于对照）。
