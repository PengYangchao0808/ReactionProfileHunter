# rph_core/steps/conformer_search/AGENTS.md

## OVERVIEW
S1 的 Unified Conformer Engine (UCE) v3.1：支持两阶段构象搜索策略。

**两阶段工作流 (v3.1 新增)**:
1. **Stage 1**: GFN0-xTB 快速粗筛 → ISOSTAT 聚类
2. **Stage 2**: GFN2-xTB 精细优化 → ISOSTAT 聚类
3. **Stage 3**: DFT OPT-SP 耦合，选全局最低能构象

**单阶段模式 (向后兼容)**: 直接使用 GFN2-xTB

## WHERE TO LOOK
- Engine: `engine.py`（`ConformerEngine`）
- QC facade: `rph_core/utils/qc_interface.py`（XTB/CREST/Gaussian factory）
- OPT/SP loop: `rph_core/utils/qc_task_runner.py`
- Geometry/log parsing: `rph_core/utils/geometry_tools.py`

## CONFIG SURFACES (v3.1)
- `config['step1']['conformer_search']['two_stage_enabled']`: 主开关 (true/false)
- `config['step1']['conformer_search']['stage1_gfn0']`: Stage 1 配置
  - `gfn_level`: GFN 等级 (0=GFN0, 1=GFN1, 2=GFN2)
  - `energy_window_kcal`: 能量窗口 (推荐 10.0)
  - `clustering.run_after`: 是否运行聚类
- `config['step1']['conformer_search']['stage2_gfn2']`: Stage 2 配置
  - `gfn_level`: GFN 等级
  - `energy_window_kcal`: 能量窗口 (推荐 3.0)
- `config['step1']['crest']`: 单阶段模式配置 (向后兼容)
- `config['executables']`: `isostat`, `shermo`, `gaussian.wrapper_path`

## DIRECTORY LAYOUT

**两阶段模式**:
```
S1_ConfGeneration/<molecule_name>/
├── xtb2/
│   ├── stage1_gfn0/           # GFN0 粗筛结果
│   │   ├── crest_conformers.xyz
│   │   └── cluster/
│   │       ├── cluster.xyz    # Stage 2 输入
│   │       └── isostat.log
│   ├── stage2_gfn2/           # GFN2 精细结果
│   │   ├── crest_ensemble.xyz
│   │   └── cluster/
│   │       ├── cluster.xyz    # 最终系综 (传递给 DFT)
│   │       └── isostat.log
│   └── ensemble.xyz           # stage2/cluster/cluster.xyz 的副本
├── cluster/                   # 单阶段模式用 (遗留)
└── dft/                       # DFT OPT/SP (Gaussian/ORCA)
```

**单阶段模式**:
```
S1_ConfGeneration/<molecule_name>/
├── xtb2/
│   ├── crest_conformers.xyz
│   └── ensemble.xyz
├── cluster/
│   ├── isomers.xyz
│   └── cluster.xyz
└── dft/
```

## GOTCHAS
- `engine.py` 里存在 `subprocess.run(..., shell=True, cwd=...)`（路径/转义敏感）；优先复用 `utils` 的 sandbox/toxic-path 逻辑。
- 两阶段模式默认启用，如需单阶段模式设置 `two_stage_enabled: false`
- GFN0 比 GFN2 快 ~10 倍，但精度较低；两阶段策略先用 GFN0 覆盖构象空间，再用 GFN2 精细优化

## ANTI-PATTERNS
- 直接在这里硬编码 g16/orca/xtb 二进制路径（应走 `config['executables']` 或 utils 统一查找）。
- 复制 `test_results/` 行为作为逻辑分支（仅用于对照）。
