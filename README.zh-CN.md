# ReactionProfileHunter

[English](README.md) | [中文](README.zh-CN.md)

<div align="center">

**产物驱动的 DFT 反应机理探索流水线 (S0-S4)**

[![Version](https://img.shields.io/badge/version-4.0.0-blue.svg)](https://github.com/PengYangchao0808/ReactionProfileHunter)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Code Style](https://img.shields.io/badge/code%20style-agents%20ready-success.svg)](AGENTS.md)

</div>

> **项目统计：** ~50k 行代码，198 个 Python 文件 | **开发指南：** 参见 [AGENTS.md](AGENTS.md) 了解代码规范

---

## 项目简介

ReactionProfileHunter (RPH) 是一个**产物驱动**的自动化 DFT 反应机理探索流水线（S0-S4），专注于有机反应的过渡态搜索与几何优化。特征提取与 ML 训练由配套的 `RPH_Postprocess` 包独立完成。

### 核心功能

- **产物驱动策略**：从产物分子出发，自动逆向搜索反应路径
- **五阶段 DFT 流水线**：机理验证 -> CENSO-LITE 构象搜索 -> PEB 逆向扫描 -> 低级别计算 -> 高精度计算
- **固定理论协议**：CREST/GFN2 + B97-3c SP (S1)、ORCA B97-3c -> r2SCAN-3c (S3)、ORCA M062X -> ORCA wB97M-V (S4)
- **多引擎支持**：Gaussian、ORCA、xTB、CREST
- **基于清单的断点恢复**：每结构失败隔离
- **WSL 实时状态查看**（`rph_watch`）
- **仅数据集入口**：不支持 SMILES 模式

### 使用场景

- 有机反应机理研究（环加成、重排、取代等）
- 过渡态结构预测与验证
- 反应数据库构建与高通量筛选
- 反应性质的机器学习研究（通过 RPH_Postprocess）

---

## 主要特性

### V4 架构 — 固定 S0-S4 流水线

V4 流水线以清晰、可重复的协议取代了 V3 的 Berny/QST2/IRC 救援链和 SMILES 入口：

| 阶段 | 协议 | 详情 |
|------|------|------|
| S0 | 可信数据集机理验证 | 读取 Reaxys 清洗后的 CSV；解析成键、原子映射、反应类型 |
| S1 | CENSO-LITE 构象搜索 | CREST/GFN2 采样 + ORCA B97-3c SP 排序 + xTB mRRHO。无 DFT OPT/FREQ。 |
| S2 | PEB 逆向扫描 | 从 S1 选定构象出发的 xTB 势能面成键扫描；生成 TS 初猜和中间体 |
| S3 | 低级别计算 | ORCA B97-3c OPT/OptTS + Freq；ORCA r2SCAN-3c SP（CPCM 丙酮） |
| S4 | 高精度计算 | ORCA M062X OPT/OptTS + Freq；ORCA wB97M-V SP（CPCM 丙酮） |

流水线流程：

```
CSV 反应记录 -> S0（机理） -> S1（CENSO-LITE） -> S2（PEB）
                                                      |
                                                      v
                                      S3（B97-3c OPT + r2SCAN-3c SP）
                                                      |
                                                      v
                                      S4（M062X OPT + wB97M-V SP）
                                                      |
                                                      v
                                      RPH_Postprocess（外部）
                                        -> 特征提取、ML 数据集
```

### 已移除的 V3 功能

以下 V3 能力不属于 V4：

- SMILES 模式入口（已移除；仅接受 CSV 数据集）
- Berny / QST2 救援 / IRC 验证链
- S1 DFT OPT 或 DFT SP（S1 仅做 CENSO-LITE）
- 双阶段构象漏斗（GFN0 -> GFN2 -> DFT）
- 路易斯酸校准子系统
- `reaction_profiles` 配置段
- DR 聚合器和 condition_feature_merger
- `QCTaskRunner` 和旧版 `CheckpointManager`
- S4 作为特征提取（现在是高精度 QC）

### 技术亮点

- **分子自治架构**：每个分子独立管理其计算目录（S1_ConfSearch/[molecule]/）
- **每结构失败隔离**：S4 任务失败不会删除或覆盖可用的 S3 输出
- **基于清单的交接**：每个阶段写入版本化清单；下一阶段读取它
- **沙盒隔离执行**：每个量子化学任务在独立沙盒中运行
- **`pipeline.state` 驱动的可恢复性**：重新运行相同命令即可继续

---

## 项目结构

```
ReactionProfileHunter/
├── rph_core/              # 核心源码 — 仅使用绝对导入
│   ├── v4_orchestrator.py # V4 流水线编排 + CLI（V4Orchestrator）
│   ├── __main__.py        # python -m rph_core 入口
│   ├── v4_watch.py        # WSL 实时状态查看器
│   ├── scheduling/        # 批处理调度与任务编排
│   ├── steps/             # S0-S4 步骤实现
│   │   ├── conformer_search/  # S1 CENSO-LITE 引擎
│   │   ├── mechanism_classifier/  # S0 机理验证
│   │   ├── step2_retro/       # S2 PEB 扫描器
│   │   ├── step3_lowlevel/    # S3 低级别 QC（B97-3c / r2SCAN-3c）
│   │   ├── step4_highlevel/   # S4 高精度 QC（M062X / wB97M-V）
│   │   └── stage_calculator.py  # S3/S4 通用的 OPT+SP+FREQ 分发
│   └── utils/             # QC/IO/检查点基础设施
│       ├── qc_interface.py     # 所有 QC 调用均经由此处
│       ├── v4_checkpoint.py    # V4 基于清单的检查点
│       ├── qc_jobs.py          # QC 作业规范映射
│       └── qc_models.py        # QC 数据模型
├── bin/                  # CLI 封装（sys.path + main()）
│   ├── rph_run           # V4 主入口
│   ├── rph_watch         # WSL 实时状态查看器
│   └── rph               # 短别名
├── config/               # 运行时配置（唯一数据源）
│   ├── defaults.yaml     # 所有配置键
│   └── templates/        # Gaussian .gjf/.com 模板
├── tests/                # pytest 测试套件（约 82 个文件）
│   ├── conftest.py       # 将仓库根目录加入 sys.path
│   └── fixtures/         # 静态测试数据
├── scripts/              # 工具脚本与 CI 工具
│   ├── ci/check_imports.py   # 导入风格检查门控
│   └── *.py / *.sh           # 数据构建、分析
├── rph_benchmark/        # 独立评估套件
│   ├── confsearch/       # 构象搜索协议基准测试
│   └── dft_theory/       # DFT 方法基准测试
├── docs/                 # 设计文档与分析报告
├── AGENTS.md             # V4 智能体编码指南
├── README.md             # 英文版
└── README.zh-CN.md       # 中文版
```

---

## 快速开始

### 系统要求

| 组件 | 要求 |
|------|------|
| 操作系统 | Linux（推荐 Ubuntu 20.04+ 或 CentOS 7+） |
| Python | 3.8 或更高版本 |
| 内存 | 至少 16 GB（推荐 64 GB） |
| CPU | 多核处理器（推荐 16 核以上） |

### 依赖的量子化学软件

V4 完整运行 S0-S4 需要以下全部软件：

- ORCA — S1 B97-3c 排序以及全部 S3/S4 DFT 任务
- xTB — S1 mRRHO 热化学、S2 PEB 扫描
- CREST — S1 构象采样

### 安装

1. **克隆仓库**
   ```bash
   git clone https://github.com/PengYangchao0808/ReactionProfileHunter.git
   cd ReactionProfileHunter
   ```

2. **安装 Python 依赖（若已配置打包文件）**
   ```bash
   pip install -e .
   ```

3. **配置量子化学软件路径**

   编辑 `config/defaults.yaml`，填入已安装软件的路径：
   ```yaml
   executables:
      orca:
       path: "/path/to/orca/orca"
       ld_library_path: "/path/to/orca"
     xtb:
       path: "/path/to/xtb/bin/xtb"
     crest:
       path: "/path/to/crest/crest"
   ```

4. **验证安装**
   ```bash
   python -m pytest -q
   ```

### 基础使用

#### 命令行入口（仅数据集模式）

V4 不接受 SMILES 输入。必须提供经过清洗的 Reaxys CSV 数据集：

```bash
bin/rph_run --csv data/trusted_reactions.csv --rx-id RXN_000001 --output ./Output/RXN_000001
```

#### 部分流水线运行

使用 `--stop-after` 仅运行选定阶段：

```bash
# 仅运行 S0 到 S1（构象搜索）
bin/rph_run --csv data/trusted_reactions.csv --rx-id RXN_000001 --output ./Output/RXN_000001 --stop-after s1

# 运行到 S3（跳过高精度 S4）
bin/rph_run --csv data/trusted_reactions.csv --rx-id RXN_000001 --output ./Output/RXN_000001 --stop-after s3
```

恢复运行默认采用严格策略。如果已完成阶段与当前科学配置不一致，RPH 会在启动
QC 前停止并报告差异，需要显式选择后续动作：

```bash
# 从 S1 重新计算，并使 S2-S4 失效
bin/rph_run --csv data/trusted_reactions.csv --rx-id RXN_000001 --output ./Output/RXN_000001 --stop-after s2 --recompute-from s1

# 保留已完成的 S0/S1 产物，从 S2 开始新计算
bin/rph_run --csv data/trusted_reactions.csv --rx-id RXN_000001 --output ./Output/RXN_000001 --stop-after s2 --resume-policy use-existing-upstream --start-from s2
```

首次运行会将配置冻结到 `run.config.json`。所有可恢复的阶段及变体状态统一存储在
`pipeline.state`；旧 `.rph/checkpoint.json` 仅做只读迁移。

#### WSL 实时状态查看

在另一个终端中监控运行进度：

```bash
bin/rph_watch --output ./Output/RXN_000001 --watch
```

#### 自定义配置

```bash
bin/rph_run --csv data/trusted_reactions.csv --rx-id RXN_000001 --output ./Output/RXN_000001 --config config/custom.yaml
```

#### Python API

```python
from pathlib import Path
from rph_core.v4_orchestrator import V4Orchestrator
from rph_core.steps.mechanism_classifier.s0_record import load_s0_reaction_record

record = load_s0_reaction_record(Path("data/trusted_reactions.csv"), "RXN_000001")
orchestrator = V4Orchestrator()
result = orchestrator.run(record, Path("./Output/RXN_000001"), stop_after="s4")
if result.get("success"):
    print(f"流水线完成。清单: {result}")
```

#### 入口方式汇总

| 方式 | 命令 | 说明 |
|------|------|------|
| CLI 脚本 | `bin/rph_run --csv <file> --rx-id <id> --output <dir>` | 自动将仓库根目录加入 sys.path |
| 短别名 | `bin/rph --csv <file> --rx-id <id> --output <dir>` | 与 rph_run 功能相同 |
| Python 模块 | `python -m rph_core --csv <file> --rx-id <id> --output <dir>` | 相同参数，依赖 Python 路径 |
| 程序化调用 | `V4Orchestrator().run()` | 完全控制，返回 dict |

---

## 详细文档

### 架构概述

V4 是基于清单交接的串行流水线：

```
  CSV 行
     |
     v
S0_Mechanism/mechanism.json -- forming_bonds、反应类型、原子映射
     |
     v
S1_ConfSearch/product/manifest.json -- selected.xyz、构象候选、能量
     |
     v
S2_PEB/manifest.json -- peb_peak（TS 初猜）、intermediate_seed、扫描轮廓
     |
     v
S3_LowLevel/manifest.json -- 每结构 opt、freq、sp 结果
     |
     v
S4_HighLevel/manifest.json -- 每结构高精度 opt、freq、sp 结果
     |
     v
pipeline.result.json -- 最终摘要
```

各阶段是幂等的。每个阶段写入其清单和签名。重新运行相同命令会复用有效检查点，仅重新执行未完成的阶段。

### 阶段间交接

| 来源 | 目标 | 产物 |
|------|------|------|
| S0 记录 | S1 | forming_bonds（映射空间）、reaction_type |
| S1 清单 | S2 | `selected.xyz`、构象能量 |
| S2 清单 | S3 | `peb_peak`（TS 初猜）、`intermediate_seed`、forming_bonds（XYZ 索引） |
| S3 清单 | S4 | 每结构 opt/sp/freq 结果、opt_xyz、状态、usable_for_ml |
| S4 清单 | RPH_Postprocess | 高精度能量、优化几何、频率数据 |

### 理论协议

```
S1: CREST/GFN2 构象采样
    + ORCA B97-3c 单点排序
    + xTB GFN1 mRRHO 热化学
    无 DFT 优化或 DFT 单点用于几何精化。

S3: ORCA B97-3c OPT（最小值）/ OptTS（TS）+ 独立 Freq（仅 TS）
    -> ORCA r2SCAN-3c SP
    溶剂: CPCM(丙酮)

S4: ORCA M062X/def2-SVP OPT / OptTS + 独立 Freq（仅 TS）
    -> ORCA wB97M-V/def2-TZVPP SP
    溶剂: CPCM(丙酮)
```

所有方法、路径、资源和超时均在 `config/defaults.yaml` 的 `theory.s3_low_level` 和 `theory.s4_high_precision` 中定义。

### 输出契约

流水线生成以下权威文件：

```
Output/<rx_id>/
├── S0_Mechanism/mechanism.json
├── S1_ConfSearch/product/manifest.json
├── S2_PEB/manifest.json
├── S3_LowLevel/manifest.json
├── S4_HighLevel/manifest.json
├── pipeline.state
├── pipeline.result.json
├── rph_v4.log
```

`forming_bonds` 在 S0 清单中是权威的，使用基于 0 的 XYZ 索引。每个 S3/S4 结构记录保留输入几何、优化几何（如有）、输出文件、能量、状态和 `usable_for_ml`。

### 输出目录示例

```
Output/RXN_000001/
├── .rph_schema.json
├── pipeline.result.json
├── pipeline.state
├── rph_v4.log
├── S0_Mechanism/
│   ├── mechanism.json
├── S1_ConfSearch/
│   ├── product/
│   │   ├── manifest.json
│   │   ├── selected.xyz
│   │   └── raw_censo/
│   └── manifest.json
├── S2_PEB/
│   ├── manifest.json
│   ├── scan_profile.json
│   ├── ts_guess.xyz
│   └── intermediate.xyz
├── S3_LowLevel/
│   ├── manifest.json
│   ├── status.json
│   ├── events.jsonl
│   ├── product/
│   │   ├── opt/
│   │   └── sp/
│   ├── product_int/
│   │   ├── opt/
│   │   └── sp/
│   └── product_ts/
│       ├── opt_ts/
│       ├── freq/
│       └── sp/
├── S4_HighLevel/
│   ├── manifest.json
│   ├── status.json
│   ├── events.jsonl
│   ├── s4.log
│   ├── product/
│   │   ├── opt/
│   │   └── sp/
│   ├── product_int/
│   │   ├── opt/
│   │   └── sp/
│   └── product_ts/
│       ├── opt_ts/
│       ├── freq/
│       └── sp/
```

### 特征提取（外部）

特征提取不属于 RPH V4 流水线。由配套的 `RPH_Postprocess` 包处理，它消费 S4 清单并生成 ML 就绪的特征数据集。

在 RPH 完成后运行：

```bash
rph-features extract --rph-run ./Output/RXN_000001 --output ./Output/RXN_000001/S4_Data
```

RPH_Postprocess 生成：
- `features_raw.csv` - 所有提取的特征
- `features_mlr.csv` - ML 就绪特征矩阵
- `feature_meta.json` - 溯源和配置元数据

这种分离使 RPH 专注于 DFT 计算，而 RPH_Postprocess 处理所有后期分析。

---

## 配置指南

### QC 可执行文件

```yaml
executables:
  gaussian:
    path: "/opt/software/gaussian/g16/g16"
    root: "/opt/software/gaussian/g16"
  orca:
    path: "/opt/software/orca/orca"
    ld_library_path: "/opt/software/orca"
  crest:
    path: "/opt/software/crest/crest"
  xtb:
    path: "/opt/software/xtb/bin/xtb"
```

### 计算资源

```yaml
resources:
  mem: "32GB"
  nproc: 16
  orca_maxcore_safety: 0.65
```

### S3 低级别理论

```yaml
theory:
  s3_low_level:
    optimization:
      engine: orca
      method: B97-3c
      solvent: acetone
      solvent_model: CPCM
      route_minimum: "Opt"
      route_ts: "OptTS"
      frequency:
        enabled_for_ts: true
        imaginary_cutoff_cm1: -50.0
    single_point:
      engine: orca
      method: r2SCAN-3c
      solvent: acetone
      solvent_model: CPCM
```

### S4 高精度理论

```yaml
theory:
  s4_high_precision:
    optimization:
      engine: orca
      method: M062X
      basis: def2-SVP
      aux_basis: def2/J
      solvent: acetone
      solvent_model: CPCM
      route_minimum: Opt
      route_ts: OptTS
      grid: DefGrid3
      scf: TightSCF
      frequency:
        enabled_for_ts: true
        imaginary_cutoff_cm1: -50.0
    single_point:
      engine: orca
      method: wB97M-V
      basis: def2-TZVPP
      aux_basis: def2/J
      solvent: acetone
      solvent_model: CPCM
```

### S1 CENSO-LITE

```yaml
step1:
  protocol: censo_lite
  censo_lite:
    crest:
      gfn_level: 2
      search_mode: imtd_smtd
      energy_window_kcal: 6.0
      solvent: acetone
    ranking:
      engine: orca
      method: B97-3c
      solvent: acetone
    xtb_thermo:
      enabled: true
      gfn_level: 1
      temperature_k: 298.15
    deduplication:
      backend: torsion_signature
      heavy_atom_rmsd_prefilter_A: 0.25
      torsion_bin_deg: 20.0
```

### PEB 扫描

```yaml
step2:
  scan:
    topology_guard_enabled: true
    scan_start_distance: 4.0
    scan_end_distance: 1.5
    scan_steps: 24
    retry_on_boundary_maximum: true
```

### 配置文件位置

- 主配置：`config/defaults.yaml`
- 模板：`config/templates/`

---

## 测试

### 运行测试

```bash
# 运行全部测试
pytest -v tests/

# V4 协议契约测试
pytest tests/test_v4_protocol_contract.py -v

# V4 检查点测试
pytest tests/test_v4_checkpoint.py -v

# V4 阶段计算器测试
pytest tests/test_v4_stage_calculator.py -v

# 测试覆盖率
pytest --cov=rph_core --cov-report=html
```

### 导入风格检查（CI 门控）

```bash
python scripts/ci/check_imports.py rph_core
```

### 注意事项

- `tests/conftest.py` 自动将仓库根目录添加到 `sys.path`，无需可编辑安装即可运行测试
- 集成测试使用 mock 量子化学计算，不需要真实的 ORCA
- V3 测试已归档到 `tests/deprecated_v3/`
- 参见 [AGENTS.md](AGENTS.md) 了解验证命令

---

## 故障排查

### 1. 可选 Gaussian 兼容后端找不到可执行文件

受支持的 V4 S0-S4 流程不再使用 Gaussian；本节仅适用于直接调用保留的兼容接口。

**错误信息**：
```
FileNotFoundError: Gaussian executable not found: /opt/software/gaussian/g16/g16
```

**解决方案**：
- 检查 `config/defaults.yaml` 中的 `executables.gaussian.path` 是否正确
- 确认已安装 Gaussian 并设置环境变量
- 尝试直接运行 `g16 < test.gjf` 测试 Gaussian 是否可用

### 2. 内存不足导致计算崩溃

**错误信息**：
```
Error termination via Lnk1e in /root/g16/g16
```

**解决方案**：
- 降低 `resources.mem` 和 `resources.nproc`
- 降低 `orca_maxcore_safety`
- 减少 `step2.scan` 中的扫描步数

### 3. CENSO-LITE 构象搜索失败

**错误信息**：
```
RuntimeError: S1 manifest has no selected candidate
```

**解决方案**：
- 增加 `step1.censo_lite.crest` 中的 CREST 采样时间或能量窗口
- 检查 ORCA 和 xTB 可执行文件是否正确配置
- 验证输入 CSV 行的产物 SMILES 是否有效

### 4. S3/S4 特定结构优化失败

当 `continue_on_structure_failure` 启用时（默认），流水线会继续处理其余结构。检查清单中的每结构状态：

- 如果 `opt_status` 不是 `complete`，检查 opt 目录中的 ORCA 输出
- 如果 `sp_status` 不是 `complete`，检查单点输出
- 如果 `ts_frequency_valid` 是 `false`，检查频率输出中的虚频模式

### 5. 流水线未正确恢复

删除 `pipeline.state` 和目标阶段清单，然后重新运行：

```bash
rm -f Output/RXN_000001/pipeline.state
bin/rph_run --csv data/trusted_reactions.csv --rx-id RXN_000001 --output ./Output/RXN_000001
```

### 日志文件

V4 将日志写入输出目录中的 `rph_v4.log`：

```bash
tail -f Output/RXN_000001/rph_v4.log
```

### 调试模式

```bash
bin/rph_run --csv data/trusted_reactions.csv --rx-id RXN_000001 --output ./Output/RXN_000001 --log-level DEBUG
```

### 获取帮助

1. 查看 [GitHub Issues](https://github.com/PengYangchao0808/ReactionProfileHunter/issues)
2. 附上完整错误日志和配置文件
3. 提供最小可复现示例

---

## 文档索引

| 文档 | 内容 |
|------|------|
| [`AGENTS.md`](AGENTS.md) | V4 编程指南 — 构建/测试命令、代码规范、约定 |
| [`rph_core/AGENTS.md`](rph_core/AGENTS.md) | 核心包架构与代码地图 |
| [`rph_core/steps/AGENTS.md`](rph_core/steps/AGENTS.md) | 步骤架构与输出契约 |
| [`rph_core/steps/conformer_search/AGENTS.md`](rph_core/steps/conformer_search/AGENTS.md) | CENSO-LITE 构象搜索 |
| [`rph_core/steps/step2_retro/AGENTS.md`](rph_core/steps/step2_retro/AGENTS.md) | PEB 扫描器 |
| [`rph_core/steps/step3_lowlevel/AGENTS.md`](rph_core/steps/step3_lowlevel/AGENTS.md) | S3 低级别 QC |
| [`rph_core/steps/step4_highlevel/AGENTS.md`](rph_core/steps/step4_highlevel/AGENTS.md) | S4 高精度 QC |
| [`rph_core/utils/AGENTS.md`](rph_core/utils/AGENTS.md) | QC/IO/检查点工具参考 |
| [`config/AGENTS.md`](config/AGENTS.md) | 配置结构说明 |
| [`tests/AGENTS.md`](tests/AGENTS.md) | 测试组织与约定 |
| [`scripts/AGENTS.md`](scripts/AGENTS.md) | 脚本与 CI 工具说明 |
| [`rph_benchmark/AGENTS.md`](rph_benchmark/AGENTS.md) | 基准测试套件文档 |
| [`docs/RPH_V4_COMPLETION_MASTER_PLAN.md`](docs/RPH_V4_COMPLETION_MASTER_PLAN.md) | V4 架构主计划 |
| [`docs/WSL_TEST_PLAN_V4.md`](docs/WSL_TEST_PLAN_V4.md) | WSL 测试计划与 T0-T6 关卡 |

---

## 学术引用

如果您在研究中使用 ReactionProfileHunter，请引用：

```bibtex
@software{reactionprofilehunter2025,
  title = {ReactionProfileHunter: Automated Reaction Mechanism Exploration},
  author = {Peng Yangchao},
  year = {2026},
  url = {https://github.com/PengYangchao0808/ReactionProfileHunter},
  version = {4.0.0}
}
```

### 相关资源

- **理论背景**：
  - Houk 课题组的计算化学方法论
  - Grimme 的 GFN-xTB 方法论文
  - wB97M-V 泛函性能评估

- **工具文档**：
  - [Gaussian 16 用户手册](https://gaussian.com/man/)
  - [ORCA 输入库](https://sites.google.com/site/orcainputlibrary/)
  - [xTB 文档](https://xtb-docs.readthedocs.io/)

---

## 贡献指南

1. **Fork 本仓库**
2. **创建特性分支**：`git checkout -b feature/amazing-feature`
3. **提交更改**：`git commit -m 'Add amazing feature'`
4. **推送到分支**：`git push origin feature/amazing-feature`
5. **提交 Pull Request**

### 代码规范

- 使用 **pathlib.Path** 处理所有路径（禁止使用字符串路径）
- 使用 **logging.getLogger(__name__) 或 LoggerMixin** 进行日志记录（核心代码中禁止 `print`）
- 遵循 **输出目录幂等性**：存在可复用输出时跳过重算
- 所有量子化学调用必须通过 `utils/qc_interface.py`
- 所有外部 QC 任务必须通过 `qc_jobs.py` 和现有的 QC 接口路由
- **仅使用绝对导入** — 禁止多点多级相对导入（例如使用 `from rph_core.utils...` 而非 `from ...utils`）
- 不得在阶段模块中添加 `subprocess.run`
- 任务降级时保留 XYZ、输出日志和失败诊断信息
- S4 任务失败不得删除或替换可用的 S3 输出
- 每个阶段写入版本化清单并可通过 `pipeline.state` 恢复
- 详见 [AGENTS.md](AGENTS.md) 了解完整编码规范

### 测试要求

- 新功能必须包含单元测试
- 所有测试必须通过：`python -m pytest -q`
- 保持测试覆盖率 > 80%

### 文档要求

- 更新相关的 `AGENTS.md` 文件
- 为新功能添加使用示例
- 更新 `rph_core/version.py` 和 README 徽章中的版本号

---

## 许可证

本项目采用 [MIT License](LICENSE) 开源协议。

---

## 联系方式

- **作者**：Peng Yangchao
- **项目主页**：https://github.com/PengYangchao0808/ReactionProfileHunter
- **问题反馈**：https://github.com/PengYangchao0808/ReactionProfileHunter/issues

---

## 致谢

- **Houk 教授**：提供双层级计算策略建议
- **Grimme 课题组**：开发 xTB 和 CREST 工具
- **Gaussian 和 ORCA 开发团队**：提供优秀的量子化学软件
- **社区贡献者**：感谢所有提交 issue 和 PR 的朋友

---

<div align="center">

**如果觉得本项目有帮助，请给我们一个 Star！**

</div>
