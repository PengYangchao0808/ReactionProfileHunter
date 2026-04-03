# ReactionProfileHunter

[English](README.md) | [中文](README.zh-CN.md)

<div align="center">

**产物驱动的反应机理探索与特征提取工具**

[![Version](https://img.shields.io/badge/version-6.2.0-blue.svg)](https://github.com/yourusername/ReactionProfileHunter)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Code Style](https://img.shields.io/badge/code%20style-agents%20ready-success.svg)](AGENTS.md)

</div>

> **项目统计：** ~50k 行代码，198 个 Python 文件 | **开发指南：** 参见 [AGENTS.md](AGENTS.md) 了解代码规范

---

## 📋 项目简介

ReactionProfileHunter (RPH) 是一个**产物驱动**的自动化反应机理探索流水线，专为有机反应的过渡态搜索、几何优化和特征提取而设计。

### 核心功能

- 🎯 **产物驱动策略**：从产物分子出发，自动逆向搜索反应路径
- 🔄 **四步流水线**：锚定/构象搜索 → 逆向扫描 → 过渡态优化 → 特征提取
- ⚡ **智能双层级计算**：xTB 预优化 → B3LYP 几何优化 → wB97X-D3BJ 高精度能量
- 🔧 **多引擎支持**：Gaussian、ORCA、xTB、CREST、Multiwfn 无缝集成
- 📊 **特征工程**：自动提取反应能垒、几何参数、电子结构等机器学习特征
- 💾 **断点恢复**：支持计算中断后从检查点继续运行

### 使用场景

- 有机反应机理研究（环加成、重排、取代等）
- 过渡态结构预测与验证
- 反应数据库构建与高通量筛选
- 反应性质的机器学习研究

---

## ✨ 主要特性

### v6.2.0 新特性 (机理感知特征增强)

- ✅ **Step1 活化特征**：热力学加权能量 (Boltzmann/Gibbs)、构象熵/柔性、离去基团几何、α-H 门控参数
- ✅ **Step2 环化特征**：动力学/热化学、TS 几何 (forming bonds)、CDFT 指标 (Fukui/Dual Descriptor/QTAIM)
- ✅ **正向扫描过渡态搜索 (NEW)**：基于 xTB 原生 `$scan` 实现，支持任意环加成反应 ([4+2], [4+3], [3+2], [5+2])
- ✅ **反应配置 profiles (NEW)**：按反应类型配置扫描参数
- ✅ **GEDT/CDFT 增强**：基于 Forming Bonds 的图分割电荷转移计算、单位锁定 (eV)、范围校验
- ✅ **Multiwfn 集成**：Tier-1 支持 (Fukui/Dual Descriptor/QTAIM)，非交互批处理，失败不阻断，缓存机制
- ✅ **严格合同**：Feature Meta 记录完整配置快照与运行溯源，缺失特征显式降级 (NaN + Warning)

### v6.1.0 特性回顾

- ✅ **简化 S4 模块**：仅负责特征提取，不再执行量子化学计算
- ✅ **固定输出格式**：`features_raw.csv`、`features_mlr.csv`、`feature_meta.json`
- ✅ **NBO 默认关闭**：减少计算开销，可按需启用
- ✅ **改进错误处理**：更清晰的日志和异常信息

### 技术亮点

- **分子自治架构**：每个分子独立管理其计算目录（`S1_ConfGeneration/[Molecule]/`）
- **OPT-SP 耦合循环**：几何优化与单点能计算智能耦合，提高收敛成功率
- **路径兼容性**：自动处理不同版本的目录结构，保证向后兼容
- **沙盒隔离执行**：每个量子化学任务在独立沙盒中运行，避免文件冲突

---

## 🚀 快速开始

### 系统要求

| 组件 | 要求 |
|------|------|
| **操作系统** | Linux (推荐 Ubuntu 20.04+, CentOS 7+) |
| **Python** | 3.8 或更高版本 |
| **内存** | 至少 16 GB（推荐 64 GB） |
| **CPU** | 多核处理器（推荐 16 核以上） |

### 依赖的量子化学软件

至少需要安装以下**一组**软件：

**选项 1：Gaussian 方案**
- [Gaussian 16](https://gaussian.com/) - 几何优化
- [xTB](https://github.com/grimme-lab/xtb) - 预优化（可选但推荐）
- [CREST](https://github.com/crest-lab/crest) - 构象搜索（可选）

**选项 2：ORCA 方案**
- [ORCA](https://orcaforum.kofo.mpg.de/) - 几何优化与单点能
- [xTB](https://github.com/grimme-lab/xtb) - 预优化（可选但推荐）
- [CREST](https://github.com/crest-lab/crest) - 构象搜索（可选）

**选项 3：混合方案（推荐）**
- Gaussian 16 - 几何优化
- ORCA - 高精度单点能计算
- xTB + CREST - 预优化与构象搜索

### 安装

1. **克隆仓库**
   ```bash
   git clone https://github.com/yourusername/ReactionProfileHunter.git
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
     gaussian:
       path: "/path/to/g16/g16"
       root: "/path/to/g16"
       profile: "/path/to/g16/g16.profile"
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

#### 方式 1：命令行直接运行（推荐新手）

```bash
# 使用 SMILES 字符串运行单个反应
bin/rph_run --smiles "C=C(C)C(=O)O" --output ./Output/rx_manual

# 指定反应类型使用正向扫描 (NEW)
bin/rph_run --smiles "C=C(C)C(=O)O" --reaction-type "[4+3]_default" --output ./Output/rx_4p3
```

支持的反应类型：
- `[5+2]_default` - 5+2 环加成 (retro_scan，默认)
- `[4+3]_default` - 4+3 环加成 (forward_scan)
- `[4+2]_default` - Diels-Alder 反应 (forward_scan)
- `[3+2]_default` - 1,3-偶极环加成 (forward_scan)

`rx_id` 来自配置（`run.single.rx_id`）或数据集的 ID 列。

#### 方式 2：配置文件驱动（推荐批量运行）

1. 编辑 `config/defaults.yaml` 中的 `run` 部分：
   ```yaml
   run:
     source: single
     single:
       product_smiles: "C=C(C)C(=O)O"
       rx_id: "rx_001"
   ```

2. 运行流水线：
   ```bash
   bin/rph_run
   ```

#### 方式 3：Python API

```python
from pathlib import Path
from rph_core.orchestrator import ReactionProfileHunter

# 初始化
hunter = ReactionProfileHunter(config_path="config/defaults.yaml")

# 运行单个反应
result = hunter.run_pipeline(
    product_smiles="C=C(C)C(=O)O",
    work_dir=Path("./Output/rx_001")
)

# 检查结果
if result.success:
    print(f"✅ 成功！特征文件: {result.features_csv}")
else:
    print(f"❌ 失败：{result.error_message}")
```

---

## 📖 详细文档

### 架构概述

ReactionProfileHunter 采用**串行四步流水线**设计：

```
S1: 锚定/构象搜索
    ↓ (product_min.xyz, precursor_min.xyz)
S2: 逆向扫描
    ↓ (ts_guess.xyz, reactant_complex.xyz)
S3: 过渡态优化/救援
    ↓ (ts_final.xyz, reactant_opt/, NBO artifacts)
S4: 特征提取与打包
    ↓ (features_raw.csv, features_mlr.csv, feature_meta.json)
```

#### 各步骤详解

| 步骤 | 功能 | 核心模块 | 关键输出 |
|------|------|----------|----------|
| **S1** | 从 SMILES 生成 3D 结构<br>构象搜索<br>DFT 优化 | `steps/anchor/`<br>`steps/conformer_search/` | `product_min.xyz`<br>`precursor_min.xyz` |
| **S2** | 逆向/正向扫描<br>生成 TS 初猜<br>识别反应物 | `steps/step2_retro/`<br>(retro_scan 或 forward_scan) | `ts_guess.xyz`<br>`reactant_complex.xyz` |
| **S3** | TS 优化与频率分析<br>反应物优化<br>救援策略 | `steps/step3_opt/` | `ts_final.xyz`<br>`reactant_opt/` |
| **S4** | 能量提取<br>几何特征计算<br>NBO/FMO 分析 | `steps/step4_features/` | `features_raw.csv`<br>`feature_meta.json` |

### 输入/输出规范

#### 必需输入

- **产物 SMILES 字符串**：例如 `"C=C(C)C(=O)O"`
- **反应 ID**（可选）：用于标识反应，自动生成或手动指定

#### 严格输出契约

RPH 遵循**严格的步骤输出契约**，每个步骤必须生成指定的文件：

**S1 输出（必需）**
- `S1_ConfGeneration/product/product_min.xyz`
- `S1_ConfGeneration/precursor/precursor_min.xyz`（可选）

**S2 输出（必需，两者缺一不可）**
- `S2_Retro/ts_guess.xyz` - TS 初猜结构
- `S2_Retro/reactant_complex.xyz` - 反应物复合物

**S3 输出（必需）**
- `S3_TS/ts_final.xyz` - 优化后的 TS 结构
- `S3_TS/reactant_sp.xyz` - 标准反应物几何
- `S3_TS/reactant_opt/standard/` 或 `S3_TS/reactant_opt/rescue/` - 反应物 OPT+Freq 运行目录

**S4 输出（必需）**
- `S4_Data/features_raw.csv` - 原始特征（包含所有计算值）
- `S4_Data/features_mlr.csv` - 机器学习特征（标准化）
- `S4_Data/feature_meta.json` - 特征元数据（版本、时间戳等）

**S4 可选输出（v5.4 NBO）**
- `S4_Data/qc_nbo.37` - NBO 分析文件（如果在 S3 中找到）

#### 输出目录示例

```
Output/rx_001/
├── S1_ConfGeneration/
│   ├── product/
│   │   └── product_min.xyz
│   └── precursor/
│       └── precursor_min.xyz
├── S2_Retro/
│   ├── ts_guess.xyz
│   └── reactant_complex.xyz
├── S3_TS/
│   ├── ts_final.xyz
│   ├── reactant_sp.xyz
│   └── reactant_opt/
│       └── standard/        # 或 rescue/
│           ├── input.gjf
│           ├── output.log
│           └── *.fchk
├── S4_Data/
│   ├── features_raw.csv
│   ├── features_mlr.csv
│   ├── feature_meta.json
│   └── qc_nbo.37           # 可选
└── rph.log
```

### 配置指南

#### 核心配置项

**1. 量子化学软件路径** (`executables`)
```yaml
executables:
  gaussian:
    path: "/root/g16/g16"
    use_wrapper: true
    wrapper_path: "./scripts/run_g16_worker.sh"
```

**2. 计算资源** (`resources`)
```yaml
resources:
  mem: "64GB"       # 总内存
  nproc: 16         # 并行核数
  orca_maxcore_safety: 0.8  # ORCA 内存安全系数
```

**3. 理论水平** (`theory`)
```yaml
  # xTB 预优化（推荐开启）
  preoptimization:
    enabled: true
    gfn_level: 2
    overlap_threshold: 1.0  # A

  # 几何优化
  optimization:
    method: B3LYP
    basis: def2-SVP
    dispersion: GD3BJ
    engine: gaussian

  # 高精度单点能
  single_point:
    method: wB97X-D3BJ
    basis: def2-TZVPP
    engine: orca
```

**4. 反应配置 profiles** (`reaction_profiles`) - 新增
```yaml
reaction_profiles:
  "[4+3]_default":
    forming_bond_count: 2
    s2_strategy: forward_scan       # 使用 xTB $scan
    scan:
      scan_start_distance: 1.8      # 初始键长 (Å)
      scan_end_distance: 3.2       # 终止键长 (Å)
      scan_steps: 20                # 扫描步数
      scan_mode: concerted          # 或 "sequential"
      scan_force_constant: 0.5      # 约束力常数

  "[4+2]_default":
    forming_bond_count: 2
    s2_strategy: forward_scan
    scan:
      scan_start_distance: 2.0
      scan_end_distance: 3.5
      scan_steps: 15

  "[5+2]_default":
    forming_bond_count: 2
    s2_strategy: retro_scan         # 传统逆向扫描

  "_universal":
    forming_bond_count: 2
    s2_strategy: forward_scan
    scan:
      scan_start_distance: 2.2
      scan_end_distance: 3.5
```

**4. 优化控制** (`optimization_control`)
```yaml
optimization_control:
  max_cycles: 100
  convergence:
    level: normal  # loose, normal, tight, verytight
  hessian:
    initial: calcfc
    recalc_every: 10
```

#### 配置文件位置

- **主配置**：`config/defaults.yaml` - 所有默认参数
- **模板**：`config/templates/` - Gaussian/ORCA 输入模板
- **环境变量**：可通过环境变量覆盖路径配置

### 高级用法

#### 批量运行

创建 TSV 文件 `reactions.tsv`：
```tsv
rx_id	product_smiles	precursor_smiles
rx_001	C=C(C)C(=O)O	C=C(C)C(=O)OC
rx_002	CC(=O)OC	CC(=O)Cl
```

在 `config/defaults.yaml` 中配置：
```yaml
run:
  source: batch
  batch:
    input_file: "reactions.tsv"
    smiles_column: "product_smiles"
```

运行：
```bash
bin/rph_run
```

#### 数据集模式

```yaml
run:
  source: dataset
  dataset:
    format: "csv"
    path: "data/reaxys_cleaned.csv"
    product_smiles_col: "product_smiles_main"
    id_col: "rx_id"
    delimiter: ","
```

#### 断点恢复

如果计算意外中断，重新运行相同命令即可继续：
```bash
bin/rph_run --smiles "C=C(C)C(=O)O" --output ./Output/rx_001
```

通过配置控制恢复策略：
```yaml
run:
  resume: true
  resume_rehydrate: true
  resume_rehydrate_policy: best_effort
```

#### 自定义参数

创建自定义配置并通过 `--config` 使用：
```bash
bin/rph_run --config config/custom.yaml
```

示例覆盖：
```yaml
theory:
  optimization:
    method: M06-2X
resources:
  nproc: 32
```

#### 启用 NBO 分析

**方式 1：配置文件**
```yaml
step3:
  reactant_opt:
    enable_nbo: true  # 在反应物 OPT+Freq 中启用 NBO
```

**NBO 文件收集规则**：
- 自动搜索 `S3_TS` 的以下子目录：
  - `nbo_analysis/`
  - `nbo/`
  - `reactant_opt/standard/`
  - `reactant_opt/rescue/`
- 识别文件扩展名：`*.37`, `*.nbo`, `*.nbo7`
- 复制到 `S4_Data/qc_nbo.37`（统一命名）

---

## 🧪 测试

### 运行测试

```bash
# 运行全部测试
pytest -v tests/

# 运行单个测试文件
pytest tests/test_s4_no_qc_execution.py -v

# 运行单个测试函数
pytest tests/test_s4_no_qc_execution.py::test_extractor_degrades_gracefully -v

# 快速 CI 检查（导入测试 + 无 QC 测试）
pytest tests/test_imports_step4_features.py tests/test_s4_no_qc_execution.py -v

# 仅运行 S4 合同测试
pytest tests/test_s4_*.py tests/test_m2_*.py tests/test_m4_*.py -v

# 测试覆盖率
pytest --cov=rph_core --cov-report=html
```

### 导入风格检查（CI 门控）

```bash
# 检查多点多级相对导入（违规时退出码 1）
python scripts/ci/check_imports.py rph_core
```

### 注意事项

- `tests/conftest.py` 自动将仓库根目录添加到 `sys.path`，无需可编辑安装即可运行测试
- 集成测试使用 mock 量子化学计算，不需要真实的 Gaussian/ORCA
- 详见 [tests/AGENTS.md](tests/AGENTS.md) 了解测试组织方式

---

## 🐛 故障排查

### 常见问题

#### 1. Gaussian 找不到可执行文件

**错误信息**：
```
FileNotFoundError: Gaussian executable not found: /root/g16/g16
```

**解决方案**：
- 检查 `config/defaults.yaml` 中的 `executables.gaussian.path` 是否正确
- 确认已安装 Gaussian 并设置环境变量
- 尝试直接运行 `g16 < test.gjf` 测试 Gaussian 是否可用

#### 2. 内存不足导致计算崩溃

**错误信息**：
```
Error termination via Lnk1e in /root/g16/g16
```

**解决方案**：
- 降低 `resources.mem` 和 `resources.nproc`
- 使用更小的基组（如 `def2-SVP` → `6-31G*`）
- 启用 ORCA 的 `maxcore` 限制

#### 3. 构象搜索找不到低能构象

**错误信息**：
```
ConformerSearchError: No conformers found below energy threshold
```

**解决方案**：
- 增加 CREST 的采样时间
- 调整能量窗口阈值
- 检查输入 SMILES 是否正确

#### 4. TS 优化失败

**错误信息**：
```
TSOptimizationError: TS optimization failed after 5 attempts
```

**解决方案**：
- 检查 `ts_guess.xyz` 是否合理（虚频数量、几何结构）
- 启用救援策略（自动）
- 手动调整 TS 初猜结构

#### 5. 特征提取缺失 NBO 数据

**错误信息**：
```
Warning: NBO file not found in S3 subdirectories
```

**解决方案**：
- 确认已启用 `step3.reactant_opt.enable_nbo: true`
- 检查 `S3_TS/reactant_opt/standard/` 中是否有 `*.37` 或 `*.nbo` 文件
- 查看 Gaussian/ORCA 输出日志，确认 NBO 计算成功

### 日志文件

所有运行日志保存在 `Output/[rx_id]/rph.log`，包含：
- 每个步骤的详细执行信息
- 量子化学计算的标准输出/错误
- 异常堆栈跟踪

查看日志：
```bash
tail -f Output/rx_001/rph.log
```

### 调试模式

启用详细日志：
```bash
bin/rph_run --smiles "..." --log-level DEBUG
```

### 获取帮助

如果问题仍未解决：
1. 查看 [GitHub Issues](https://github.com/yourusername/ReactionProfileHunter/issues)
2. 附上完整错误日志和配置文件
3. 提供最小可复现示例

---

## 📚 进阶阅读

### 文档索引

| 文档 | 内容 |
|------|------|
| [`AGENTS.md`](AGENTS.md) | **智能化编程指南** —— 构建/测试命令、代码规范、约定 |
| [`rph_core/steps/AGENTS.md`](rph_core/steps/AGENTS.md) | 各步骤模块架构说明 |
| [`config/AGENTS.md`](config/AGENTS.md) | 配置文件结构说明 |
| [`tests/AGENTS.md`](tests/AGENTS.md) | 测试组织与约定 |
| [`scripts/AGENTS.md`](scripts/AGENTS.md) | 脚本与 CI 工具说明 |
| [`docs/DUAL_LEVEL_STRATEGY_SUMMARY.md`](docs/DUAL_LEVEL_STRATEGY_SUMMARY.md) | 双层级计算策略详解 |
| [`docs/QUICK_START_DUAL_LEVEL.md`](docs/QUICK_START_DUAL_LEVEL.md) | 双层级计算快速上手 |
| [`docs/BUGFIX_STEP3_QCRESULT.md`](docs/BUGFIX_STEP3_QCRESULT.md) | S3 QCResult 重构说明 |
| [`docs/S4_FEATURES_SUMMARY.md`](docs/S4_FEATURES_SUMMARY.md) | S4 特征提取功能总结 |
| [`docs/TESTING_GUIDE.md`](docs/TESTING_GUIDE.md) | 测试指南与环境准备 |

### 学术引用

如果您在研究中使用 ReactionProfileHunter，请引用：

```bibtex
@software{reactionprofilehunter2024,
  title = {ReactionProfileHunter: Automated Reaction Mechanism Exploration},
  author = {Your Name},
  year = {2024},
  url = {https://github.com/yourusername/ReactionProfileHunter},
  version = {6.2.0}
}
```

### 相关资源

- **理论背景**：
  - Houk 课题组的计算化学方法论（推荐阅读）
  - Grimme 的 GFN-xTB 方法论文
  - wB97X-D3 泛函性能评估

- **工具文档**：
  - [Gaussian 16 用户手册](https://gaussian.com/man/)
  - [ORCA 输入库](https://sites.google.com/site/orcainputlibrary/)
  - [xTB 文档](https://xtb-docs.readthedocs.io/)

---

## 🤝 贡献指南

我们欢迎各种形式的贡献！

### 如何贡献

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
- **仅使用绝对导入** —— 禁止多点多级相对导入（例如使用 `from rph_core.utils...` 而非 `from ...utils`）
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

## 📄 许可证

本项目采用 [MIT License](LICENSE) 开源协议。

---

## 📧 联系方式

- **作者**：QCcalc Team
- **邮箱**：your.email@example.com
- **项目主页**：https://github.com/yourusername/ReactionProfileHunter
- **问题反馈**：https://github.com/yourusername/ReactionProfileHunter/issues

---

## 🙏 致谢

- **Houk 教授**：提供双层级计算策略建议
- **Grimme 课题组**：开发 xTB 和 CREST 工具
- **Gaussian 和 ORCA 开发团队**：提供优秀的量子化学软件
- **社区贡献者**：感谢所有提交 issue 和 PR 的朋友

---

<div align="center">

**如果觉得本项目有帮助，请给我们一个 ⭐ Star！**

</div>
