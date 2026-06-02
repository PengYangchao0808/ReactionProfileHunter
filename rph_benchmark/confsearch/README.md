# S1 Conformer Search Benchmark

## 概述

本模块用于比较 RPH Step1 四套协议 (`ext` / `full` / `lite` / `zero`) 的构象搜索质量与计算成本。

**运行边界**：仅执行 S1，跳过 S0/S2/S3/S4。

## 脚本说明

| 脚本 | 定位 | 用法 |
|------|------|------|
| `run.sh` | 正式 benchmark 编排器 | 多反应 x 多协议矩阵运行，自动调用 evaluate.py |
| `evaluate.py` | 后评估（纯文件读取） | 解析 provenance + conformer energies，输出 cross-protocol 报告 |
| `smoke_test.sh` | 单反应快速冒烟测试 | 通过 `sed` 切换 defaults.yaml 中的协议，逐轮运行验证 |

## 快速开始

```bash
cd ReactionProfileHunter

# 跑指定反应
bash benchmark/confsearch/run.sh --rx-id 1,3,8

# 跑指定反应和协议
bash benchmark/confsearch/run.sh --rx-id 1 --protocols ext,lite

# 跑全部反应
bash benchmark/confsearch/run.sh --all

# 单反应冒烟测试
bash benchmark/confsearch/smoke_test.sh
bash benchmark/confsearch/smoke_test.sh lite zero
```

## run.sh 参数

```
--rx-id ID           从 CSV 选择反应 (可逗号分隔或重复指定)
--all                运行 CSV 中所有反应
--protocols LIST     协议列表 (默认: ext,full,lite,zero)
--reaction-type TYPE 反应类型 (默认: [4+3]_default)
--output DIR         输出根目录 (默认: PROJECT_ROOT/Output/benchmark)
--force-clean        运行前删除已有输出 (默认)
--no-clean           保留已有输出 (断点续算模式)
```

## 输出结构

```
Output/benchmark/
└── rx{id}/{protocol}/
    ├── RXN_*/S1_ConfGeneration/
    │   ├── product/product_min.xyz
    │   ├── product/product_global_min.xyz
    │   ├── precursor/precursor_global_min.xyz
    │   └── provenance.json
    └── evaluation/
        ├── summary.json
        └── summary.md
```

## evaluate.py 输出

`evaluate.py` 在每个反应的 `evaluation/` 子目录下生成：

- `summary.json` — 机器可读的完整指标
- `summary.md` — 人类可读的 markdown 报告，包含：
  - Per-Protocol Summary (状态、funnel 统计、能量)
  - Per-Molecule Detail (每个分子的 DFT count、SP 能量、Gibbs 能量)
  - Cross-Protocol Comparison (与 ext 的 delta E、delta G、DFT count ratio)
  - Protocol Behavior Summary (handoff mode、fallback、selection mode)

## 已知限制

项目路径可能包含空格和方括号（如 `[4+3]`），导致 ORCA/Gaussian 走 sandbox 路径（先复制到 `/tmp` 再拷回），带来 I/O 开销。可通过 `--output` 指定到安全路径缓解：

```bash
bash benchmark/confsearch/run.sh --rx-id 1 --output /tmp/rph_benchmark
```

## 依赖

- Python 3.8+（标准库即可，无额外 pip 包）
- RPH 管线配置的 QC 后端（ORCA/Gaussian + xTB + CREST）
- `data/reaxys_cleaned.csv` 数据集
