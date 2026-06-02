# 从 rph_output_backup 提取特征的完整方案

> 日期: 2026-05-23  
> 数据源: `external_data/rph_output_backup/` (12 个 RXN 目录)

---

## 1. 数据现状

### 1.1 目录结构

```
external_data/rph_output_backup/
├── RXN_0b71b9e9/
│   ├── reaction_manifest.json   ← RXN 级元数据 (version, condition_ids, theory_signature)
│   ├── S0_Mechanism/            ← S0: 机理分类结果
│   ├── conditions/              ← 条件级热力学 + 特征合并
│   │   ├── COND_8/
│   │   ├── COND_9/
│   │   └── COND_10/
│   │       ├── condition_manifest.json  ← 条件元数据 (temperature, solvent, yield, dr)
│   │       ├── merged_features.csv       ← 该条件下的合并特征
│   │       ├── condition_features_mlr.csv ← ML-ready 特征
│   │       ├── thermo_calculation.json    ← 热力学计算结果
│   │       └── branches/
│   │           ├── BR_MAJOR/  ← 主路径热力学
│   │           └── BR_DR_001/ ← DR 路径热力学
│   ├── branches/
│   │   ├── BR_MAJOR/
│   │   │   ├── branch_manifest.json
│   │   │   ├── S0_Mechanism/   ← 该路径的 S0
│   │   │   ├── S1_ConfGeneration/ ← S1: 构象搜索结果 (xyz, crest, cluster)
│   │   │   ├── S2_Retro/       ← S2: 逆合成扫描 (ts_guess, intermediate)
│   │   │   ├── S3_TS/          ← S3: TS 优化 (ts_final, sp_matrix)
│   │   │   │   ├── ts_final.xyz
│   │   │   │   ├── mechanism_meta.json
│   │   │   │   ├── artifacts_index.json
│   │   │   │   ├── sp_matrix_metadata.json
│   │   │   │   └── ts_opt/berny/  ← Gaussian 优化产物
│   │   │   └── S4_Data/        ← ★ 特征已提取
│   │   │       ├── features_raw.csv    (62 列)
│   │   │       ├── features_mlr.csv    (22 列)
│   │   │       ├── feature_meta.json
│   │   │       └── deployable_features.csv
│   │   └── BR_DR_001/          ← DR 路径
│   │       ├── branch_manifest.json
│   │       ├── S0_Mechanism/
│   │       ├── S1_ConfGeneration/
│   │       ├── S2_Retro/
│   │       ├── S3_TS/          ← S3 部分完成 (含 ts_final.xyz)
│   │       └── ★ 缺少 S4_Data ← 需要提取特征
│   └── pipeline.state
└── small_molecules/             ← S1 小分子缓存
```

### 1.2 完成度统计

| 分支 | S3 完成 | S4 完成 | 缺失 |
|------|---------|---------|------|
| BR_MAJOR | 12/12 ✅ | 12/12 ✅ | 无 |
| BR_DR_001 | 9/12 ✅ S3 | 0/12 ❌ | 全部缺少 S4, 3 个缺少 S3 |

### 1.3 已有特征数据质量

`features_raw.csv` (BR_MAJOR) 包含 60+ 特征列:
- `thermo.*`: dG_activation, dG_reaction, dE_activation, method, solvent
- `geom.*`: natoms_ts, r1, r2, asynch, rg_ts, close_contacts, dr
- `qc.*`: has_gibbs, sp_report_validated, forming_bonds_valid, sample_weight
- `s1_*`: E_avg_weighted, E_span, Nconf_eff, Sconf, dG_act, tau
- `s2_*`: eps_homo, eps_lumo, eta, gedt_value, mu, omega, q_fragment
- `ts.*`: dipole_debye, imag1_cm1_abs, n_imag

`condition_manifest.json` 包含:
- temperature, solvent, catalyst, yield_pct, dr_major, dr_minor, ee_pct

---

## 2. 总体策略

分两阶段执行:

```
Phase A: 补全 BR_DR_001 的 S4 缺漏
  └─ 对 9 个有 S3 的 DR 分支运行 rph-features extract
  └─ 跳过 3 个无 S3 的分支

Phase B: 整合为 ML 数据集
  └─ 合并 BR_MAJOR + BR_DR_001 的 feature CSV
  └─ 用 condition_manifest 注入 yield/dr/temperature 标签
  └─ 输出统一 ML 数据集到 rph_features_output/
```

---

## 3. Phase A: 补全 S4 特征

### 3.1 问题

`rph-features extract` 期望的输入是单层 S0–S3 目录结构:

```
预期                       实际 (分支结构)
────                        ────────
<work_dir>/                  branches/BR_MAJOR/
├── S1_ConfGeneration/       ├── S1_ConfGeneration/
├── S2_Retro/                ├── S2_Retro/
└── S3_TS/                   └── S3_TS/
```

需要编写一个适配器脚本，对每个分支创建符号链接或直接传递分支路径。

### 3.2 适配方案 (批处理脚本)

```bash
scripts/extract_backup_features.sh
```

核心逻辑:

```python
# scripts/backup_feature_extract.py
"""
对 external_data/rph_output_backup/ 中缺失 S4 的分支执行特征提取。
用法:
    python scripts/backup_feature_extract.py
    python scripts/backup_feature_extract.py --single RXN_0b71b9e9
"""

import subprocess, sys, json
from pathlib import Path

BACKUP = Path("external_data/rph_output_backup")
FEATURES_OUT = Path("rph_features_output")
RPH_FEATURES = sys.executable + " -m rph_features.cli"

def main():
    rxns = [sys.argv[2]] if len(sys.argv) > 2 and sys.argv[1] == "--single" \
           else sorted(BACKUP.glob("RXN_*"))

    for rxn_dir in rxns:
        rxn_id = rxn_dir.name
        for branch_dir in sorted((rxn_dir / "branches").iterdir()):
            branch_id = branch_dir.name
            s4_dir = branch_dir / "S4_Data"
            s3_exist = (branch_dir / "S3_TS" / "ts_final.xyz").exists()

            if s4_dir.exists():
                print(f"✓ {rxn_id}/{branch_id}: S4 already exists, skipping")
                continue
            if not s3_exist:
                print(f"✗ {rxn_id}/{branch_id}: No S3 ts_final.xyz, cannot extract")
                continue

            # 输出路径
            out_dir = FEATURES_OUT / rxn_id / branch_id
            out_dir.mkdir(parents=True, exist_ok=True)

            # 调用 rph-features extract
            cmd = [
                sys.executable, "-m", "rph_features.cli", "extract",
                "--rph-run", str(branch_dir.resolve()),
                "--output", str(out_dir.resolve()),
            ]
            print(f"→ {rxn_id}/{branch_id}: extracting features...")
            result = subprocess.run(cmd, capture_output=True, text=True)
            print(result.stdout)
            if result.returncode != 0:
                print(f"  ERROR: {result.stderr[:500]}")
            else:
                print(f"  ✓ Done → {out_dir}")

    # Phase B: 整合数据集
    consolidate_dataset(FEATURES_OUT)

def consolidate_dataset(features_root: Path):
    """整合所有分支和条件特征为统一 ML 数据集"""
    import pandas as pd

    all_features = []
    for rxn_dir in features_root.iterdir():
        if not rxn_dir.is_dir():
            continue
        for branch_dir in rxn_dir.iterdir():
            features_csv = branch_dir / "features" / "features_mlr.csv"
            if features_csv.exists():
                df = pd.read_csv(features_csv)
                df["reaction_id"] = rxn_dir.name
                df["branch_id"] = branch_dir.name
                all_features.append(df)

    # 合并 BR_MAJOR 的原有 S4 数据
    for rxn_dir in BACKUP.glob("RXN_*"):
        orig_s4 = rxn_dir / "branches" / "BR_MAJOR" / "S4_Data" / "features_mlr.csv"
        if orig_s4.exists():
            df = pd.read_csv(orig_s4)
            df["reaction_id"] = rxn_dir.name
            df["branch_id"] = "BR_MAJOR"
            all_features.append(df)

    if not all_features:
        print("No features found to consolidate")
        return

    dataset = pd.concat(all_features, ignore_index=True)

    # 注入条件标签 (从 condition_manifest.json)
    labels = []
    for rxn_dir in BACKUP.glob("RXN_*"):
        for cond_dir in (rxn_dir / "conditions").iterdir():
            manifest_file = cond_dir / "condition_manifest.json"
            if manifest_file.exists():
                m = json.loads(manifest_file.read_text())
                labels.append(m)
    if labels:
        label_df = pd.DataFrame(labels)
        dataset = dataset.merge(label_df, on="reaction_id", how="left")

    out_path = features_root / "unified_ml_dataset.csv"
    dataset.to_csv(out_path, index=False)
    print(f"✅ Unified dataset: {out_path} ({len(dataset)} samples, {len(dataset.columns)} cols)")

if __name__ == "__main__":
    main()
```

### 3.3 先行验证 (先试一个 RXN)

```bash
python scripts/backup_feature_extract.py --single RXN_0b71b9e9
```

验证输出：
```
rph_features_output/
└── RXN_0b71b9e9/
    └── BR_DR_001/
        ├── features/features_raw.csv
        ├── features/features_mlr.csv
        ├── features/feature_meta.json
        └── provenance/manifest.json
```

### 3.4 全面运行

```bash
# 后台运行，输出到日志
nohup python scripts/backup_feature_extract.py > logs/backup_extract.log 2>&1 &
```

---

## 4. Phase B: 整合为 ML 数据集

### 4.1 数据源

| 来源 | 路径 | 内容 |
|------|------|------|
| BR_MAJOR 原有 S4 | `branches/BR_MAJOR/S4_Data/features_mlr.csv` | 12 条主路径特征 (已有) |
| BR_DR_001 新提取 | `rph_features_output/RXN_*/BR_DR_001/features/features_mlr.csv` | ~9 条 DR 路径特征 (新提取) |
| Condition 标签 | `conditions/COND_*/condition_manifest.json` | yield, dr, temperature, solvent |

### 4.2 合并逻辑

```
rph_features_output/
└── unified_ml_dataset.csv    ← 最终 ML 数据集
```

Schema:

| 列 | 来源 | 说明 |
|----|------|------|
| reaction_id | RXN 目录名 | 关联标识 |
| condition_id | condition_manifest | 条件 ID |
| branch_id | BR_MAJOR/BR_DR_001 | 分支标识 |
| thermo.dE_activation | features_mlr.csv | 活化能 |
| geom.r_avg | features_mlr.csv | 平均键距 |
| ... 其他特征 | features_mlr.csv | 全部 MLR 特征 |
| yield_pct | condition_manifest | 产率 (目标变量) |
| dr_major | condition_manifest | 非对映选择性 (目标变量) |
| temperature_K | condition_manifest | 温度 (特征) |
| solvent | condition_manifest | 溶剂 (分类特征) |

### 4.3 使用 rph_ml 继续

```bash
# 诊断特征质量
rph-features diagnose --features-root ./rph_features_output

# 选择特征
rph-ml select-features \
    --dataset ./rph_features_output/unified_ml_dataset.csv \
    --output ./rph_ml_output/selected/

# 训练产率预测模型
rph-ml train-yield \
    --dataset ./rph_features_output/unified_ml_dataset.csv \
    --output ./rph_ml_output/models/

# 可视化结果
rph-ml visualize \
    --results ./rph_ml_output/models/ \
    --output ./rph_ml_output/plots/
```

---

## 5. 实施步骤明细

### Step 1: 创建批处理脚本

```bash
cat > scripts/backup_feature_extract.py << 'PYEOF'
# 上面 3.2 节的 Python 脚本
PYEOF
chmod +x scripts/backup_feature_extract.py
```

### Step 2: 确保 rph_features 可导入

```bash
# 确保 rph_features 包在 Python path 中
# 方式 A: pip install (推荐)
pip install -e ./rph_features

# 方式 B: 临时添加到 PYTHONPATH
export PYTHONPATH="$PYTHONPATH:$(pwd)/rph_features"
```

### Step 3: 单例验证

```bash
python scripts/backup_feature_extract.py --single RXN_0b71b9e9
```

验证输出文件存在且 CSV 格式正确:
```bash
head -3 rph_features_output/RXN_0b71b9e9/BR_DR_001/features/features_raw.csv
python -c "import pandas as pd; df=pd.read_csv('rph_features_output/RXN_0b71b9e9/BR_DR_001/features/features_raw.csv'); print(df.shape, df.columns.tolist()[:5])"
```

### Step 4: 全面运行

```bash
mkdir -p logs
nohup python scripts/backup_feature_extract.py > logs/backup_extract.log 2>&1 &
```

监控进度:
```bash
tail -f logs/backup_extract.log
```

### Step 5: 整合数据集

数据集整合在批处理脚本的 `consolidate_dataset()` 函数中自动完成。
也可手动执行:

```bash
python scripts/backup_feature_extract.py --consolidate-only
```

### Step 6: rph_ml 建模

```bash
# 确保 rph_ml 可导入
pip install -e ./rph_ml

# 诊断
rph-ml-diagnose --features-root ./rph_features_output

# 训练
rph-ml-train-yield \
    --dataset ./rph_features_output/unified_ml_dataset.csv \
    --target yield_pct \
    --output ./rph_ml_output/models/

# 可视化
rph-ml-visualize \
    --results ./rph_ml_output/models/ \
    --output ./rph_ml_output/plots/
```

---

## 6. 可能遇到的问题及对策

| 问题 | 对策 |
|------|------|
| BR_DR_001 S3 只有 ts_final.xyz 缺少 .log/.fchk | `rph-features extract` 会自然降级 (DEGRADED 状态)，在 feature_meta 中标记缺失 |
| rph_features 还未完全独立，缺少 rph_core utils | 确保 `rph_core` 在 PYTHONPATH 中 (pip install -e . 或在项目根目录运行) |
| 分支目录包含空格/特殊字符 | path_compat.is_toxic_path 已处理此问题，但最好确认路径安全性 |
| 新提取的特征与 BR_MAJOR 原有特征列名不一致 | `consolidate_dataset()` 会自动对齐列名，缺失列用 NaN 填充 |
| 3 个 BR_DR_001 缺少 S3 ts_final.xyz | 跳过，在日志中记录为 incomplete |

---

## 7. 预期输出

```
rph_features_output/
├── RXN_0b71b9e9/
│   └── BR_DR_001/
│       ├── features/features_raw.csv
│       ├── features/features_mlr.csv
│       ├── features/feature_meta.json
│       └── provenance/manifest.json
├── ...
├── RXN_f5d5c7b9/
│   └── BR_DR_001/
│       └── features/...       (9 个 DR 分支，预计 6-9 个成功)
├── unified_ml_dataset.csv      ← 最终 ML 数据集
│   (约 21-24 样本: 12 BR_MAJOR + 9 BR_DR_001)
│
rph_ml_output/
├── models/                     ← 训练好的模型
├── plots/                      ← 可视化图表
└── selected/features.txt       ← 选中的特征列表
```

---

## 8. 快速启动

一句话执行 (全部自动化):

```bash
# 1. 确保包可导入
pip install -e ./rph_features
pip install -e ./rph_ml

# 2. 一键执行: 特征提取 + 整合 + 建模
python scripts/backup_feature_extract.py
rph-ml train-yield --dataset rph_features_output/unified_ml_dataset.csv --output ./rph_ml_output/models/
```
