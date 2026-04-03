# ReactionProfileHunter 真实QC计算验证方案

**版本**: v6.2.0  
**日期**: 2026-03-11  
**验证目标**: 使用真实数据集中的一个[4+3]环加成反应进行端到端QC计算验证

---

## 一、验证目标

### 1.1 验证范围
- **完整S1-S4流水线**: 从SMILES到特征提取
- **双水平计算**: xTB → B3LYP/def2-SVP → wB97X-D3BJ/def2-TZVPP
- **反应类型**: [4+3] 环加成 (分子内)
- **计算时间预估**: 2-6小时 (取决于硬件)

### 1.2 选择示例反应

从 `data/reaxys_cleaned.csv` 中选择 **rx_id=9422028**:

```csv
rx_id: 9422028
产物SMILES: O=C1C[C@H]2C=C[C@]3(C[C@H]4COC(=O)N4[C@H]13)O2
前体SMILES: C=C=CN1C(=O)OCC1Cc1ccco1
反应类型: [4+3] 环加成 (allenamide类型)
置信度: HIGH
溶剂: DCM (二氯甲烷)
温度: -45°C
```

**选择理由**:
- 单产物反应 (is_multi_product=False)
- 高质量映射 (mapping OK)
- 分子内反应 (INTRA_TYPE_I)
- 适中的分子大小 (~20原子)
- 文献有收率数据支持

---

## 二、环境准备

### 2.1 确认QC软件安装

```bash
# 检查已安装的软件
which orca xtb crest g16

# 预期输出:
# /opt/software/orca/orca
# /opt/software/xtb/bin/xtb
# /opt/software/crest/crest
# /opt/software/gaussian/g16/g16

# 测试版本
orca --version
xtb --version
crest --version
g16 < /dev/null 2>&1 | head -5
```

### 2.2 工作目录设置

```bash
# 进入项目根目录
cd /mnt/e/Calculations/AI4S_ML_Studys/[4+3]\ Mechain\ learning/ReactionProfileHunter/ReactionProfileHunter_20260121

# 创建输出根目录
mkdir -p ./rph_output

# 验证Python环境
python -c "import rph_core; print(f'RPH版本: {rph_core.__version__}')"
```

### 2.3 资源配置确认

当前系统配置:
- **内存**: 建议 ≥32GB (配置为64GB)
- **CPU核数**: 建议 ≥8核 (配置为16核)
- **磁盘空间**: 确保 ≥10GB 可用空间

---

## 三、配置修改

### 3.1 创建验证专用配置

创建 `config/validation.yaml`:

```yaml
# ReactionProfileHunter 验证配置
# 基于 defaults.yaml 修改，降低计算成本以加速验证

# ==================== 执行路径配置 ====================
executables:
  orca:
    path: "/opt/software/orca/orca"
    ld_library_path: "/opt/software/orca"
  gaussian:
    path: "/opt/software/gaussian/g16/g16"
    root: "/opt/software/gaussian/g16"
    profile: "/opt/software/gaussian/g16/g16.profile"
    use_wrapper: false  # 验证时关闭wrapper简化流程
  xtb:
    path: "/opt/software/xtb/bin/xtb"
  crest:
    path: "/opt/software/crest/crest"
  shermo:
    path: "Shermo"

# ==================== 资源配置 ====================
resources:
  mem: "32GB"
  nproc: 8
  orca_maxcore_safety: 0.8

# ==================== 理论水平配置 ====================
theory:
  # xTB预优化
  preoptimization:
    enabled: true
    gfn_level: 2
    solvent: acetone
    nproc: 4
    overlap_threshold: 1.0
    opt_level: crude  # 快速模式

  # 几何优化 - 降低级别以加速验证
  optimization:
    method: B3LYP
    basis: def2-SVP
    dispersion: GD3BJ
    engine: gaussian
    nproc: 8
    mem: 32GB
    solvent: acetone

  # 高精度单点能 - 验证时可先用较小基组
  single_point:
    method: wB97X-D3BJ
    basis: def2-SVP  # 验证时临时降低，生产用def2-TZVPP
    engine: orca
    nproc: 8
    maxcore: 1000
    solvent: acetone
    fallback_to_gaussian: true

# ==================== 优化控制 ====================
optimization_control:
  timeout:
    enabled: true
    default_seconds: 3600  # 1小时超时
  
  oscillation:
    window_size: 10
    energy_tolerance: 0.0001
    max_oscillation_count: 3
  
  hessian:
    initial: calcfc
    recalc_every: 10
  
  convergence:
    level: normal

# ==================== Step1 构象搜索 ====================
step1:
  conformer_search:
    two_stage_enabled: false  # 验证时关闭两阶段，加速
    # 单阶段GFN2快速搜索
    stage1_gfn0:
      enabled: false
    stage2_gfn2:
      enabled: true
      gfn_level: 2

# ==================== Step3 反应物优化 ====================
step3:
  reactant_opt:
    enable_nbo: false  # 验证时关闭NBO以加速
    
# ==================== 运行配置 ====================
run:
  source: single
  single:
    product_smiles: "O=C1C[C@H]2C=C[C@]3(C[C@H]4COC(=O)N4[C@H]13)O2"
    rx_id: "rx_9422028"
    precursor_smiles: "C=C=CN1C(=O)OCC1Cc1ccco1"
    leaving_group_key: null
  
  output_root: "./rph_output"
  resume: false  # 验证时从头开始
  resume_rehydrate: false

# ==================== 反应档案 ====================
reaction_profiles:
  "[4+3]_default":
    forming_bond_count: 2
    s2_strategy: forward_scan
    scan:
      scan_start_distance: 1.8
      scan_end_distance: 3.2
      scan_steps: 15  # 验证时减少步数
      scan_mode: concerted
      scan_force_constant: 0.5
```

---

## 四、执行步骤

### 4.1 方式一: 使用CLI (推荐)

```bash
# 进入项目目录
cd /path/to/ReactionProfileHunter

# 运行验证
bin/rph_run \
  --config config/validation.yaml \
  --smiles "O=C1C[C@H]2C=C[C@]3(C[C@H]4COC(=O)N4[C@H]13)O2" \
  --output ./rph_output/rx_9422028

# 实时监控日志
tail -f ./rph_output/rx_9422028/rph.log
```

### 4.2 方式二: 使用Python API

创建 `run_validation.py`:

```python
#!/usr/bin/env python3
"""ReactionProfileHunter 真实QC验证脚本"""

import sys
import logging
from pathlib import Path

# 添加项目根目录到路径
sys.path.insert(0, str(Path(__file__).parent))

from rph_core.orchestrator import ReactionProfileHunter
from rph_core.utils.log_manager import setup_logging

def main():
    """运行验证计算"""
    
    # 配置日志
    setup_logging(
        level=logging.INFO,
        log_file=Path("./rph_output/rx_9422028/rph.log")
    )
    
    logger = logging.getLogger(__name__)
    logger.info("=" * 60)
    logger.info("RPH 真实QC计算验证开始")
    logger.info("=" * 60)
    
    # 初始化
    config_path = Path("config/validation.yaml")
    hunter = ReactionProfileHunter(config_path=config_path)
    
    # 产物和前体
    product_smiles = "O=C1C[C@H]2C=C[C@]3(C[C@H]4COC(=O)N4[C@H]13)O2"
    precursor_smiles = "C=C=CN1C(=O)OCC1Cc1ccco1"
    
    # 运行流水线
    result = hunter.run_pipeline(
        product_smiles=product_smiles,
        work_dir=Path("./rph_output/rx_9422028"),
        precursor_smiles=precursor_smiles,
        leaving_group_key=None,
        skip_steps=[]
    )
    
    # 检查结果
    if result.success:
        logger.info("=" * 60)
        logger.info("✅ 验证成功完成!")
        logger.info(f"产物: {result.product_smiles}")
        logger.info(f"前体: {result.precursor_smiles}")
        logger.info(f"输出目录: {result.work_dir}")
        
        # 检查特征文件
        features_csv = result.work_dir / "S4_Data" / "features_mlr.csv"
        if features_csv.exists():
            logger.info(f"特征文件: {features_csv}")
            
        logger.info("=" * 60)
        return 0
    else:
        logger.error("=" * 60)
        logger.error("❌ 验证失败")
        logger.error(f"错误: {result.error_message}")
        logger.error("=" * 60)
        return 1

if __name__ == "__main__":
    sys.exit(main())
```

执行:
```bash
python run_validation.py
```

---

## 五、分步执行与监控

### 5.1 Step 1: 构象搜索 (预计10-30分钟)

```bash
# 仅运行S1
bin/rph_run --config config/validation.yaml --skip-steps step2 step3 step4

# 检查输出
ls -la ./rph_output/rx_9422028/S1_ConfGeneration/
# 预期:
# S1_ConfGeneration/product/product_min.xyz
# S1_ConfGeneration/precursor/precursor_min.xyz
```

**成功标志**:
- `product_min.xyz` 存在且合理 (~20原子)
- `precursor_min.xyz` 存在且合理
- 无原子重叠警告

### 5.2 Step 2: 逆向扫描 (预计20-60分钟)

```bash
# 仅运行S2
bin/rph_run --config config/validation.yaml --skip-steps step3 step4

# 检查输出
ls -la ./rph_output/rx_9422028/S2_Retro/
# 预期:
# S2_Retro/ts_guess.xyz
# S2_Retro/reactant_complex.xyz
```

**成功标志**:
- TS猜测结构有合理的成键距离
- 反应物复合物能量合理

### 5.3 Step 3: TS优化 (预计1-3小时)

```bash
# 运行S3
bin/rph_run --config config/validation.yaml --skip-steps step4

# 检查输出
ls -la ./rph_output/rx_9422028/S3_TS/
# 预期:
# S3_TS/ts_final.xyz
# S3_TS/reactant_opt/standard/*.log
```

**成功标志**:
- TS有且仅有1个虚频
- 反应物优化收敛
- 能量合理 (检查 `.sp` 文件)

### 5.4 Step 4: 特征提取 (预计5-10分钟)

```bash
# 运行S4
bin/rph_run --config config/validation.yaml

# 检查输出
ls -la ./rph_output/rx_9422028/S4_Data/
# 预期:
# S4_Data/features_raw.csv
# S4_Data/features_mlr.csv
# S4_Data/feature_meta.json
```

---

## 六、结果验证清单

### 6.1 文件完整性检查

```bash
#!/bin/bash
# validate_output.sh

WORK_DIR="./rph_output/rx_9422028"

echo "=== RPH 输出验证 ==="

# S1检查
if [ -f "$WORK_DIR/S1_ConfGeneration/product/product_min.xyz" ]; then
    echo "✅ S1 产物优化完成"
    wc -l "$WORK_DIR/S1_ConfGeneration/product/product_min.xyz"
else
    echo "❌ S1 产物缺失"
fi

# S2检查
if [ -f "$WORK_DIR/S2_Retro/ts_guess.xyz" ]; then
    echo "✅ S2 TS猜测完成"
else
    echo "❌ S2 TS猜测缺失"
fi

# S3检查
if [ -f "$WORK_DIR/S3_TS/ts_final.xyz" ]; then
    echo "✅ S3 TS优化完成"
    
    # 检查虚频
    TS_LOG=$(find "$WORK_DIR/S3_TS" -name "*TS*.log" | head -1)
    if [ -f "$TS_LOG" ]; then
        IMAG=$(grep -c "imaginary frequencies" "$TS_LOG" 2>/dev/null || echo "0")
        echo "   虚频数量: $IMAG (期望: 1)"
    fi
else
    echo "❌ S3 TS优化缺失"
fi

# S4检查
if [ -f "$WORK_DIR/S4_Data/features_mlr.csv" ]; then
    echo "✅ S4 特征提取完成"
    echo "   特征数量:"
    head -1 "$WORK_DIR/S4_Data/features_mlr.csv" | tr ',' '\n' | wc -l
else
    echo "❌ S4 特征缺失"
fi

echo "=== 验证完成 ==="
```

### 6.2 能量合理性检查

```python
# check_energies.py
import json
from pathlib import Path

def check_energies(work_dir: Path):
    """检查能量值合理性"""
    
    # 读取mech_index.json
    mech_file = work_dir / "S4_Data" / "mech_index.json"
    if not mech_file.exists():
        print("❌ mech_index.json 不存在")
        return
    
    with open(mech_file) as f:
        data = json.load(f)
    
    print("=== 能量检查 ===")
    
    # 检查关键能量
    energies = data.get("energies", {})
    
    # 产物能量
    e_prod = energies.get("product_l2_sp")
    if e_prod:
        print(f"✅ 产物能量 (L2): {e_prod:.6f} Hartree")
    else:
        print("❌ 产物能量缺失")
    
    # 反应物能量
    e_react = energies.get("reactant_l2_sp")
    if e_react:
        print(f"✅ 反应物能量 (L2): {e_react:.6f} Hartree")
    else:
        print("❌ 反应物能量缺失")
    
    # TS能量
    e_ts = energies.get("ts_l2_sp")
    if e_ts:
        print(f"✅ TS能量 (L2): {e_ts:.6f} Hartree")
        
        # 检查活化能
        if e_react:
            ea = (e_ts - e_react) * 627.509  # kcal/mol
            print(f"   活化能: {ea:.2f} kcal/mol")
            
            if 0 < ea < 100:
                print("   ✅ 活化能在合理范围")
            else:
                print("   ⚠️ 活化能异常")
    else:
        print("❌ TS能量缺失")
    
    # 反应热
    if e_prod and e_react:
        delta = (e_prod - e_react) * 627.509
        print(f"✅ 反应热: {delta:.2f} kcal/mol")

if __name__ == "__main__":
    check_energies(Path("./rph_output/rx_9422028"))
```

### 6.3 特征完整性检查

```bash
# 检查CSV特征
cat ./rph_output/rx_9422028/S4_Data/features_mlr.csv | head -5

# 检查JSON元数据
python -m json.tool ./rph_output/rx_9422028/S4_Data/feature_meta.json | head -50
```

---

## 七、时间预估与资源监控

### 7.1 各阶段时间预估 (8核/32GB)

| 步骤 | 预估时间 | 主要计算 | 监控命令 |
|------|----------|----------|----------|
| S1构象搜索 | 10-30 min | CREST GFN2 | `htop` 查看crest进程 |
| S1 DFT优化 | 20-60 min | B3LYP/def2-SVP | `tail -f *.log` |
| S2逆向扫描 | 20-60 min | xTB $scan | 检查xtb输出 |
| S3 TS优化 | 30-90 min | B3LYP TS search | `grep "E(RB3LYP)" *.log` |
| S3反应物优化 | 30-60 min | B3LYP opt+freq | 检查收敛状态 |
| S3 L2单点 | 20-40 min | wB97X-D3BJ | ORCA输出 |
| S4特征提取 | 5-10 min | 解析文件 | 无QC计算 |
| **总计** | **2-5小时** | — | — |

### 7.2 资源监控脚本

```bash
# monitor.sh - 在另一个终端运行
while true; do
    clear
    echo "=== RPH 资源监控 ==="
    echo "时间: $(date)"
    echo ""
    
    # CPU和内存
    echo "CPU使用:"
    ps aux | grep -E "(orca|g16|xtb|crest)" | grep -v grep || echo "无QC进程"
    echo ""
    
    # 磁盘使用
    echo "磁盘使用:"
    du -sh ./rph_output/rx_9422028/* 2>/dev/null || echo "目录尚未创建"
    echo ""
    
    # 当前步骤
    echo "当前步骤:"
    ls -td ./rph_output/rx_9422028/*/ 2>/dev/null | head -1
    echo ""
    
    # 最新日志
    echo "最新日志:"
    tail -5 ./rph_output/rx_9422028/rph.log 2>/dev/null || echo "暂无日志"
    
    sleep 10
done
```

---

## 八、故障排查

### 8.1 常见问题

| 问题 | 症状 | 解决方案 |
|------|------|----------|
| **原子重叠** | S1 Gaussian崩溃 | 启用xTB预优化 (已默认启用) |
| **TS不收敛** | S3反复失败 | 检查S2猜测质量，调整TS参数 |
| **内存不足** | ORCA/Gaussian killed | 减少nproc，增加maxcore |
| **虚频过多** | TS有>1虚频 | 调整初始猜测，使用QST2 rescue |
| **负活化能** | TS能量低于反应物 | 检查结构连接性，重新优化反应物 |

### 8.2 日志分析

```bash
# 检查错误
 grep -i "error\|fail\|kill" ./rph_output/rx_9422028/rph.log

# 检查警告
 grep -i "warn" ./rph_output/rx_9422028/rph.log

# 查看最后100行
 tail -100 ./rph_output/rx_9422028/rph.log
```

---

## 九、验证成功标准

### 9.1 必需输出文件

```
rph_output/rx_9422028/
├── S1_ConfGeneration/
│   ├── product/product_min.xyz
│   └── precursor/precursor_min.xyz
├── S2_Retro/
│   ├── ts_guess.xyz
│   └── reactant_complex.xyz
├── S3_TS/
│   ├── ts_final.xyz
│   ├── reactant_sp.xyz
│   └── reactant_opt/standard/
│       ├── reactant_complex.log (收敛)
│       └── reactant_complex.fchk
└── S4_Data/
    ├── features_raw.csv
    ├── features_mlr.csv
    └── feature_meta.json
```

### 9.2 合理性检查

| 检查项 | 期望结果 | 验证方法 |
|--------|----------|----------|
| TS虚频 | 恰好1个 | 检查日志 |
| 活化能 | 5-50 kcal/mol | 能量差计算 |
| 反应热 | -30 to +30 kcal/mol | 能量差计算 |
| 成键距离 | 1.5-2.5 Å (TS) | 结构可视化 |
| 特征数量 | >50个 | CSV列数 |

---

## 十、后续分析

### 10.1 结果可视化

```python
# visualize.py
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# 读取特征
df = pd.read_csv("./rph_output/rx_9422028/S4_Data/features_mlr.csv")
print("提取的特征:")
print(df.T)

# 能量剖面 (如果有多个反应)
# plt.plot([0, 1, 2], [e_react, e_ts, e_prod], 'o-')
# plt.xlabel('Reaction Coordinate')
# plt.ylabel('Energy (kcal/mol)')
```

### 10.2 对比文献

```bash
# 提取关键值
echo "活化能: $(python check_energies.py | grep "活化能")"
echo "反应热: $(python check_energies.py | grep "反应热")"

# 与文献值对比 (rx_9422028文献收率等数据在CSV中)
```

---

## 十一、清理与归档

```bash
# 验证完成后清理临时文件
find ./rph_output/rx_9422028 -name "*.tmp" -delete
find ./rph_output/rx_9422028 -name "*.rwf" -delete  # Gaussian大文件

# 压缩归档
tar -czvf rx_9422028_validation.tar.gz ./rph_output/rx_9422028/

# 或删除 (节省空间)
# rm -rf ./rph_output/rx_9422028
```

---

**准备就绪**: 执行 `bin/rph_run --config config/validation.yaml` 开始验证  
**预计时间**: 2-5小时  
**监控**: `tail -f ./rph_output/rx_9422028/rph.log`
