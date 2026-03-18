# Outcome 脚本使用指南

## 概述

`Outcome-20251023.sh` 是一个用于计算化学数据处理的 Bash 脚本，主要功能是从量子化学计算输出文件中提取单点能量和热力学数据，并生成最终的 Gibbs 自由能结果。

## 主要功能

- **自动目录检测**：智能识别项目根目录和输入输出路径
- **单点能量提取**：从 SP 计算输出文件中提取 FINAL SINGLE POINT ENERGY
- **热力学计算**：使用 Shermo 软件计算热力学性质
- **数据归档**：支持旧数据迁移和归档
- **多格式输出**：生成 TXT、CSV 和 Excel 格式的最终结果

## 系统要求

- **Bash shell**：必须使用 Bash 运行此脚本
- **Shermo 软件**：需要在 PATH 中可用
- **Python 3**（可选）：用于生成 Excel 文件

## 主要流程

### 1. 目录结构检测

```bash
# 自动向上查找包含 OPT/SP 的父目录（最多 5 层）
PROJECT_ROOT="$(detect_project_root "$OUTCOME_DIR")"

# 输入目录检测（大小写不敏感）
SP_DIR="$PROJECT_ROOT/$SP_BASE/Finished"  # 或裸目录
OPT_DIR="$PROJECT_ROOT/$OPT_BASE/Finished" # 或裸目录

# 输出目录固定
SUM_DIR="$OUTCOME_DIR/SUM"
```

### 2. 预处理和旧数据归档

```bash
# 环境变量控制
INGEST_LEGACY=1      # 启用旧数据归档（默认）
FORCE_RECALC=0       # 强制重新计算（默认关闭）

# 归档内容：
# - 旧 AA 文件（项目名.txt）
# - 旧 SUM 目录内容
# - 旧 FINAL 结果文件
```

### 3. 单点能量提取（SP 阶段）

**处理文件**：`SP_DIR/*.out`（排除 SMD 和模板文件）

**关键步骤**：
- 过滤 SMD 溶剂模型文件（`*.smd_SP.out`）
- 提取 "FINAL SINGLE POINT ENERGY" 值
- 规范化基名（去除 `_SP`、`_OPT` 后缀）
- 输出到 `AA_FILE`（项目名.txt）

### 4. 热力学计算（Shermo 阶段）

**处理文件**：`OPT_DIR/*.out`（排除 SMD 和模板文件）

**关键步骤**：
- 程序类型检测（ORCA/Gaussian）
- 频率计算验证
- 创建规范化临时副本
- 调用 Shermo 注入 SP 能量
- 输出到 `SUM_DIR/*_SUM.out`

### 5. 最终结果生成

**输入**：`SUM_DIR/*_SUM.out`

**输出格式**：
- **TXT**：`项目名_FINAL.txt` - 可读格式
- **CSV**：`项目名_FINAL.csv` - 逗号分隔值
- **XLSX**：`项目名_FINAL.xlsx` - Excel 格式（需要 Python3）

**数据提取**：匹配 "thermal correction to G" 或 "Sum of electronic Free Energies"

## 环境变量配置

| 变量名 | 默认值 | 说明 |
|--------|--------|------|
| `INGEST_LEGACY` | 1 | 启用旧数据归档 |
| `INGEST_LINK` | 0 | 使用软链接而非复制 |
| `INGEST_MOVE` | 1 | 移动而非复制旧文件 |
| `FORCE_RECALC` | 0 | 强制重新计算所有阶段 |
| `STRICT_FREQ_CHECK` | 1 | 严格频率计算检查 |
| `PROJECT_ROOT` | - | 手动指定项目根目录 |

## 文件过滤规则

### 排除的文件类型：
- **SMD 溶剂模型**：任何包含 `.smd` 的文件名
- **模板文件**：`Temp.out`、`Temp_OPT.out`
- **无频率计算**：未检测到频率块的文件

### 基名规范化：
- `base_SP` → `base`
- `base_OPT` → `base`
- 其他保持不变

## 典型使用方式

```bash
# 基本运行
cd /path/to/Outcome
bash Outcome-20251023.sh

# 强制重新计算
FORCE_RECALC=1 bash Outcome-20251023.sh

# 指定项目根目录
PROJECT_ROOT=/custom/path bash Outcome-20251023.sh

# 禁用旧数据归档
INGEST_LEGACY=0 bash Outcome-20251023.sh
```

## 输出文件说明

### 中间文件：
- `项目名.txt` - 单点能量映射文件
- `SUM/*_SUM.out` - Shermo 计算结果

### 最终结果：
- `项目名_FINAL.txt` - 可读的 Gibbs 自由能列表
- `项目名_FINAL.csv` - 结构化数据（名称,数值）
- `项目名_FINAL.xlsx` - Excel 格式数据

## 错误处理

脚本使用 `set -euo pipefail` 确保严格错误处理：

- **目录不存在**：优雅跳过或退出
- **程序未安装**：检查 Shermo 可用性
- **文件解析失败**：记录警告并继续
- **数据不匹配**：跳过无法处理的文件

## 调试技巧

```bash
# 查看路径检测
echo "项目根目录: $PROJECT_ROOT"
echo "SP目录: $SP_DIR"
echo "OPT目录: $OPT_DIR"

# 检查环境变量
env | grep -E "(INGEST|FORCE|STRICT)"

# 验证文件过滤
find "$SP_DIR" -name "*.out" | grep -v smd
```

## 注意事项

1. **文件命名**：避免使用特殊字符，确保基名唯一
2. **计算程序**：支持 ORCA 和 Gaussian 输出格式
3. **编码**：自动处理 CRLF/LF 行尾问题
4. **性能**：大文件集建议使用 `FORCE_RECALC=0` 避免重复计算
5. **备份**：重要数据建议在运行前备份

这个脚本为计算化学工作流提供了完整的自动化解决方案，从原始输出文件到最终的热力学数据报告。
