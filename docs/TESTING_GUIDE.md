# 测试准备与运行流程指南

**项目**: ReactionProfileHunter v6.2  
**最后更新**: 2025-03-11  
**适用范围**: 开发环境、CI/CD 环境

---

## 🎯 快速开始 (Quick Start)

如果你只想快速验证环境是否正常工作，执行：

```bash
# 进入项目目录
cd /mnt/e/Calculations/AI4S_ML_Studys/[4+3] Mechain learning/ReactionProfileHunter/ReactionProfileHunter_20260121

# 运行最快的基础测试（约 5-10 秒）
pytest tests/test_imports_step4_features.py -v
```

如果显示所有测试通过 ✅，说明环境准备就绪！

---

## 📋 准备工作清单

### 1. 环境要求

#### Python 版本
```bash
python --version  # 需要 3.8+
```

#### 必需的 Python 包

**核心依赖**（运行测试必需）：
```bash
pip install pytest pytest-cov
```

**项目依赖**（完整功能）：
```bash
# 方式 1：如果项目有 setup.py 或 pyproject.toml
pip install -e .

# 方式 2：手动安装核心依赖
pip install numpy scipy pandas pyyaml networkx rdkit

# 方式 3：使用 requirements.txt（如果存在）
pip install -r requirements.txt
```

#### 项目结构验证

确保 `tests/conftest.py` 存在且包含：
```python
import sys
from pathlib import Path

project_root = Path(__file__).resolve().parents[1]
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))
```

这允许测试直接导入 `rph_core` 而无需安装包。

---

### 2. 量子化学软件（可选）

**重要**: 大部分测试 **不需要** 真实量子化学软件！

测试使用 **mock** 技术模拟 Gaussian/ORCA/xTB 的行为。

#### 不需要安装的测试 ✅
- `tests/test_imports_*.py` - 仅检查导入
- `tests/test_s4_no_qc_execution.py` - S4 特征提取（无 QC）
- `tests/test_mock_*.py` - Mock 集成测试
- `tests/test_m2_*.py` - S4 合同测试
- `tests/test_m4_*.py` - M4 机制测试

#### 需要真实 QC 软件的测试（默认跳过）
- `tests/test_qc_interface_v52.py` - QC 接口（真实运行）
- `tests/test_orca_interface.py` - ORCA 接口
- `tests/test_sandbox_toxic_paths.py` - 沙盒测试（需要 xTB）

这些测试默认被 `@pytest.mark.skipif(True, ...)` 装饰器跳过。

#### 如果你想运行真实 QC 测试

1. 安装所需软件：
   - Gaussian 16
   - ORCA 5.0+
   - xTB 6.5+
   - CREST 2.12+

2. 编辑 `config/defaults.yaml`：
   ```yaml
   executables:
     gaussian:
       path: "/your/path/to/g16"
     orca:
       path: "/your/path/to/orca"
     xtb:
       path: "/your/path/to/xtb"
   ```

3. 修改测试文件中的 skip 条件：
   ```python
   # 从：
   @pytest.mark.skipif(True, reason="Requires real ORCA environment")
   
   # 改为：
   @pytest.mark.skipif(shutil.which("orca") is None, reason="Requires ORCA")
   ```

---

## 🚀 测试运行流程

### 阶段 1：CI Gate（快速检查，~10秒）

这是每次提交前必须通过的测试：

```bash
# 运行导入检查和基础功能测试
pytest tests/test_imports_step4_features.py tests/test_s4_no_qc_execution.py -v
```

**预期输出**:
```
tests/test_imports_step4_features.py::test_import_core PASSED
tests/test_imports_step4_features.py::test_import_extractors PASSED
tests/test_s4_no_qc_execution.py::test_extractor_degrades_gracefully PASSED
...
========================= X passed in X.XXs ==========================
```

---

### 阶段 2：导入风格检查（CI 强制）

```bash
# 检查是否使用了禁止的多级相对导入
python ci/check_imports.py rph_core
```

**预期输出**:
```
✅ No forbidden import patterns found
📁 Scanned 198 Python files
✅ PASSED: Import style check
```

**如果失败**:
```
❌ IMPORT VIOLATIONS FOUND
📄 rph_core/steps/step4_features/extractors/geometry.py
   Line 15: from ...utils.file_io
   ⚠️  from ...utils (resolves to rph_core.steps.utils)
   ✅ FIX: Use 'from rph_core.utils...' instead
```

**修复**: 将 `from ...utils.file_io import X` 改为 `from rph_core.utils.file_io import X`

---

### 阶段 3：模块测试（~1-2分钟）

按功能模块运行测试：

#### S4 特征提取测试
```bash
pytest tests/test_s4_*.py -v
```

#### S3 优化测试
```bash
pytest tests/test_s3_checkpoint.py tests/test_m3_*.py -v
```

#### S2 逆向扫描测试
```bash
pytest tests/test_retro_scanner_v52.py tests/test_step2_path_compat.py -v
```

#### M2/M4 机制测试
```bash
pytest tests/test_m2_*.py tests/test_m4_*.py -v
```

---

### 阶段 4：完整测试套件（~5-10分钟）

```bash
# 运行所有测试
pytest -v tests/
```

**注意**: 某些测试可能会被跳过（标记为 `s`），这是正常的，表示这些测试需要真实 QC 环境。

---

### 阶段 5：覆盖率报告（可选）

```bash
pytest --cov=rph_core --cov-report=html
```

生成后打开 `htmlcov/index.html` 查看详细覆盖率报告。

---

## 🔧 常见测试场景

### 场景 1：开发新功能时

```bash
# 1. 先确保基础检查通过
pytest tests/test_imports_step4_features.py -v
python ci/check_imports.py rph_core

# 2. 运行相关模块测试
pytest tests/test_s4_*.py -v  # 如果修改 S4

# 3. 运行完整测试
pytest -v tests/
```

### 场景 2：调试特定失败测试

```bash
# 运行单个测试函数
pytest tests/test_s4_no_qc_execution.py::test_extractor_degrades_gracefully -v

# 添加更多输出
pytest tests/test_s4_no_qc_execution.py::test_extractor_degrades_gracefully -vvs

# 在失败时进入 PDB 调试
pytest tests/test_s4_no_qc_execution.py::test_extractor_degrades_gracefully --pdb
```

### 场景 3：CI/CD 流水线

```yaml
# .github/workflows/test.yml 示例
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.12'
      
      - name: Install dependencies
        run: |
          pip install pytest pytest-cov
          pip install numpy scipy pandas pyyaml networkx rdkit
      
      - name: Run import smoke tests
        run: pytest tests/test_imports_step4_features.py tests/test_s4_no_qc_execution.py -v
      
      - name: Run import style check
        run: python ci/check_imports.py rph_core
      
      - name: Run all tests
        run: pytest -v tests/
```

### 场景 4：提交前预检查

创建 `.git/hooks/pre-commit`：

```bash
#!/bin/bash

echo "🧪 Running pre-commit checks..."

# 1. 导入风格检查
echo "📋 Checking import style..."
python ci/check_imports.py rph_core
if [ $? -ne 0 ]; then
    echo "❌ Import style check failed"
    exit 1
fi

# 2. 快速测试
echo "🚀 Running quick tests..."
pytest tests/test_imports_step4_features.py tests/test_s4_no_qc_execution.py -v
if [ $? -ne 0 ]; then
    echo "❌ Quick tests failed"
    exit 1
fi

echo "✅ All checks passed!"
exit 0
```

使钩子生效：
```bash
chmod +x .git/hooks/pre-commit
```

---

## 🐛 故障排除

### 问题 1：`ModuleNotFoundError: No module named 'rph_core'`

**原因**: Python 无法找到项目模块

**解决**:
```bash
# 方法 1：使用 pytest 运行（推荐）
pytest tests/test_imports_step4_features.py -v
# conftest.py 会自动添加项目根目录到 sys.path

# 方法 2：手动添加 PYTHONPATH
export PYTHONPATH=/mnt/e/Calculations/AI4S_ML_Studys/[4+3] Mechain learning/ReactionProfileHunter/ReactionProfileHunter_20260121:$PYTHONPATH
python -m pytest tests/test_imports_step4_features.py -v

# 方法 3：安装为可编辑包
pip install -e .
pytest tests/test_imports_step4_features.py -v
```

### 问题 2：导入风格检查失败

**错误**:
```
from ...utils.file_io import read_xyz
```

**修复**:
```python
# ❌ 错误 - 多级相对导入
from ...utils.file_io import read_xyz

# ✅ 正确 - 绝对导入
from rph_core.utils.file_io import read_xyz
```

### 问题 3：测试报告 `No module named 'pytest_cov'`

**解决**:
```bash
pip install pytest-cov
```

### 问题 4：某些测试被跳过

这是正常行为！被跳过的测试通常是：
- 需要真实 QC 软件（Gaussian/ORCA/xTB）
- 需要特定环境配置

查看跳过的原因：
```bash
pytest tests/test_qc_interface_v52.py -v --tb=short
```

输出示例：
```
tests/test_qc_interface_v52.py::TestQCInterface::test_real_gaussian SKIPPED (Requires real Gaussian environment)
```

---

## 📊 测试组织结构

```
tests/
├── conftest.py                      # pytest 配置，添加项目根目录到 sys.path
├── test_imports_step4_features.py   # 快速导入检查（CI gate）
├── test_s4_no_qc_execution.py       # S4 无 QC 测试（CI gate）
├── test_s4_*.py                     # S4 特征提取测试
├── test_m2_*.py                     # M2 合同/模式测试
├── test_m3_*.py                     # M3 模板/QC 测试
├── test_m4_*.py                     # M4 机制测试
├── test_mock_*.py                   # Mock 集成测试
├── test_qc_interface_v52.py         # QC 接口测试（真实软件，默认跳过）
├── test_orca_interface.py           # ORCA 接口测试（真实软件，默认跳过）
└── tmp_v2_2_test/                   # 提交的 S1 测试夹具（不要删除！）
    └── da_reaction/
        └── S1_Anchor/               # Diels-Alder 示例数据
```

---

## ✅ 预提交检查清单

提交代码前，确保：

- [ ] `pytest tests/test_imports_step4_features.py tests/test_s4_no_qc_execution.py -v` 通过
- [ ] `python ci/check_imports.py rph_core` 通过
- [ ] 新增代码遵循 AGENTS.md 中的代码规范
- [ ] 如果修改了 S4，运行 `pytest tests/test_s4_*.py -v`
- [ ] 文档已更新（README/AGENTS.md）

---

## 📚 相关文档

| 文档 | 内容 |
|------|------|
| [`AGENTS.md`](AGENTS.md) | 代码规范和开发指南 |
| [`tests/AGENTS.md`](tests/AGENTS.md) | 测试组织详情 |
| [`ci/AGENTS.md`](ci/AGENTS.md) | CI 集成指南 |
| [`MODIFICATION_REPORT_20250311.md`](MODIFICATION_REPORT_20250311.md) | 本次修改报告 |

---

**维护者**: QCcalc Team  
**最后更新**: 2025-03-11
