# 📋 ReactionProfileHunter v2.1 实施计划

**文档日期**: 2026-01-10  
**分析依据**: PROMOTE.md v2.1 串行架构 + 现有代码审查  
**目标**: 识别文档与代码的差距，制定可执行的修改计划

---

## 🚫 明确排除项 (Out of Scope)

以下功能不在本次实施计划范围内，需求未验证或优先级较低：

- ❌ NEB (Nudged Elastic Band) 方法支持
- ❌ 多参考态方法支持 (CASSCF, CASPT2 等)
- ❌ 自动机器学习模型训练
- ❌ Web UI 或图形界面
- ❌ 实时计算结果可视化
- ❌ 分布式计算集群支持
- ❌ 其他 QC 软件 (NWChem, Turbomole 等) 接口

---

## 🎯 第一个里程碑 (4h) - ORCAInterface 核心验证

### 目标
ORCAInterface 能生成输入文件 + 解析输出文件 (Mock 测试)

### 验证步骤
```bash
# 1. 创建骨架文件
mkdir -p rph_core/utils
touch rph_core/utils/orca_interface.py

# 2. 创建测试文件
mkdir -p tests
touch tests/test_orca_interface.py

# 3. 运行测试
pytest tests/test_orca_interface.py -v

# 4. 预期结果: 3个测试用例全部通过
```

### 测试用例清单
| 测试用例 | 功能 | 验证命令 | 预期输出 |
|:---|---|---|---|
| `test_generate_input_m062x()` | 生成 M062X 输入 | pytest -v -k m062x | 包含 `! M062X def2-TZVPP def2/J RIJCOSX` |
| `test_generate_input_with_solvent()` | 生成溶剂模型输入 | pytest -v -k solvent | 包含 `%cpcm` 块和 `smd true` |
| `test_parse_output_energy()` | 解析输出能量 | pytest -v -k parse | 从 fixture 解析出 `energy = -XXX.XXXXXX` |

### 骨架代码 (4小时内可完成)
```python
# rph_core/utils/orca_interface.py
from pathlib import Path
from dataclasses import dataclass
from typing import Optional
import re

@dataclass
class QCResult:
    """量子化学计算结果"""
    energy: float
    converged: bool
    coordinates: Optional[list] = None
    error_message: Optional[str] = None

class ORCAInterface:
    """ORCA 接口 - 高精度单点能计算"""
    
    def __init__(self, 
                 method: str = "M062X",
                 basis: str = "def2-TZVPP",
                 aux_basis: str = "def2/J",
                 nprocs: int = 16,
                 maxcore: int = 4000,
                 solvent: str = "acetone"):
        self.method = method
        self.basis = basis
        self.aux_basis = aux_basis
        self.nprocs = nprocs
        self.maxcore = maxcore
        self.solvent = solvent
    
    def _generate_input(self, xyz_file: Path, output_dir: Path) -> Path:
        """生成 ORCA 输入文件"""
        route = f"! {self.method} {self.basis} {self.aux_basis} RIJCOSX tightSCF"
        route += " noautostart miniprint nopop"
        
        cpcm_block = ""
        if self.solvent and self.solvent.upper() != "NONE":
            cpcm_block = f"""
%cpcm
   smd true
   SMDsolvent "{self.solvent}"
end
"""
        
        inp_content = f"""{route}
%maxcore {self.maxcore}
%pal nprocs {self.nprocs} end
{cpcm_block}
* xyzfile 0 1 {xyz_file.name}
"""
        inp_file = output_dir / f"{xyz_file.stem}.inp"
        inp_file.write_text(inp_content)
        return inp_file
    
    def _parse_output(self, out_file: Path) -> QCResult:
        """解析 ORCA 输出文件"""
        content = out_file.read_text()
        
        if "ORCA TERMINATED NORMALLY" not in content:
            return QCResult(energy=0.0, converged=False, 
                          error_message="ORCA 未正常终止")
        
        energy_match = re.search(r"FINAL SINGLE POINT ENERGY\s+([\-\d\.]+)", content)
        energy = float(energy_match.group(1)) if energy_match else 0.0
        
        return QCResult(energy=energy, converged=True)
```

```python
# tests/test_orca_interface.py
import pytest
from pathlib import Path
from rph_core.utils.orca_interface import ORCAInterface, QCResult

@pytest.fixture
def sample_output(tmp_path):
    """模拟 ORCA 输出文件"""
    out_file = tmp_path / "orca.out"
    out_file.write_text("""
ORCA CALCULATION DONE
FINAL SINGLE POINT ENERGY      -123.45678901
ORCA TERMINATED NORMALLY
""")
    return out_file

@pytest.fixture
def sample_xyz(tmp_path):
    """模拟 XYZ 文件"""
    xyz_file = tmp_path / "test.xyz"
    xyz_file.write_text("""3
comment
C    0.0  0.0  0.0
H    1.0  0.0  0.0
H   -1.0  0.0  0.0
""")
    return xyz_file

def test_generate_input_m062x(sample_xyz, tmp_path):
    """测试 M062X 输入生成"""
    orca = ORCAInterface(method="M062X", basis="def2-TZVPP")
    inp_file = orca._generate_input(sample_xyz, tmp_path)
    
    content = inp_file.read_text()
    assert "! M062X def2-TZVPP def2/J RIJCOSX" in content
    assert "%maxcore 4000" in content
    assert "%pal nprocs 16" in content

def test_generate_input_with_solvent(sample_xyz, tmp_path):
    """测试溶剂模型输入生成"""
    orca = ORCAInterface(solvent="acetone")
    inp_file = orca._generate_input(sample_xyz, tmp_path)
    
    content = inp_file.read_text()
    assert "%cpcm" in content
    assert 'SMDsolvent "acetone"' in content
    assert "smd true" in content

def test_parse_output_energy(sample_output):
    """测试能量解析"""
    orca = ORCAInterface()
    result = orca._parse_output(sample_output)
    
    assert result.converged is True
    assert result.energy == -123.45678901
    assert result.error_message is None
```

---

## 📊 任务依赖关系图

### 显式依赖链
```
ORCAInterface._generate_input()     ← 无依赖，可立即开发
           ↓
ORCAInterface._parse_output()       ← 无依赖，可立即开发
           ↓
ORCAInterface.single_point()        ← 需要上面两个方法
           ↓
SPMatrixBuilder._run_sp()           ← 需要 ORCAInterface.single_point()
           ↓
SPMatrixBuilder.run()               ← 需要 _run_sp()
           ↓
orchestrator.run_pipeline()         ← 需要 SPMatrixBuilder.run()
           ↓
FeatureMiner.run()                  ← 需要 orchestrator 传递 SPMatrixReport
```

### 关键依赖结论
- ✅ **可并行开发**: `ORCAInterface._generate_input()` 和 `_parse_output()` 可以同时开发
- ✅ **最小阻塞点**: 只要 `ORCAInterface.single_point()` 完成后，`SPMatrixBuilder` 就可以立即开始
- ✅ **S3.5 失败降级**: S3.5 失败不阻塞 S4，S4 将使用 L1 能量 (精度降低)

---

## 📋 总体差距矩阵

| 模块 | PROMOTE.md 要求 | 现有代码状态 | 差距级别 |
|:---:|---|---|:---:|
| **S1 ProductAnchor** | L1 Opt + L2 SP (ORCA) | 仅 CREST + XTB | 🔴 严重 |
| **S2 RetroScanner** | 完整 forming_bonds 传递 | 部分实现 | 🟡 中等 |
| **S3 TSOptimizer** | Berny + QST2 + IRC | 结构完整，IRC 未触发 | 🟢 轻微 |
| **S3.5 SP_Matrix** | ORCA RIJCOSX 批量 SP | ❌ 完全缺失 | 🔴 严重 |
| **S4 FeatureMiner** | 四模块解耦 + ORCA 能量 | 单文件，Gaussian 能量 | 🔴 严重 |
| **utils/qc_interface** | ORCA 接口支持 | ❌ 仅 XTB/Gaussian | 🔴 严重 |
| **Orchestrator** | S3.5 编排 + Batch | 无 S3.5，无 Batch | 🔴 严重 |

---

## 1. utils/qc_interface.py 详细实施计划

### 1.1 Session 拆分 (共 18h)

#### Day 1 - Session 1 (2h): `_generate_input()` + 测试
**目标**: 实现输入文件生成并测试

**任务**:
1. 创建 `ORCAInterface` 类骨架
2. 实现 `_generate_input()` 方法
3. 实现 `_is_double_hybrid()` 辅助方法
4. 编写 `test_generate_input_m062x()` 测试

**验收命令**:
```bash
pytest tests/test_orca_interface.py::test_generate_input_m062x -v
```

**预期输出**:
```
tests/test_orca_interface.py::test_generate_input_m062x PASSED [100%]
```

**验收标准**:
- [x] 生成的 `.inp` 文件包含 `! M062X def2-TZVPP def2/J RIJCOSX`
- [x] 包含 `%maxcore 4000`
- [x] 包含 `%pal nprocs 16`

---

#### Day 1 - Session 2 (2h): `_parse_output()` + 测试
**目标**: 实现输出文件解析并测试

**任务**:
1. 实现 `_parse_output()` 方法
2. 添加错误处理 (检查正常终止)
3. 编写 `test_parse_output_energy()` 测试
4. 创建多个 fixture 输出文件 (正常/失败/能量解析失败)

**验收命令**:
```bash
pytest tests/test_orca_interface.py::test_parse_output_energy -v
pytest tests/test_orca_interface.py::test_parse_output_failed -v
```

**预期输出**:
```
tests/test_orca_interface.py::test_parse_output_energy PASSED
tests/test_orca_interface.py::test_parse_output_failed PASSED
```

**验收标准**:
- [x] 从 fixture 输出解析出 `energy = -123.45678901`
- [x] 失败时返回 `converged=False`
- [x] 失败时包含 `error_message`

---

#### Day 1 - Session 3 (2h): `_find_orca_binary()` + `_run_orca()` 骨架
**目标**: 实现 ORCA 可执行文件定位和运行骨架

**任务**:
1. 实现 `_find_orca_binary()` 方法
2. 实现 `_run_orca()` 骨架 (使用 subprocess)
3. Mock 测试 (无需真实 ORCA 环境)
4. 添加超时和错误处理

**验收命令**:
```bash
pytest tests/test_orca_interface.py::test_find_orca_binary -v
pytest tests/test_orca_interface.py::test_run_orca_mock -v
```

**预期输出**:
```
tests/test_orca_interface.py::test_find_orca_binary SKIPPED (no ORCA)
tests/test_orca_interface.py::test_run_orca_mock PASSED
```

**验收标准**:
- [x] 能从环境变量或配置文件找到 ORCA 路径
- [x] Mock 测试验证 subprocess 调用正确
- [x] 无 ORCA 环境时优雅降级

---

#### Day 1 - Session 4 (2h): `single_point()` 整合 + Mock 测试
**目标**: 整合所有方法实现完整单点能计算

**任务**:
1. 实现 `single_point()` 方法
2. 整合输入生成 → 运行 → 解析流程
3. 编写 `test_single_point_mock()` 端到端测试
4. 添加日志记录

**验收命令**:
```bash
pytest tests/test_orca_interface.py::test_single_point_mock -v
```

**预期输出**:
```
tests/test_orca_interface.py::test_single_point_mock PASSED
```

**验收标准**:
- [x] 端到端测试通过
- [x] 返回 `QCResult` 对象
- [x] 日志记录完整

---

#### Day 2 - Session 5 (2h): SMD 溶剂模型 + 测试
**目标**: 实现完整 SMD 溶剂模型支持

**任务**:
1. 完善 `_generate_input()` 中的 CPCM/SMD 块
2. 支持多种溶剂 (acetone, water, toluene 等)
3. 编写溶剂模型测试用例
4. 测试 `solvent=None` 行为

**验收命令**:
```bash
pytest tests/test_orca_interface.py::test_solvent_model -v
pytest tests/test_orca_interface.py::test_solvent_none -v
```

**预期输出**:
```
tests/test_orca_interface.py::test_solvent_model PASSED
tests/test_orca_interface.py::test_solvent_none PASSED
```

**验收标准**:
- [x] 生成的输入包含正确的 `%cpcm` 块
- [x] `solvent=None` 时不生成溶剂块
- [x] 支持至少 5 种常见溶剂

---

#### Day 2 - Session 6 (2h): 双杂化泛函自动匹配 + 测试
**目标**: 实现双杂化泛函自动添加 `/C` 辅助基组

**任务**:
1. 完善 `_is_double_hybrid()` 方法
2. 识别 PWPB95, DSD-PBEP86, B2PLYP 等双杂化泛函
3. 自动添加 `def2-TZVPP/C` 辅助基组
4. 编写双杂化泛函测试

**验收命令**:
```bash
pytest tests/test_orca_interface.py::test_double_hybrid -v
```

**预期输出**:
```
tests/test_orca_interface.py::test_double_hybrid PASSED
```

**验收标准**:
- [x] PWPB95 自动添加 `def2-TZVPP/C`
- [x] B3LYP (非双杂化) 不添加 `/C`
- [x] 路由行格式正确

---

#### Day 2 - Session 7 (1h): 真实 ORCA 集成测试
**目标**: 在真实 ORCA 环境中验证 (如果可用)

**任务**:
1. 创建 `test_real_sp.py` (标记为 `requires_orca`)
2. 使用简单分子 (如 H2O) 进行真实计算
3. 验证输出能量合理性
4. 测试超时和异常处理

**验收命令**:
```bash
pytest tests/test_orca_interface.py::test_real_sp -v -m requires_orca
```

**预期输出** (如果有 ORCA):
```
tests/test_orca_interface.py::test_real_sp PASSED
```

**预期输出** (如果无 ORCA):
```
tests/test_orca_interface.py::test_real_sp SKIPPED (requires ORCA binary)
```

**验收标准**:
- [x] 真实计算返回合理能量
- [x] 输入/输出文件格式正确
- [x] 异常情况正确处理

---

### 1.4 可执行版验收标准

#### 1. 输入生成测试 (无需 ORCA)
```bash
pytest tests/test_orca_interface.py::test_generate_input -v
```

**预期输出**: 生成的 `.inp` 文件包含:
```
! M062X def2-TZVPP def2/J RIJCOSX tightSCF
%maxcore 4000
* xyzfile 0 1 input.xyz
```

#### 2. 输出解析测试 (无需 ORCA)
```bash
pytest tests/test_orca_interface.py::test_parse_output -v
```

**预期输出**: 从 fixture 输出解析出 `energy = -123.45678901`

#### 3. 集成测试 (需要 ORCA)
```bash
pytest tests/test_orca_interface.py::test_real_sp -v -m requires_orca
```

**预期输出**: 返回 `QCResult(energy=-76.367, converged=True)`

---

## 2. step1_anchor/anchor_manager.py 修改计划

### 2.1 问题诊断

**PROMOTE.md 要求** (Section 2.2, 3.2):
- S1 输出包含 L2 SP 能量: `e_product_l2 = run_orca_sp(...)`
- 产物锚定后立即执行高精度 SP

**现有代码**: 
```python
# Step 1.4: DFT 优化 (可选)
if self.dft_config.get('method'):
    self.logger.info("DFT 优化启用，执行高质量精修...")
    # 这里可以调用 GaussianInterface 进行优化
    # 暂时保留 XTB 结果作为输出

self.logger.info(f"Step 1 完成: {product_min_xyz}")
return product_min_xyz  # ❌ 缺失 L2 能量返回
```

### 2.2 Session 拆分 (共 7h)

#### Session 1 (2h): 新增 `_run_l2_sp()` 方法
**任务**:
1. 导入 `ORCAInterface`
2. 实现 `_run_l2_sp()` 方法
3. 添加错误处理和日志
4. 编写测试

**验收命令**:
```bash
pytest tests/test_anchor_manager.py::test_run_l2_sp -v
```

**验收标准**:
- [x] 能调用 ORCAInterface 执行单点能计算
- [x] 失败时抛出清晰的异常
- [x] 返回正确的能量值

---

#### Session 2 (1h): 修改 `run()` 返回值
**任务**:
1. 修改 `run()` 方法返回 `(product_xyz, e_product_l2)` 元组
2. 更新调用方代码
3. 确保向后兼容

**验收命令**:
```bash
pytest tests/test_anchor_manager.py::test_run_returns_tuple -v
```

**验收标准**:
- [x] 返回类型为 `Tuple[Path, float]`
- [x] 第一元素为产物 XYZ 路径
- [x] 第二元素为 L2 能量 (float)

---

#### Session 3 (3h): 补全 DFT 优化分支
**任务**:
1. 实现可选的 L1 Gaussian Opt
2. 集成到现有工作流
3. 验证优化和 SP 的正确顺序

**验收命令**:
```bash
pytest tests/test_anchor_manager.py::test_dft_optimization -v
```

**验收标准**:
- [x] L1 Opt → L2 SP 流程正确
- [x] 优化结果用于后续 SP 计算

---

#### Session 4 (1h): 保存元数据
**任务**:
1. 实现 `_save_metadata()` 方法
2. 保存 L2 能量到 JSON
3. 保存计算方法和基组信息

**验收命令**:
```bash
pytest tests/test_anchor_manager.py::test_save_metadata -v
```

**验收标准**:
- [x] `s1_metadata.json` 文件存在
- [x] 包含 `e_product_l2` 字段
- [x] 包含 `l2_method` 字段

---

### 2.3 详细实现规范

```python   
# anchor_manager.py 修改

from rph_core.utils.qc_interface import ORCAInterface  # 新增导入

class ProductAnchor(LoggerMixin):
    
    def __init__(self, config: dict):
        # ... 现有代码 ...
        
        # 新增: L2 SP 配置
        self.l2_config = config.get('high_level_sp', {})
        self.orca = ORCAInterface(
            method=self.l2_config.get('method', 'M062X'),
            basis=self.l2_config.get('basis', 'def2-TZVPP'),
            nprocs=self.l2_config.get('nprocs', 16),
            maxcore=self.l2_config.get('maxcore', 4000),
            solvent=self.l2_config.get('solvent', 'acetone')
        )
    
    def run(self, smiles: str, output_dir: Path) -> Tuple[Path, float]:
        """
        执行产物锚定工作流
        
        Returns:
            (product_min_xyz, e_product_l2) - 产物XYZ路径和L2能量
        """
        # ... 现有 Step 1.1-1.4 代码 ...
        
        # Step 1.5: [NEW] L2 高精度单点能 (ORCA)
        self.logger.info("执行 L2 高精度单点能计算 (ORCA)...")
        e_product_l2 = self._run_l2_sp(product_min_xyz, output_dir / "L2_SP")
        
        self.logger.info(f"Step 1 完成: {product_min_xyz}")
        self.logger.info(f"  L2 能量: {e_product_l2:.8f} Hartree")
        
        # 保存元数据
        self._save_metadata(output_dir, product_min_xyz, e_product_l2)
        
        return product_min_xyz, e_product_l2
    
    def _run_l2_sp(self, xyz_file: Path, output_dir: Path) -> float:
        """执行 L2 级别 ORCA 单点能计算"""
        result = self.orca.single_point(xyz_file, output_dir)
        
        if not result.converged:
            self.logger.error(f"L2 SP 计算失败: {result.error_message}")
            raise RuntimeError(f"L2 SP 失败: {result.error_message}")
        
        return result.energy
    
    def _save_metadata(self, output_dir: Path, xyz_path: Path, l2_energy: float):
        """保存 S1 元数据到 JSON"""
        import json
        meta = {
            "product_xyz": str(xyz_path),
            "e_product_l2": l2_energy,
            "l2_method": f"{self.l2_config.get('method')}/{self.l2_config.get('basis')}"
        }
        meta_file = output_dir / "s1_metadata.json"
        meta_file.write_text(json.dumps(meta, indent=2))
```

### 2.4 可执行版验收标准

```bash
pytest tests/test_anchor_manager.py::test_run_returns_tuple -v
pytest tests/test_anchor_manager.py::test_save_metadata -v
```

**预期输出**:
- `ProductAnchor.run()` 返回 `(Path, float)` 元组
- `s1_metadata.json` 包含:
  ```json
  {
    "product_xyz": "/path/to/product.xyz",
    "e_product_l2": -123.45678901,
    "l2_method": "M062X/def2-TZVPP"
  }
  ```

---

## 3. 新增 step3_5_sp/ 目录 (全新模块)

### 3.1 问题诊断

**PROMOTE.md 要求** (Section 3.3.5):
- Step 3.5 对 S1-S3 的所有关键节点进行 ORCA L2 SP
- 覆盖矩阵: Product, Reactant, TS_Final, Frag_A, Frag_B
- 输出 SP_Matrix_Report 供 S4 读取

**现有代码**: ❌ 完全缺失

### 3.2 新建文件清单

```
rph_core/steps/
└── step3_5_sp/           # [NEW]
    ├── __init__.py
    ├── sp_matrix.py      # 主控: SPMatrixBuilder 类
    └── sp_report.py      # 报告: SPMatrixReport 数据结构
```

### 3.3 Session 拆分 (共 8h)

#### Session 1 (2h): 创建 `SPMatrixReport` 数据结构
**任务**:
1. 创建 `sp_report.py`
2. 定义 `SPMatrixReport` dataclass
3. 实现 `to_dict()` 方法
4. 编写测试

**验收命令**:
```bash
pytest tests/test_sp_report.py::test_sp_report_creation -v
```

**验收标准**:
- [x] `SPMatrixReport` 包含所有必需字段
- [x] `to_dict()` 返回正确的字典格式
- [x] 序列化/反序列化正常

---

#### Session 2 (3h): 实现 `SPMatrixBuilder` 核心逻辑
**任务**:
1. 创建 `sp_matrix.py`
2. 实现 `SPMatrixBuilder.__init__()`
3. 实现 `_run_sp()` 方法
4. 实现 `run()` 主方法
5. 编写 Mock 测试

**验收命令**:
```bash
pytest tests/test_sp_matrix.py::test_run_sp -v
pytest tests/test_sp_matrix.py::test_run -v
```

**验收标准**:
- [x] 能正确调用 ORCAInterface
- [x] 能处理 5 个节点的计算
- [x] 返回正确的 SPMatrixReport

---

#### Session 3 (3h): 实现复用 S1 能量逻辑
**任务**:
1. 实现 `e_product_l2` 复用逻辑
2. 添加日志记录
3. 实现 `_save_report()` 方法
4. 编写集成测试

**验收命令**:
```bash
pytest tests/test_sp_matrix.py::test_reuse_product_energy -v
```

**验收标准**:
- [x] 当 `e_product_l2` 不为 None 时跳过计算
- [x] 报告正确保存为 JSON
- [x] 日志清晰记录复用行为

---

### 3.4 详细实现规范

```python
# step3_5_sp/sp_report.py

from dataclasses import dataclass
from typing import Dict, Optional

@dataclass
class SPMatrixReport:
    """SP 矩阵报告"""
    e_product: float       # L2 产物能量
    e_reactant: float      # L2 底物能量
    e_ts_final: float      # L2 TS 能量
    e_frag_a_ts: float     # L2 片段A在TS几何下能量
    e_frag_b_ts: float     # L2 片段B在TS几何下能量
    e_frag_a_relaxed: Optional[float] = None  # 可选: 松弛后能量
    e_frag_b_relaxed: Optional[float] = None
    
    def to_dict(self) -> Dict[str, float]:
        return {
            "e_product_l2": self.e_product,
            "e_reactant_l2": self.e_reactant,
            "e_ts_final_l2": self.e_ts_final,
            "e_frag_a_ts_l2": self.e_frag_a_ts,
            "e_frag_b_ts_l2": self.e_frag_b_ts,
            "e_frag_a_relaxed_l2": self.e_frag_a_relaxed,
            "e_frag_b_relaxed_l2": self.e_frag_b_relaxed
        }
```

```python
# step3_5_sp/sp_matrix.py

from pathlib import Path
from typing import Optional
import logging
import json

from rph_core.utils.qc_interface import ORCAInterface
from rph_core.steps.step3_5_sp.sp_report import SPMatrixReport
from rph_core.utils.log_manager import LoggerMixin

logger = logging.getLogger(__name__)

class SPMatrixBuilder(LoggerMixin):
    """
    SP 矩阵构建器 (Step 3.5)
    
    职责:
    1. 收集 S1-S3 产生的关键几何节点
    2. 对每个节点执行 ORCA L2 SP
    3. 生成 SPMatrixReport 供 S4 使用
    """
    
    def __init__(self, config: dict):
        self.config = config
        
        # 初始化 ORCA 接口
        self.orca = ORCAInterface(
            method=config.get('method', 'M062X'),
            basis=config.get('basis', 'def2-TZVPP'),
            aux_basis=config.get('aux_basis', 'def2/J'),
            nprocs=config.get('nprocs', 16),
            maxcore=config.get('maxcore', 4000),
            solvent=config.get('solvent', 'acetone')
        )
        
        self.logger.info(f"SPMatrixBuilder 初始化: {self.orca.method}/{self.orca.basis}")
    
    def run(
        self,
        product_xyz: Path,
        reactant_xyz: Path,
        ts_final_xyz: Path,
        frag_a_xyz: Path,
        frag_b_xyz: Path,
        output_dir: Path,
        e_product_l2: Optional[float] = None  # S1 可能已计算
    ) -> SPMatrixReport:
        """
        执行 SP 矩阵构建
        
        Args:
            product_xyz: 产物结构 (来自 S1)
            reactant_xyz: 底物结构 (来自 S2)
            ts_final_xyz: TS 结构 (来自 S3)
            frag_a_xyz: 片段A在TS几何下
            frag_b_xyz: 片段B在TS几何下
            output_dir: 输出目录
            e_product_l2: 如果 S1 已计算，可直接传入避免重复
        
        Returns:
            SPMatrixReport 对象
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        self.logger.info("Step 3.5 开始: 高精度 SP 矩阵构建")
        
        # 1. Product L2 SP (可能已在 S1 计算)
        if e_product_l2 is not None:
            self.logger.info("  [Product] 使用 S1 已计算的 L2 能量")
        else:
            self.logger.info("  [Product] 执行 L2 SP...")
            e_product_l2 = self._run_sp(product_xyz, output_dir / "product_sp")
        
        # 2. Reactant L2 SP
        self.logger.info("  [Reactant] 执行 L2 SP...")
        e_reactant_l2 = self._run_sp(reactant_xyz, output_dir / "reactant_sp")
        
        # 3. TS_Final L2 SP
        self.logger.info("  [TS_Final] 执行 L2 SP...")
        e_ts_l2 = self._run_sp(ts_final_xyz, output_dir / "ts_sp")
        
        # 4. Fragment A L2 SP (at TS geometry)
        self.logger.info("  [Frag_A] 执行 L2 SP...")
        e_frag_a_l2 = self._run_sp(frag_a_xyz, output_dir / "frag_a_sp")
        
        # 5. Fragment B L2 SP (at TS geometry)
        self.logger.info("  [Frag_B] 执行 L2 SP...")
        e_frag_b_l2 = self._run_sp(frag_b_xyz, output_dir / "frag_b_sp")
        
        report = SPMatrixReport(
            e_product=e_product_l2,
            e_reactant=e_reactant_l2,
            e_ts_final=e_ts_l2,
            e_frag_a_ts=e_frag_a_l2,
            e_frag_b_ts=e_frag_b_l2
        )
        
        # 保存报告
        self._save_report(report, output_dir / "sp_matrix_report.json")
        
        self.logger.info("Step 3.5 完成: SP 矩阵构建成功")
        return report
    
    def _run_sp(self, xyz_file: Path, output_dir: Path) -> float:
        """执行单个 ORCA SP 计算
        
        支持自动降级到 Gaussian (如果配置启用)
        """
        try:
            result = self.orca.single_point(xyz_file, output_dir)
            if not result.converged:
                raise RuntimeError(f"ORCA SP 失败: {xyz_file} - {result.error_message}")
            return result.energy
            
        except Exception as orca_error:
            # 检查是否启用 Gaussian fallback
            if self.config.get('fallback_to_gaussian', True):
                self.logger.warning(f"ORCA 不可用，降级使用 Gaussian: {orca_error}")
                try:
                    from rph_core.utils.qc_interface import GaussianInterface
                    self.gaussian = GaussianInterface(
                        method=self.config.get('fallback_method', 'B3LYP'),
                        basis=self.config.get('fallback_basis', 'def2-TZVP'),
                        nprocs=self.config.get('nprocs', 16)
                    )
                    result = self.gaussian.single_point(xyz_file, output_dir)
                    if not result.converged:
                        raise RuntimeError(f"Gaussian SP 也失败: {xyz_file}")
                    self.logger.info(f"  [Fallback] 成功使用 Gaussian 完成 SP")
                    return result.energy
                except Exception as gaussian_error:
                    self.logger.error(f"Gaussian fallback 也失败: {gaussian_error}")
                    raise
            else:
                # 未启用 fallback，直接抛出异常
                raise
    
    def _save_report(self, report: SPMatrixReport, output_path: Path):
        """保存 SP 矩阵报告"""
        output_path.write_text(json.dumps(report.to_dict(), indent=2))
        self.logger.info(f"  SP 矩阵报告已保存: {output_path}")
```

### 3.5 可执行版验收标准

```bash
pytest tests/test_sp_matrix.py::test_sp_matrix_builder_run -v
```

**预期输出**:
- `SPMatrixBuilder.run()` 成功对 5 个节点执行 ORCA SP
- `sp_matrix_report.json` 包含所有 L2 能量:
  ```json
  {
    "e_product_l2": -123.45678901,
    "e_reactant_l2": -456.78901234,
    "e_ts_final_l2": -345.67890123,
    "e_frag_a_ts_l2": -234.56789012,
    "e_frag_b_ts_l2": -123.45678901,
    "e_frag_a_relaxed_l2": null,
    "e_frag_b_relaxed_l2": null
  }
  ```
- 能正确复用 S1 已计算的产物 L2 能量

---

## 4. step4_features/ 模块重构计划

### 4.1 问题诊断

**PROMOTE.md 要求** (Section 3.4):
- 四模块解耦: `electronic.py`, `fmo_reactivity.py`, `steric_geometry.py`, `entropy.py`
- 使用 S3.5 的 L2 能量计算畸变能
- 支持 HBDE 计算

**现有代码**:
- `feature_miner.py`: 单文件实现，依赖 Gaussian 能量
- `fragment_extractor.py`: 已实现，但使用 Gaussian
- `distortion_calculator.py`: 已实现基础公式

### 4.2 Session 拆分 (共 13h)

#### Session 1 (3h): 新增 `electronic.py` (HBDE 计算)
**任务**:
1. 创建 `electronic.py`
2. 实现 `calculate_hbde()` 函数
3. 实现 `calculate_nics_index()` 函数
4. 编写测试

**验收命令**:
```bash
pytest tests/test_electronic.py::test_calculate_hbde -v
```

**验收标准**:
- [x] HBDE 计算公式正确: `(E_cation + E_anion - E_precursor) * 627.509`
- [x] 返回单位为 kcal/mol
- [x] 测试用例覆盖正常边界情况

---

#### Session 2-3 (4h): 新增 `steric_geometry.py` (畸变能 + Sterimol)
**任务**:
1. 创建 `steric_geometry.py`
2. 集成 `distortion_calculator.py`
3. 实现基于 L2 能量的畸变能计算
4. (可选) 集成 morfeus 计算 Sterimol/Vbur

**验收命令**:
```bash
pytest tests/test_steric_geometry.py::test_distortion_energy_l2 -v
```

**验收标准**:
- [x] 正确计算 E_distortion_A, E_distortion_B
- [x] 正确计算 E_interaction
- [x] 支持基于 SPMatrixReport 的计算

---

#### Session 4 (2h): 修改 `feature_miner.py` 读取 SP 矩阵报告
**任务**:
1. 修改 `run()` 方法签名，添加 `sp_matrix_report` 参数
2. 实现 `_extract_features_from_sp_matrix()` 方法
3. 实现降级逻辑 (无 SP 矩阵时使用 L1 能量)
4. 编写测试

**验收命令**:
```bash
pytest tests/test_feature_miner.py::test_extract_features_with_sp_matrix -v
pytest tests/test_feature_miner.py::test_degradation_to_l1 -v
```

**验收标准**:
- [x] 有 SPMatrixReport 时使用 L2 能量
- [x] 无 SPMatrixReport 时使用 L1 能量 (降级)
- [x] 输出 CSV 包含 `*_L2` 后缀的高精度特征

---

#### Session 5-6 (4h): 其他模块实现 (可选)
**任务**:
1. 新增 `fmo_reactivity.py` (HOMO/LUMO 指数)
2. 新增 `entropy.py` (振动熵提取)
3. 集成测试

---

### 4.3 关键修改: feature_miner.py

```python
# feature_miner.py 关键修改

from rph_core.steps.step3_5_sp.sp_report import SPMatrixReport

class FeatureMiner(LoggerMixin):
    
    def run(
        self,
        ts_final: Path,
        reactant: Path,
        product: Path,
        output_dir: Path,
        forming_bonds: Optional[Tuple[Tuple[int, int], Tuple[int, int]]] = None,
        fragment_indices: Optional[Tuple[List[int], List[int]]] = None,
        sp_matrix_report: Optional[SPMatrixReport] = None  # [NEW] S3.5 报告
    ) -> Path:
        """
        执行特征提取
        
        [NEW] 优先使用 sp_matrix_report 中的 L2 能量
        """
        # ...
        
        if sp_matrix_report:
            self.logger.info("使用 S3.5 高精度 SP 矩阵计算特征")
            features = self._extract_features_from_sp_matrix(
                sp_matrix_report, ts_final, reactant, product, forming_bonds
            )
        else:
            self.logger.warning("未提供 SP 矩阵，使用优化级别能量（精度较低）")
            features = self._extract_features(
                ts_final, reactant, product, forming_bonds
            )
        
        # ...
    
    def _extract_features_from_sp_matrix(
        self,
        sp_report: SPMatrixReport,
        ts_file: Path,
        reactant_file: Path,
        product_file: Path,
        forming_bonds
    ) -> Dict[str, float]:
        """基于 L2 SP 矩阵提取高精度特征"""
        features = {}
        
        # 活化能 (L2 级别)
        features["dG_activation_L2"] = self.dist_calc.calculate_activation_energy(
            sp_report.e_ts_final, sp_report.e_reactant
        )
        
        # 反应能 (L2 级别)
        features["dG_reaction_L2"] = self.dist_calc.calculate_reaction_energy(
            sp_report.e_product, sp_report.e_reactant
        )
        
        # 双片段畸变能 (L2 级别)
        distortion_results = self.dist_calc.calculate_distortion_interaction(
            e_ts=sp_report.e_ts_final,
            e_fragment_a_at_ts=sp_report.e_frag_a_ts,
            e_fragment_b_at_ts=sp_report.e_frag_b_ts,
            e_fragment_a_relaxed=sp_report.e_frag_a_relaxed or sp_report.e_frag_a_ts,
            e_fragment_b_relaxed=sp_report.e_frag_b_relaxed or sp_report.e_frag_b_ts
        )
        features.update({
            "E_distortion_A_L2": distortion_results["e_distortion_a"],
            "E_distortion_B_L2": distortion_results["e_distortion_b"],
            "E_distortion_total_L2": distortion_results["e_distortion_total"],
            "E_interaction_L2": distortion_results["e_interaction"]
        })
        
        # 几何特征 (仍从 TS 结构提取)
        # ...
        
        return features
```

### 4.4 可执行版验收标准

```bash
pytest tests/test_feature_miner.py::test_extract_features_with_sp_matrix -v
```

**预期输出**: CSV 文件包含以下 L2 精度特征:
```
dG_activation_L2,dG_reaction_L2,E_distortion_A_L2,E_distortion_B_L2,E_interaction_L2
15.234, -8.567, 5.123, 3.456, 6.655
```

---

## 5. orchestrator.py 修改计划

### 5.1 问题诊断

**PROMOTE.md 要求**:
- 编排 S3.5 (SP Matrix)
- 支持批量处理 (`run_batch()`)
- 支持断点续传

**现有代码**: 
- 无 S3.5 调用
- 无 batch 支持
- 无 checkpoint 机制

### 5.2 Session 拆分 (共 11h)

#### Session 1 (1h): 新增 S3.5 引擎属性
**任务**:
1. 添加 `s35_engine` property
2. 实现懒加载
3. 编写测试

**验收命令**:
```bash
pytest tests/test_orchestrator.py::test_s35_engine_property -v
```

**验收标准**:
- [x] 首次访问时正确初始化 SPMatrixBuilder
- [x] 多次访问返回同一实例

---

#### Session 2-3 (2h): 在 `run_pipeline()` 中插入 S3.5 调用
**任务**:
1. 在 S3 完成后插入 S3.5 代码
2. 实现错误处理 (S3.5 失败不阻塞 S4)
3. 实现降级逻辑标记
4. 编写测试

**验收命令**:
```bash
pytest tests/test_orchestrator.py::test_pipeline_with_s35 -v
```

**验收标准**:
- [x] S3.5 在 S3 完成后自动执行
- [x] S3.5 失败时正确降级
- [x] 降级标记正确传递到 S4

---

#### Session 4 (1h): 修改 S4 调用以传递 `sp_matrix_report`
**任务**:
1. 修改 S4 调用代码
2. 传递 `sp_matrix_report` 参数
3. 验证数据流

**验收命令**:
```bash
pytest tests/test_orchestrator.py::test_s4_receives_sp_matrix -v
```

**验收标准**:
- [x] S4 正确接收 `SPMatrixReport`
- [x] 数据传递完整

---

#### Session 5-6 (4h): 新增 `run_batch()` 方法
**任务**:
1. 实现 `run_batch()` 方法
2. 支持多 SMILES 批量处理
3. 实现并发控制
4. 编写测试

**并发控制策略**:
- 使用 `concurrent.futures.ProcessPoolExecutor`
- 默认 `max_workers = min(4, cpu_count())`
- 每个任务独立工作目录，避免文件冲突
- 使用 `tqdm` 显示进度条
- 单个任务失败不影响其他任务

**验收命令**:
```bash
pytest tests/test_orchestrator.py::test_run_batch -v
```

**验收标准**:
- [x] 能处理多个 SMILES
- [x] 每个任务独立运行
- [x] 错误不阻塞其他任务
- [x] 支持配置并发数

---

#### Session 7 (3h): 实现 JSON checkpoint 机制
**任务**:
1. 实现 `_save_checkpoint()` 方法
2. 实现 `_load_checkpoint()` 方法
3. 实现断点续传逻辑
4. 编写测试

**验收命令**:
```bash
pytest tests/test_orchestrator.py::test_checkpoint_save_load -v
```

**验收标准**:
- [x] checkpoint 文件正确保存
- [x] 能从 checkpoint 恢复执行
- [x] 已完成步骤被跳过

---

### 5.3 关键修改

```python
# orchestrator.py 关键修改

class ReactionProfileHunter:
    
    @property
    def s35_engine(self):
        """Step 3.5 引擎 (懒加载)"""
        if self._s35_engine is None:
            from rph_core.steps.step3_5_sp import SPMatrixBuilder
            self._s35_engine = SPMatrixBuilder(self.config.get('step3_5', {}))
            self.logger.debug("✓ Step 3.5 (SPMatrixBuilder) 已初始化")
        return self._s35_engine
    
    def run_pipeline(self, product_smiles: str, work_dir: Path, skip_steps: list = None):
        # ... 现有 S1-S3 代码 ...
        
        # === Step 3.5: 高精度 SP 矩阵 [NEW] ===
        sp_matrix_report = None
        if 's3.5' not in skip_steps and result.ts_final_xyz:
            try:
                self.logger.info(">>> Step 3.5: 构建高精度 SP 矩阵...")
                
                # 需要先提取片段
                from rph_core.steps.step4_features.fragment_extractor import FragmentExtractor
                frag_extractor = FragmentExtractor(self.config.get('step4', {}))
                
                frag_dir = work_dir / "S35_SP" / "fragments"
                frag_a_xyz, frag_b_xyz = frag_extractor.extract_fragment_xyz(
                    result.ts_final_xyz,
                    result.forming_bonds,  # 需要 forming_bonds 确定切分位置
                    frag_dir
                )
                
                sp_matrix_report = self.s35_engine.run(
                    product_xyz=result.product_xyz,
                    reactant_xyz=result.reactant_xyz,
                    ts_final_xyz=result.ts_final_xyz,
                    frag_a_xyz=frag_a_xyz,
                    frag_b_xyz=frag_b_xyz,
                    output_dir=work_dir / "S35_SP",
                    e_product_l2=result.e_product_l2  # 复用 S1 的 L2 能量
                )
                
                self.logger.info("    ✓ SP 矩阵构建完成")
                
            except Exception as e:
                result.error_step = "Step3.5_SPMatrix"
                result.error_message = str(e)
                self.logger.error(f"Step 3.5 失败: {e}", exc_info=True)
                # SP 失败不阻塞 S4，但 S4 将使用低精度能量
        
        # === Step 4: 特征挖掘 ===
        if 's4' not in skip_steps and result.ts_final_xyz:
            try:
                self.logger.info(">>> Step 4: 提取物理有机特征...")
                features_csv = self.s4_engine.run(
                    ts_final=result.ts_final_xyz,
                    reactant=result.reactant_xyz,
                    product=result.product_xyz,
                    output_dir=work_dir / "S4_Data",
                    forming_bonds=result.forming_bonds,
                    sp_matrix_report=sp_matrix_report  # [NEW] 传递 L2 能量
                )
                # ...
```

### 5.4 可执行版验收标准

```bash
pytest tests/test_orchestrator.py::test_pipeline_with_s35 -v
```

**预期输出**:
- S3.5 在 S3 完成后自动执行
- S4 能接收并使用 `SPMatrixReport`
- 输出目录结构符合 `S35_SP/` 规范

---

## 6. config/defaults.yaml 补全计划

### 6.1 问题诊断

**当前问题**: 嵌套太深，且没区分必须 vs 可选

**修改**: 简化配置，明确标记必须配置项

### 6.2 简化版配置

```yaml
# ========================
# Step 3.5: High-Precision SP Matrix
# ========================
step3_5:
  # 必须配置 (无默认值)
  orca_binary: null  # 必须设置！例如: /opt/orca/orca
  
  # 可选配置 (有默认值)
  method: M062X           # 默认
  basis: def2-TZVPP       # 默认
  aux_basis: def2/J       # 默认
  nprocs: 16              # 默认
  maxcore: 4000           # 默认
  solvent: acetone        # 默认，设为 null 则不加溶剂
  
  # 双杂化泛函额外配置
  double_hybrid_c_basis: def2-TZVPP/C  # 自动检测双杂化泛函时使用
  
# ========================
# ORCA 全局配置
# ========================
orca:
  # 必须配置
  binary_path: null  # 或环境变量 $ORCA_BIN
  
  # 可选配置
  scratch_dir: /tmp/orca_scratch  # 默认临时目录
  
# ========================
# Step 1: Product Anchor (更新)
# ========================
step1:
  crest:
    gfn_level: 2
    solvent: acetone
    energy_window: 6.0
    threads: 16
  dft:
    method: B3LYP
    basis: def2-SVP
    dispersion: GD3BJ
    verify_no_imaginary: true
  high_level_sp:  # L2 SP 配置
    method: M062X
    basis: def2-TZVPP
    nprocs: 16
    maxcore: 4000
    solvent: acetone
```

### 6.3 配置验证脚本

```python
# utils/config_validator.py

def validate_config(config: dict) -> Tuple[bool, List[str]]:
    """
    验证配置文件的完整性
    
    Returns:
        (is_valid, error_messages)
    """
    errors = []
    
    # 检查 ORCA 必须配置
    if 'orca' in config:
        if config['orca'].get('binary_path') is None:
            errors.append("orca.binary_path 未设置 (必须配置)")
    
    # 检查 Step 3.5 ORCA 路径
    if 'step3_5' in config:
        if config['step3_5'].get('orca_binary') is None:
            errors.append("step3_5.orca_binary 未设置 (必须配置)")
    
    return len(errors) == 0, errors
```

**验收命令**:
```bash
python -c "from rph_core.utils.config_validator import validate_config; import yaml; config = yaml.safe_load(open('config/defaults.yaml')); print(validate_config(config))"
```

**预期输出**:
```
(False, ['step3_5.orca_binary 未设置 (必须配置)'])
```

---

## 7. 错误处理与降级策略

### 7.1 S3.5 失败降级逻辑

```python
# orchestrator.py 中应该有的逻辑

# S3.5 失败后的处理
if sp_matrix_report is None:
    self.logger.warning("S3.5 失败，S4 将使用 L1 能量 (精度降低)")
    result.sp_quality = "L1"  # 标记精度级别
    result.e_product_l2 = None
else:
    result.sp_quality = "L2"
    result.e_product_l2 = sp_matrix_report.e_product

# S4 中应该有的逻辑
if sp_quality == "L1":
    features["energy_quality"] = "L1_ONLY"
    features["dG_activation_L2"] = None  # 明确标记不可用
    features["dG_reaction_L2"] = None
    features["E_distortion_A_L2"] = None
    features["E_distortion_B_L2"] = None
    features["E_interaction_L2"] = None
```

### 7.2 降级策略清单

| 场景 | 降级行为 | 标记字段 | 影响特征 |
|:---|---|:---:|:---|
| ORCA 不可用 | 跳过 S3.5，使用 L1 能量 | `sp_quality="L1"` | `*_L2` = `None` |
| S3.5 部分失败 | 使用已完成的能量 | `sp_quality="L1_PARTIAL"` | 失败节点 `*_L2` = `None` |
| ORCA 计算超时 | 记录错误，跳过该节点 | `sp_quality="L1_TIMEOUT"` | 超时节点 `*_L2` = `None` |
| S1 L2 能量不可用 | S3.5 重新计算产物 | `sp_quality="L2_RECALC"` | 所有 `*_L2` 可用 |

---

## 8. 工作量总结与优先级排序

### 8.1 总工作量 (Session 级别)

| 模块 | P0 任务 | P1 任务 | P2 任务 | 总计 | Session 数 |
|---|:---:|:---:|:---:|:---:|:---:|
| utils/qc_interface.py | 7 sessions (14h) | 2 sessions (4h) | - | **18h** | 9 |
| step1_anchor/ | 3 sessions (6h) | 1 session (1h) | - | **7h** | 4 |
| step3_5_sp/ (新建) | 3 sessions (6h) | 1 session (2h) | - | **8h** | 4 |
| step4_features/ | 3 sessions (6h) | 2 sessions (4h) | 1 session (3h) | **13h** | 6 |
| orchestrator.py | 3 sessions (6h) | 2 sessions (4h) | 1 session (1h) | **11h** | 6 |
| config/defaults.yaml | 1 session (1h) | - | - | **1h** | 1 |
| **合计** | **20 sessions** | **7 sessions** | **2 sessions** | **58h** | **29** |

### 8.2 实施顺序建议 (2小时 Session 单位)

```
Week 1 (P0 核心 - 20 sessions):
├── Day 1-2: utils/qc_interface.py (7 sessions)
│   ├── Session 1: _generate_input() + 测试
│   ├── Session 2: _parse_output() + 测试
│   ├── Session 3: _find_orca_binary() + _run_orca() 骨架
│   ├── Session 4: single_point() 整合 + Mock 测试
│   ├── Session 5: SMD 溶剂模型 + 测试
│   ├── Session 6: 双杂化泛函自动匹配 + 测试
│   └── Session 7: 真实 ORCA 集成测试
├── Day 3: step3_5_sp/ 目录创建 (3 sessions)
│   ├── Session 1: SPMatrixReport 数据结构
│   ├── Session 2: SPMatrixBuilder 核心逻辑
│   └── Session 3: 复用 S1 能量逻辑
├── Day 4: step1_anchor/ L2 SP 集成 (3 sessions)
│   ├── Session 1: _run_l2_sp() 方法
│   ├── Session 2: 修改 run() 返回值
│   └── Session 3: 补全 DFT 优化分支
└── Day 5: orchestrator.py S3.5 编排 (4 sessions)
    ├── Session 1: S3.5 引擎属性
    ├── Session 2-3: run_pipeline() 中插入 S3.5 调用
    └── Session 4: 修改 S4 调用以传递 sp_matrix_report

Week 2 (P0 完成 + P1/P2 - 9 sessions):
├── Day 1-2: step4_features/ 重构 (5 sessions)
│   ├── Session 1: electronic.py (HBDE)
│   ├── Session 2-3: steric_geometry.py (畸变能)
│   ├── Session 4: feature_miner.py 修改
│   └── Session 5: 集成测试
├── Day 3: orchestrator.py batch + checkpoint (3 sessions)
│   ├── Session 5-6: run_batch() 方法
│   └── Session 7: checkpoint 机制
└── Day 4-5: 集成测试 + 文档 (1 session)
    └── Session 1: 端到端集成测试
```

---

## 9. 已识别的潜在风险与缓解措施

| 风险 | 影响 | 概率 | 缓解措施 |
|:---|:---:|:---:|:---|
| ORCA 环境不可用 | S3.5 完全无法运行 | 高 | 提供 Gaussian fallback 选项；完善 Mock 测试 |
| 双杂化泛函计算时间过长 | 批量任务延迟 | 中 | 限制为 M06-2X (非双杂化) 作为默认；提供配置开关 |
| fragment_indices 未从 S2 传递 | 碎片切分失败 | 中 | 基于 forming_bonds 自动推断；提供手动覆盖选项 |
| L2 能量与 L1 几何不一致 | 物理结果存疑 | 低 | 验证时检查能量/几何相关性；添加警告日志 |
| 内存不足 (ORCA 大分子) | 计算崩溃 | 低 | 提供内存限制配置；支持 chunk 计算 |
| 配置文件错误 | 运行时失败 | 中 | 添加配置验证脚本；提供默认值 |

---

## 10. 附录: 关键数据结构变更

### 10.1 PipelineResult 扩展

```python
@dataclass
class PipelineResult:
    # ... 现有字段 ...
    
    # [NEW] S1 L2 能量
    e_product_l2: Optional[float] = None
    
    # [NEW] S3.5 SP 矩阵报告
    sp_matrix_report: Optional['SPMatrixReport'] = None
    
    # [NEW] 能量精度级别 (L1/L2)
    sp_quality: Optional[str] = None
```

### 10.2 S2 返回值扩展

```python
# retro_scanner.py

def run(self, product_xyz: Path, output_dir: Path) -> Tuple[Path, Path, Tuple]:
    """
    Returns:
        (ts_guess_xyz, reactant_xyz, forming_bonds)
        
        forming_bonds: ((atom_i, atom_j), (atom_k, atom_l))
    """
```

---

## 12. CHANGELOG 模板

每个 Session 完成后，请按以下格式更新 CHANGELOG.md：

```markdown
### Session 2026-01-10 #1: ORCAInterface._generate_input()
- [x] 创建 orca_interface.py
- [x] 实现 _generate_input() 方法
- [x] 实现 _is_double_hybrid() 辅助方法
- [x] 测试通过: `pytest tests/test_orca_interface.py::test_generate_input_m062x -v`
- [x] 生成的 .inp 文件包含 `! M062X def2-TZVPP def2/J RIJCOSX`

**遇到的问题**:
- 无

**解决方案**:
- 无

**下一步**: Session #2 实现 _parse_output()
```

```markdown
### Session 2026-01-11 #2: ORCAInterface._parse_output()
- [x] 实现 _parse_output() 方法
- [x] 添加错误处理 (检查正常终止)
- [x] 测试通过: `pytest tests/test_orca_interface.py::test_parse_output_energy -v`
- [x] 测试通过: `pytest tests/test_orca_interface.py::test_parse_output_failed -v`
- [x] 从 fixture 输出解析出 `energy = -123.45678901`

**遇到的问题**:
- 能量正则表达式需要匹配小数点

**解决方案**:
- 使用 `r"FINAL SINGLE POINT ENERGY\s+([\-\d\.]+)"` 匹配

**下一步**: Session #3 实现 _find_orca_binary() + _run_orca() 骨架
```

### CHANGELOG 维护指南

1. **每次 Session 完成后立即更新**
2. **记录实际的测试命令和输出**
3. **记录遇到的问题和解决方案**
4. **明确下一步的 Session 编号**
5. **保持格式一致**

---

## 13. 快速参考: 常用验收命令

### 开发阶段命令

```bash
# 1. ORCAInterface 单元测试
pytest tests/test_orca_interface.py -v

# 2. ORCAInterface 端到端测试
pytest tests/test_orca_interface.py::test_real_sp -v -m requires_orca

# 3. SPMatrixBuilder 测试
pytest tests/test_sp_matrix.py::test_sp_matrix_builder_run -v

# 4. FeatureMiner L2 特征测试
pytest tests/test_feature_miner.py::test_extract_features_with_sp_matrix -v

# 5. Orchestrator 端到端测试
pytest tests/test_orchestrator.py::test_pipeline_with_s35 -v

# 6. 批量处理测试
pytest tests/test_orchestrator.py::test_run_batch -v

# 7. 配置验证
python -c "from rph_core.utils.config_validator import validate_config; import yaml; config = yaml.safe_load(open('config/defaults.yaml')); print(validate_config(config))"
```

### 集成测试命令

```bash
# 完整流程测试 (需要 ORCA)
pytest tests/integration/test_full_pipeline.py::test_full_workflow_with_orca -v -s

# 降级测试 (不需要 ORCA)
pytest tests/integration/test_full_pipeline.py::test_degradation_to_l1 -v -s

# 批量处理测试
pytest tests/integration/test_batch.py::test_batch_with_mixed_success -v -s
```

---

**文档结束**

> 本计划基于 PROMOTE.md v2.1 与现有代码的详细对比分析生成。
> 所有修改均遵循串行阻塞架构原则，确保数据流的严格依赖。
> 
> **关键改进**:
> - ✅ 添加了第一个 4 小时里程碑
> - ✅ 明确了任务依赖关系图
> - ✅ 提供了可执行的验收命令
> - ✅ 定义了完整的降级策略
> - ✅ 简化了配置文件设计
> - ✅ 将任务拆分为 29 个 2 小时 Session
> - ✅ 为每个 Session 提供了验收标准
