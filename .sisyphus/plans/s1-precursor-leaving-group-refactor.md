# S1 Precursor + Leaving Group Conformer Search Refactor Plan

**版本**: v6.2.0  
**创建时间**: 2026-01-31  
**计划类型**: Major Architecture Refactor  

---

## TL;DR

> **Quick Summary**: 将 S1 从单一产物锚定扩展为完整的分子系综生成阶段，支持前体（precursor）、产物（product）、离去基团（leaving group）三类分子的构象搜索。实现小分子缓存机制，避免重复计算。
> 
> **Deliverables**: 
> - 重命名 `S1_Product` → `S1_ConfGeneration`
> - 新增前体和离去基团的构象搜索
> - 小分子全局缓存目录 `SmallMolecules/`
> - 修改 orchestrator 支持多分子输入
> 
> **Estimated Effort**: Large (5-7 days)  
> **Parallel Execution**: YES - 4 waves  
> **Critical Path**: Task 1 → Task 5 → Task 9 → Task 12

---

## Context

### Original Request
用户报告程序当前只处理产物（product），但实际需求是：
1. 对 CSV 中的**前体（precursor）**和**离去小分子（leaving group）**都进行计算
2. 小分子应在输出目录中单独设置一个全局缓存文件夹
3. 在每个反应计算前首先判定是否存在对应的小分子，如果存在则跳过计算
4. 将 S1 定位修改为 `ConfGeneration`（构象生成），对前体和产物都进行构象搜索

### CSV Data Structure (from `reaxys_cleaned.csv`)
```csv
rx_id,precursor_smiles,product_smiles_main,ylide_leaving_group,leaving_group,...
39717847,"C=CC(=O)CCCC1=CC(=O)COC1OC(C)=O","O=C1C=C2CCCC(=O)[C@@H]3C[C@H]1O[C@H]23","AcOH","AcOH",...
```

**Key Columns**:
- `rx_id`: Reaction ID
- `precursor_smiles`: 前体 SMILES
- `product_smiles_main`: 产物 SMILES
- `ylide_leaving_group`: 离去基团（优先级高）
- `leaving_group`: 离去基团（备选）

### Current S1 Architecture (v6.1)

**Current Flow**:
```
ReactionProfileHunter.run()
  └─> s1_engine.run({"product": product_smiles})
      └─> ConformerEngine.run(smiles)
          └─> S1_Product/product_min.xyz
```

**Current Limitations**:
- ❌ Only processes **product** SMILES
- ❌ Directory name `S1_Product/` implies product-only
- ❌ No handling of precursor or leaving groups
- ❌ No caching mechanism for small molecules

### Research Findings

#### 1. S1 Architecture (from bg_322cea79)
- **`AnchorPhase`** (`rph_core/steps/anchor/handler.py`): 已支持多分子输入 `Dict[str, str]`
- **`ConformerEngine`** (`rph_core/steps/conformer_search/engine.py`): 可复用处理任意 SMILES
- **Current Output**: `S1_Product/[molecule_name]/`
  - `xtb2/`: CREST artifacts
  - `cluster/`: isostat clustering
  - `dft/`: DFT OPT logs
  - `product_min.xyz`: Final global minimum

#### 2. Conformer Search Reusability (from bg_c236ecc5)
- **Workflow**: RDKit Embed → CREST → Isostat Cluster → DFT OPT-SP → Boltzmann Weighting
- **Caching**: ConformerEngine checks for `global_min.xyz` before re-running
- **Key Functions**:
  - `ConformerEngine.run(smiles)`: Main entry point
  - `CRESTInterface.run_conformer_search()`: CREST wrapper
  - `run_isostat()`: Clustering

#### 3. Small Molecule Detection
**Strategy**: Use SMILES canonicalization + molecular formula as cache key
- Heavy atoms < 10 → classify as "small molecule"
- Example: `AcOH` (acetic acid) → `C2H4O2` → cache as `SmallMolecules/C2H4O2_AcOH/`

---

## Work Objectives

### Core Objective
将 S1 步骤从**产物锚定**重构为**完整构象生成阶段**，支持前体、产物、离去基团的统一处理，并实现小分子全局缓存机制。

### Concrete Deliverables
1. **目录重命名**: `S1_Product/` → `S1_ConfGeneration/`
2. **多分子支持**: 
   - `S1_ConfGeneration/product/product_min.xyz`
   - `S1_ConfGeneration/precursor/precursor_min.xyz`
   - `SmallMolecules/C2H4O2_AcOH/molecule_min.xyz` (全局缓存)
3. **Orchestrator 修改**: 从 CSV 读取 precursor 和 leaving group
4. **S2/S3/S4 适配**: 更新路径引用

### Definition of Done
- [ ] S1 可接受 precursor/product/leaving_group 三类 SMILES
- [ ] 小分子缓存机制工作（检测存在 → 跳过计算）
- [ ] 所有步骤的路径引用更新为 `S1_ConfGeneration`
- [ ] 集成测试通过（单个反应完整流程）
- [ ] 文档更新（README + AGENTS.md）

### Must Have
- 前体构象搜索（复用 ConformerEngine）
- 离去基团处理（SMILES → 构象搜索）
- 小分子缓存逻辑（基于分子式 + SMILES）
- 目录结构重命名（S1_Product → S1_ConfGeneration）

### Must NOT Have (Guardrails)
- ❌ 不要在 S1 中实现反应中心检测（属于 S2）
- ❌ 不要改变 ConformerEngine 的核心逻辑（仅复用）
- ❌ 不要硬编码离去基团白名单（应该通用处理任意 SMILES）
- ❌ 不要在单个反应目录中重复缓存小分子

---

## Verification Strategy

### Test Infrastructure Assessment
- **Infrastructure exists**: YES (pytest in `tests/`)
- **User wants tests**: YES (after implementation)
- **Framework**: pytest
- **QA approach**: Automated tests + Manual verification

### Test Strategy
每个 TODO 包含自动化测试用例验证：
1. **Unit Tests**: 测试小分子检测、缓存查找逻辑
2. **Integration Tests**: 测试完整反应流程（precursor + product + leaving group）
3. **Regression Tests**: 确保现有产物锚定功能不受影响

---

## Execution Strategy

### Parallel Execution Waves

```
Wave 1 (Start Immediately):
├── Task 1: 实现小分子检测工具函数
├── Task 2: 实现缓存管理器
└── Task 3: 扫描代码中所有 S1_Product 引用

Wave 2 (After Wave 1):
├── Task 4: 创建 SmallMoleculeCache 类
├── Task 5: 修改 orchestrator 读取 precursor 和 leaving group
└── Task 6: 修改 AnchorPhase 支持多分子类型

Wave 3 (After Wave 2):
├── Task 7: 重命名所有 S1_Product → S1_ConfGeneration
├── Task 8: 更新 S2 的路径引用
└── Task 9: 更新 S3 的路径引用

Wave 4 (After Wave 3):
├── Task 10: 更新 S4 的路径引用
├── Task 11: 更新测试文件
└── Task 12: 集成测试 + 文档更新

Critical Path: Task 1 → Task 5 → Task 9 → Task 12
Parallel Speedup: ~50% faster than sequential
```

### Dependency Matrix

| Task | Depends On | Blocks | Can Parallelize With |
|------|------------|--------|---------------------|
| 1 | None | 2, 4 | 3 |
| 2 | 1 | 4, 5 | 3 |
| 3 | None | 7 | 1, 2 |
| 4 | 1, 2 | 5, 6 | None |
| 5 | 2, 4 | 6, 7 | None |
| 6 | 4, 5 | 7 | None |
| 7 | 3, 6 | 8, 9, 10 | None |
| 8 | 7 | 12 | 9, 10 |
| 9 | 7 | 12 | 8, 10 |
| 10 | 7 | 12 | 8, 9 |
| 11 | 7 | 12 | 8, 9, 10 |
| 12 | 8, 9, 10, 11 | None | None (final) |

---

## TODOs

### Wave 1: 基础工具和调研

- [x] 1. 实现小分子检测工具函数

  **What to do**:
  - 在 `rph_core/utils/molecule_utils.py` 中创建函数
  - `is_small_molecule(smiles: str) -> bool`: 判断是否为小分子（< 10 重原子）
  - `get_molecule_key(smiles: str) -> str`: 生成缓存 key（分子式 + 规范化 SMILES）
  - `canonicalize_smiles(smiles: str) -> str`: 使用 RDKit 规范化 SMILES
  
  **Must NOT do**:
  - 不要硬编码小分子白名单（应基于原子数动态判断）
  
  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []
  - **Reason**: 简单工具函数实现，无需特殊技能
  
  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 2, 3)
  - **Blocks**: Tasks 2, 4
  - **Blocked By**: None (can start immediately)
  
  **References**:
  - `rph_core/steps/conformer_search/engine.py`: RDKit usage patterns
  - RDKit Documentation: `Chem.MolFromSmiles`, `Chem.MolToSmiles`, `Chem.Descriptors.HeavyAtomCount`
  
  **Acceptance Criteria**:
  ```bash
  # Unit test
  python -m pytest tests/test_molecule_utils.py::test_is_small_molecule -v
  # Assert: AcOH → True, full precursor → False
  
  # Function test
  python -c "from rph_core.utils.molecule_utils import is_small_molecule, get_molecule_key; \
    assert is_small_molecule('CC(=O)O') == True; \
    assert get_molecule_key('CC(=O)O') == 'C2H4O2_CC(=O)O'"
  ```
  
  **Commit**: YES
  - Message: `feat(utils): add small molecule detection and cache key generation`
  - Files: `rph_core/utils/molecule_utils.py`, `tests/test_molecule_utils.py`
  - Pre-commit: `pytest tests/test_molecule_utils.py`

- [x] 2. 实现小分子缓存管理器

  **What to do**:
  - 在 `rph_core/utils/small_molecule_cache.py` 中创建 `SmallMoleculeCache` 类
  - 方法：
    - `__init__(cache_root: Path)`: 初始化缓存根目录 `SmallMolecules/`
    - `get_or_create(smiles: str, name: str) -> Path`: 查找缓存或创建新目录
    - `exists(smiles: str) -> bool`: 检查缓存是否存在
    - `get_path(smiles: str) -> Optional[Path]`: 获取缓存路径
  - 缓存结构: `SmallMolecules/{molecular_formula}_{name}/molecule_min.xyz`
  
  **Must NOT do**:
  - 不要在缓存类中直接调用 ConformerEngine（只管理路径）
  
  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []
  - **Reason**: 简单的文件路径管理类
  
  **Parallelization**:
  - **Can Run In Parallel**: YES (after Task 1)
  - **Parallel Group**: Wave 1 (with Tasks 1, 3)
  - **Blocks**: Tasks 4, 5
  - **Blocked By**: Task 1 (depends on `molecule_utils`)
  
  **References**:
  - Task 1: `molecule_utils.get_molecule_key()`
  - `rph_core/utils/checkpoint_manager.py`: 目录管理模式
  
  **Acceptance Criteria**:
  ```bash
  # Unit test
  python -m pytest tests/test_small_molecule_cache.py -v
  
  # Integration test
  python -c "from rph_core.utils.small_molecule_cache import SmallMoleculeCache; \
    from pathlib import Path; \
    cache = SmallMoleculeCache(Path('./test_cache')); \
    path = cache.get_or_create('CC(=O)O', 'AcOH'); \
    assert path.exists(); \
    assert cache.exists('CC(=O)O') == True"
  ```
  
  **Commit**: YES
  - Message: `feat(utils): add SmallMoleculeCache for global small molecule caching`
  - Files: `rph_core/utils/small_molecule_cache.py`, `tests/test_small_molecule_cache.py`
  - Pre-commit: `pytest tests/test_small_molecule_cache.py`

- [x] 3. 扫描代码中所有 S1_Product 引用

  **What to do**:
  - 使用 `grep` 或 `ast-grep` 查找所有包含 `S1_Product` 的文件
  - 记录到临时文件 `.sisyphus/s1_product_references.txt`
  - 分类：
    - **目录名**: `S1_Product/` 字面量
    - **变量名**: `s1_product_dir`, `S1_DIR_ALIASES`
    - **文档**: README, AGENTS.md, docstrings
  
  **Must NOT do**:
  - 不要在此任务中修改代码（仅扫描和记录）
  
  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []
  - **Reason**: 简单的代码搜索任务
  
  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1, 2)
  - **Blocks**: Task 7
  - **Blocked By**: None (can start immediately)
  
  **References**:
  - Existing grep results (already executed in this session)
  
  **Acceptance Criteria**:
  ```bash
  # Grep command executed
  grep -r "S1_Product" --include="*.py" --include="*.md" rph_core/ tests/ > .sisyphus/s1_product_references.txt
  
  # Verify output
  cat .sisyphus/s1_product_references.txt | wc -l
  # Assert: > 50 references found
  ```
  
  **Commit**: NO (documentation only)

### Wave 2: 核心逻辑修改

- [x] 4. 在 AnchorPhase 中集成 SmallMoleculeCache

  **What to do**:
  - 修改 `rph_core/steps/anchor/handler.py`
  - 在 `AnchorPhase.__init__()` 中初始化 `self.small_mol_cache = SmallMoleculeCache(work_dir.parent / "SmallMolecules")`
  - 在 `run()` 方法中，对每个分子检查：
    ```python
    if is_small_molecule(smiles):
        cache_path = self.small_mol_cache.get_or_create(smiles, name)
        if (cache_path / "molecule_min.xyz").exists():
            logger.info(f"✓ Small molecule {name} found in cache, skipping")
            # Load cached result
            return cached_result
    # Otherwise, run ConformerEngine as usual
    ```
  
  **Must NOT do**:
  - 不要改变非小分子的处理逻辑
  - 不要在 AnchorPhase 中硬编码小分子阈值（从 Task 1 获取）
  
  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []
  - **Reason**: 简单的集成逻辑，调用已有工具
  
  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Sequential (Wave 2 start)
  - **Blocks**: Tasks 5, 6
  - **Blocked By**: Tasks 1, 2
  
  **References**:
  - Task 1: `is_small_molecule()`
  - Task 2: `SmallMoleculeCache`
  - `rph_core/steps/anchor/handler.py:100-120`: Current `run()` loop
  
  **Acceptance Criteria**:
  ```bash
  # Unit test (mock)
  python -m pytest tests/test_anchor_cache_integration.py -v
  
  # Integration test
  # Run AnchorPhase twice with same small molecule
  # Assert: Second run skips computation
  ```
  
  **Commit**: YES
  - Message: `feat(anchor): integrate SmallMoleculeCache for small molecule reuse`
  - Files: `rph_core/steps/anchor/handler.py`, `tests/test_anchor_cache_integration.py`
  - Pre-commit: `pytest tests/test_anchor_cache_integration.py`

- [ ] 5. 修改 orchestrator 读取 precursor 和 leaving group

  **What to do**:
  - 修改 `rph_core/orchestrator.py` 中的 `run_single_reaction()` 方法
  - 从 CSV 读取：
    ```python
    precursor_smiles = reaction_data.get('precursor_smiles')
    leaving_group = reaction_data.get('ylide_leaving_group') or reaction_data.get('leaving_group')
    ```
  - 构建分子字典：
    ```python
    molecules = {
        "product": product_smiles,
        "precursor": precursor_smiles,
        "leaving_group": leaving_group if leaving_group and leaving_group != "" else None
    }
    # Filter None values
    molecules = {k: v for k, v in molecules.items() if v}
    ```
  - 传递给 S1: `s1_result = self.s1_engine.run(molecules)`
  
  **Must NOT do**:
  - 不要改变现有的 product-only 分支（保留向后兼容）
  - 不要在 orchestrator 中实现构象搜索逻辑
  
  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []
  - **Reason**: 简单的数据提取和传递
  
  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Sequential (after Task 4)
  - **Blocks**: Tasks 6, 7
  - **Blocked By**: Tasks 2, 4
  
  **References**:
  - CSV structure (from bg_d3a44ced)
  - `rph_core/utils/task_builder.py`: CSV parsing logic
  - `rph_core/orchestrator.py:280-300`: Current S1 invocation
  
  **Acceptance Criteria**:
  ```bash
  # Unit test
  python -m pytest tests/test_orchestrator_multi_molecule.py -v
  
  # Functional test
  # Run with CSV row containing precursor + leaving group
  # Assert: S1 receives all 3 SMILES
  # Assert: orchestrator.logger shows "Processing 3 molecules"
  ```
  
  **Commit**: YES
  - Message: `feat(orchestrator): read and pass precursor + leaving group to S1`
  - Files: `rph_core/orchestrator.py`, `tests/test_orchestrator_multi_molecule.py`
  - Pre-commit: `pytest tests/test_orchestrator_multi_molecule.py`

- [ ] 6. 修改 AnchorPhase 输出结构支持多分子

  **What to do**:
  - 修改 `AnchorPhaseResult` dataclass，确保 `anchored_molecules` 可以包含多个条目
  - 修改 `AnchorPhase.run()` 返回值，包含：
    ```python
    {
        "product": {"xyz": Path, "e_sp": float, ...},
        "precursor": {"xyz": Path, "e_sp": float, ...},
        "leaving_group": {"xyz": Path, "e_sp": float, ...}
    }
    ```
  - 在 `ConformerEngine` 内部，每个分子创建子目录：
    - `S1_ConfGeneration/product/`
    - `S1_ConfGeneration/precursor/`
    - 小分子仍使用全局缓存 `SmallMolecules/`
  
  **Must NOT do**:
  - 不要改变 ConformerEngine 的核心搜索逻辑
  - 不要在 AnchorPhase 中硬编码分子类型（应动态处理）
  
  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []
  - **Reason**: 数据结构调整，逻辑简单
  
  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Sequential (after Task 5)
  - **Blocks**: Task 7
  - **Blocked By**: Tasks 4, 5
  
  **References**:
  - `rph_core/steps/anchor/handler.py:27-33`: `AnchorPhaseResult` dataclass
  - `rph_core/steps/conformer_search/engine.py`: Directory structure creation
  
  **Acceptance Criteria**:
  ```bash
  # Integration test
  python -m pytest tests/test_anchor_multi_molecule_output.py -v
  
  # Structure verification
  # Run AnchorPhase with {"product": "...", "precursor": "..."}
  # Assert: S1_ConfGeneration/product/product_min.xyz exists
  # Assert: S1_ConfGeneration/precursor/precursor_min.xyz exists
  ```
  
  **Commit**: YES (groups with Task 7)
  - Will commit together with directory rename

### Wave 3: 目录重命名和路径更新

- [ ] 7. 重命名所有 S1_Product → S1_ConfGeneration

  **What to do**:
  - 使用 Task 3 的扫描结果，逐个文件替换：
    - **Orchestrator**: `s1_work_dir = work_dir / "S1_ConfGeneration"`
    - **Constants**: `S1_DIR_ALIASES = ["S1_ConfGeneration", "S1_Product"]` (保留向后兼容)
    - **Tests**: 更新所有测试中的路径字面量
    - **Docs**: 更新 README.md 和 AGENTS.md
  - 使用 `ast-grep` 确保替换完整
  
  **Must NOT do**:
  - 不要在单次提交中混入逻辑修改（纯重命名）
  - 不要删除 `S1_Product` 的向后兼容支持（保留在 aliases 中）
  
  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []
  - **Reason**: 批量文本替换，机械性操作
  
  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Sequential (critical path)
  - **Blocks**: Tasks 8, 9, 10, 11
  - **Blocked By**: Tasks 3, 6
  
  **References**:
  - `.sisyphus/s1_product_references.txt` (from Task 3)
  - `rph_core/steps/step4_features/mech_packager.py:37`: `S1_DIR_ALIASES` pattern
  
  **Acceptance Criteria**:
  ```bash
  # Verification
  grep -r "S1_Product" --include="*.py" rph_core/ | grep -v "S1_DIR_ALIASES" | grep -v "# backward compat"
  # Assert: No results (all replaced or in compat layer)
  
  # Regression test
  python -m pytest tests/ -k "s1" -v
  # Assert: All S1-related tests pass
  ```
  
  **Commit**: YES
  - Message: `refactor(s1): rename S1_Product to S1_ConfGeneration (with compat)`
  - Files: `rph_core/**/*.py`, `tests/**/*.py`, `README.md`, `rph_core/steps/*/AGENTS.md`
  - Pre-commit: `pytest tests/ -k "s1"`

- [ ] 8. 更新 S2 的路径引用

  **What to do**:
  - 修改 `rph_core/steps/step2_retro/retro_scanner.py`
  - 更新注释和日志中的 `S1_Product` 引用
  - 确保 S2 可以接受新的输入路径：
    - Product: `S1_ConfGeneration/product/product_min.xyz`
    - Precursor: `S1_ConfGeneration/precursor/precursor_min.xyz` (用于备选逻辑)
  
  **Must NOT do**:
  - 不要改变 S2 的核心算法（仅路径适配）
  - 不要在 S2 中添加离去基团处理（S2 只关心 product → TS guess）
  
  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []
  - **Reason**: 简单的路径更新，无算法改动
  
  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 3 (with Tasks 9, 10, 11)
  - **Blocks**: Task 12
  - **Blocked By**: Task 7
  
  **References**:
  - `rph_core/steps/step2_retro/retro_scanner.py:86-112`: Product path loading
  - S1 output contract (Task 6)
  
  **Acceptance Criteria**:
  ```bash
  # Unit test
  python -m pytest tests/test_step2_path_compat.py -v
  
  # Integration test
  # Run S2 with new S1_ConfGeneration structure
  # Assert: S2 finds product_min.xyz correctly
  ```
  
  **Commit**: YES
  - Message: `fix(s2): update path references for S1_ConfGeneration`
  - Files: `rph_core/steps/step2_retro/*.py`, `tests/test_step2_path_compat.py`
  - Pre-commit: `pytest tests/test_step2_path_compat.py`

- [ ] 9. 更新 S3 的路径引用

  **What to do**:
  - 修改 `rph_core/steps/step3_opt/ts_optimizer.py`
  - 更新日志和注释中的 `S1_Product` 引用
  - 确保 S3 可以找到 product 能量用于能垒计算
  
  **Must NOT do**:
  - 不要改变 S3 的 TS 优化逻辑
  - 不要在 S3 中添加 precursor 相关逻辑（S3 只关心 TS）
  
  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []
  - **Reason**: 简单的路径更新
  
  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 3 (with Tasks 8, 10, 11)
  - **Blocks**: Task 12
  - **Blocked By**: Task 7
  
  **References**:
  - `rph_core/steps/step3_opt/ts_optimizer.py`: Product energy usage
  
  **Acceptance Criteria**:
  ```bash
  # Unit test
  python -m pytest tests/test_step3_path_compat.py -v
  
  # Verify no hardcoded S1_Product paths remain
  grep "S1_Product" rph_core/steps/step3_opt/*.py
  # Assert: No results (except comments)
  ```
  
  **Commit**: YES
  - Message: `fix(s3): update path references for S1_ConfGeneration`
  - Files: `rph_core/steps/step3_opt/*.py`, `tests/test_step3_path_compat.py`
  - Pre-commit: `pytest tests/test_step3_path_compat.py`

- [ ] 10. 更新 S4 的路径引用

  **What to do**:
  - 修改 `rph_core/steps/step4_features/mech_packager.py`
  - 更新 `S1_DIR_ALIASES = ["S1_ConfGeneration", "S1_Product"]`
  - 修改 `_resolve_product_asset()` 查找路径：
    - 优先: `S1_ConfGeneration/product/product_min.xyz`
    - 备选: `S1_Product/product_min.xyz` (向后兼容)
  - 添加 `_resolve_precursor_asset()` 用于未来扩展
  
  **Must NOT do**:
  - 不要改变 S4 的特征提取逻辑
  - 不要在 S4 中直接处理 precursor（当前不需要）
  
  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []
  - **Reason**: 路径适配，逻辑简单
  
  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 3 (with Tasks 8, 9, 11)
  - **Blocks**: Task 12
  - **Blocked By**: Task 7
  
  **References**:
  - `rph_core/steps/step4_features/mech_packager.py:37`: Current `S1_DIR_ALIASES`
  - `rph_core/steps/step4_features/mech_packager.py:290-310`: `_resolve_product_asset()`
  
  **Acceptance Criteria**:
  ```bash
  # Unit test
  python -m pytest tests/test_s4_path_compat.py -v
  
  # Integration test
  # Run S4 with S1_ConfGeneration structure
  # Assert: S4 finds product artifacts correctly
  # Assert: features_raw.csv generated successfully
  ```
  
  **Commit**: YES
  - Message: `fix(s4): update S1 path aliases for S1_ConfGeneration`
  - Files: `rph_core/steps/step4_features/mech_packager.py`, `tests/test_s4_path_compat.py`
  - Pre-commit: `pytest tests/test_s4_path_compat.py`

- [ ] 11. 更新所有测试文件的路径

  **What to do**:
  - 批量替换 `tests/` 中的 `S1_Product` → `S1_ConfGeneration`
  - 更新 mock 数据路径
  - 确保所有 fixture 使用新的目录名
  
  **Must NOT do**:
  - 不要改变测试逻辑（仅路径更新）
  
  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []
  - **Reason**: 批量文本替换
  
  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 3 (with Tasks 8, 9, 10)
  - **Blocks**: Task 12
  - **Blocked By**: Task 7
  
  **References**:
  - Task 3 scan results
  
  **Acceptance Criteria**:
  ```bash
  # Run all tests
  python -m pytest tests/ -v
  # Assert: All tests pass
  
  # Verify replacements
  grep -r "S1_Product" tests/ | grep -v "compat"
  # Assert: Minimal results (only in compat tests)
  ```
  
  **Commit**: YES
  - Message: `test: update all test fixtures for S1_ConfGeneration`
  - Files: `tests/**/*.py`
  - Pre-commit: `pytest tests/`

### Wave 4: 集成测试和文档

- [ ] 12. 端到端集成测试 + 文档更新

  **What to do**:
  - **集成测试**:
    - 创建 `tests/test_e2e_precursor_leaving_group.py`
    - 使用真实 CSV 行运行完整流程
    - 验证目录结构：
      ```
      rph_output/rx_39717847/
      ├── S1_ConfGeneration/
      │   ├── product/product_min.xyz
      │   └── precursor/precursor_min.xyz
      ├── SmallMolecules/
      │   └── C2H4O2_AcOH/molecule_min.xyz
      ├── S2_Retro/...
      ├── S3_TransitionAnalysis/...
      └── S4_Data/...
      ```
  - **文档更新**:
    - 更新 `README.md`: 说明 S1 现在处理 precursor + leaving group
    - 更新 `rph_core/steps/anchor/AGENTS.md`: 新的 IO contract
    - 更新 `config/defaults.yaml`: 添加注释说明 precursor 和 leaving group 字段
  - **性能验证**:
    - 运行两个相同 leaving group 的反应
    - 验证第二次跳过小分子计算
  
  **Must NOT do**:
  - 不要在集成测试中使用真实 QC 计算（使用 mock）
  
  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []
  - **Reason**: 整合已有功能，编写测试
  
  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Sequential (final integration)
  - **Blocks**: None (completion task)
  - **Blocked By**: Tasks 8, 9, 10, 11
  
  **References**:
  - All previous tasks
  - `tests/test_integration.py`: 现有集成测试模式
  
  **Acceptance Criteria**:
  ```bash
  # E2E test
  python -m pytest tests/test_e2e_precursor_leaving_group.py -v -s
  # Assert: All assertions pass
  # Assert: SmallMolecules/ cache hit on second run
  
  # Documentation check
  grep -i "precursor" README.md
  grep -i "leaving.group" README.md
  # Assert: Both found and explained
  
  # Full test suite
  python -m pytest tests/ -v
  # Assert: 100% pass rate
  ```
  
  **Commit**: YES
  - Message: `test(e2e): add integration test for precursor + leaving group workflow`
  - Files: `tests/test_e2e_precursor_leaving_group.py`, `README.md`, `rph_core/steps/anchor/AGENTS.md`, `config/defaults.yaml`
  - Pre-commit: `pytest tests/test_e2e_precursor_leaving_group.py`

---

## Commit Strategy

| After Task | Message | Files | Verification |
|------------|---------|-------|--------------|
| 1 | `feat(utils): add small molecule detection and cache key generation` | `rph_core/utils/molecule_utils.py`, `tests/test_molecule_utils.py` | `pytest tests/test_molecule_utils.py` |
| 2 | `feat(utils): add SmallMoleculeCache for global small molecule caching` | `rph_core/utils/small_molecule_cache.py`, `tests/test_small_molecule_cache.py` | `pytest tests/test_small_molecule_cache.py` |
| 4 | `feat(anchor): integrate SmallMoleculeCache for small molecule reuse` | `rph_core/steps/anchor/handler.py`, `tests/test_anchor_cache_integration.py` | `pytest tests/test_anchor_cache_integration.py` |
| 5 | `feat(orchestrator): read and pass precursor + leaving group to S1` | `rph_core/orchestrator.py`, `tests/test_orchestrator_multi_molecule.py` | `pytest tests/test_orchestrator_multi_molecule.py` |
| 7 | `refactor(s1): rename S1_Product to S1_ConfGeneration (with compat)` | `rph_core/**/*.py`, `tests/**/*.py`, `README.md`, `AGENTS.md` | `pytest tests/ -k "s1"` |
| 8 | `fix(s2): update path references for S1_ConfGeneration` | `rph_core/steps/step2_retro/*.py`, `tests/test_step2_path_compat.py` | `pytest tests/test_step2_path_compat.py` |
| 9 | `fix(s3): update path references for S1_ConfGeneration` | `rph_core/steps/step3_opt/*.py`, `tests/test_step3_path_compat.py` | `pytest tests/test_step3_path_compat.py` |
| 10 | `fix(s4): update S1 path aliases for S1_ConfGeneration` | `rph_core/steps/step4_features/mech_packager.py`, `tests/test_s4_path_compat.py` | `pytest tests/test_s4_path_compat.py` |
| 11 | `test: update all test fixtures for S1_ConfGeneration` | `tests/**/*.py` | `pytest tests/` |
| 12 | `test(e2e): add integration test for precursor + leaving group workflow` | `tests/test_e2e_precursor_leaving_group.py`, `README.md`, docs | `pytest tests/test_e2e_precursor_leaving_group.py` |

---

## Success Criteria

### Verification Commands
```bash
# 1. Precursor 构象搜索工作
ls -la rph_output/rx_39717847/S1_ConfGeneration/precursor/precursor_min.xyz
# Expected: file exists

# 2. 小分子缓存工作
ls -la SmallMolecules/C2H4O2_AcOH/molecule_min.xyz
# Expected: file exists

# 3. 缓存复用工作（运行两次相同反应）
# First run: "Running conformer search for leaving_group..."
# Second run: "✓ Small molecule AcOH found in cache, skipping"

# 4. 所有测试通过
pytest tests/ -v
# Expected: 100% pass rate

# 5. 目录重命名完成
find rph_output/ -name "S1_Product" -type d
# Expected: no results (or only in old runs for compat)
```

### Final Checklist
- [ ] ✅ Precursor conformer search 工作正常
- [ ] ✅ Leaving group conformer search 工作正常
- [ ] ✅ 小分子缓存机制正常（检测、跳过、复用）
- [ ] ✅ 目录从 `S1_Product` 重命名为 `S1_ConfGeneration`
- [ ] ✅ S2/S3/S4 路径引用全部更新
- [ ] ✅ 向后兼容性保留（S1_DIR_ALIASES）
- [ ] ✅ 所有单元测试通过
- [ ] ✅ 集成测试通过（完整反应流程）
- [ ] ✅ 文档更新完成（README + AGENTS.md）

---

## Notes

### 向后兼容性
- 保留 `S1_DIR_ALIASES = ["S1_ConfGeneration", "S1_Product"]`，确保旧的 S1_Product 目录仍可被 S4 读取
- 小分子缓存是**新增功能**，不影响现有产物锚定流程

### 性能优化
- 小分子缓存预期减少 **30-50%** 的重复计算（对于共用离去基团的反应）
- 前体构象搜索增加计算时间，但对反应机理理解至关重要

### 未来扩展
- [ ] S4 可扩展为提取 precursor 相关特征（如前体扭曲能）
- [ ] 小分子缓存可扩展为支持自定义阈值（当前硬编码 < 10 重原子）
- [ ] 可添加 `--skip-precursor` CLI 参数用于仅产物计算模式
