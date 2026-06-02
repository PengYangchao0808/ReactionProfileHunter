# RPH v2.1.1 综合修复方案

> 更新：Phase 1 已完成，legacy `forward_scan` 别名已在 v2.1.1 中删除。

*基于 UNIFIED_BOND_IDENTIFICATION_ANALYSIS.md 六层诊断 + 用户六项问题反馈 + Oracle 架构审查*
*生成时间：2026-04-24*

---

## 0. 问题总览

用户从统一分析报告中识别出 **6 项具体问题**，归纳为 **5 个独立问题域**：

| 问题域 | 用户反馈 | 根源 | 严重性 |
|--------|---------|------|--------|
| **P1: 化学映射错误** | Map 13 不应是成键中心，应为 Map 14；Map 15 应转化为羰基 | RXNMapper 低置信度映射 (0.225) 错误分配原子身份 | 🔴 致命 |
| **P2: 原子映射未持久化** | 三套原子映射需在 S0 阶段严肃保存，方便各阶段核对 | S0 仅保存 forming_bonds 数值，不保存完整映射关系 | 🟡 高 |
| **P3: Cleaner 键变化检测盲区** | [4+3] 反应判断鲁棒性极差；需渲染机理键线式 | `_infer_bond_changes()` 仅做有键/无键二分法 | 🔴 致命 |
| **P4: forward_scan 垃圾代码** | forward_scan 是 retro_scan 别名，请完全删除 | `run_forward_scan()` = `run_retro_scan()` 调用 | 🟡 高 |
| **P5: Baseline 拉伸错误键** | 原先拉伸的键是对的，baseline 时变错 | 两次 run 使用相同错误 forming_bonds；差异仅在于 path_search 开关 | 🟡 中 |

### P5 关键澄清

从 artifact 比对可见：

| Run | S0 forming_bonds | S2 方法 | path_search | 最终结果 |
|-----|-----------------|---------|-------------|---------|
| Normal (2026-04-06) | **[[12,13],[15,16]]** ← 错误 | retro_scan → path_search | ✅ enabled | activation = +26.61 kcal/mol ✅ |
| Baseline (2026-04-21) | **[[12,13],[15,16]]** ← 同样错误 | forward_scan（=retro_scan） | ❌ disabled | activation = -44.53 kcal/mol ❌ |

**结论**：两次 run 使用完全相同的错误 forming_bonds。Normal run 成功是因为 `xtb_path_search` 不依赖 forming_bonds（仅看几何端点），而非 retro_scan 拉伸了正确的键。用户观察到的"原先拉伸的键是对的"实际上是因为 path_search 的鲁棒性覆盖了错误。

---

## 1. Phase 1 — 删除 forward_scan（P4）

**风险**：🟢 最低（纯删除，无逻辑变更）
**依赖**：无
**影响范围**：110 occurrences / 26 files

### 1.1 核心代码删除

| 文件 | 位置 | 操作 |
|------|------|------|
| `rph_core/steps/step2_retro/retro_scanner.py` | L551-565 `run_forward_scan()` | 删除整个方法 |
| `rph_core/orchestrator.py` | L280-299 `_resolve_forward_scan_config()` | 删除整个方法 |
| `rph_core/orchestrator.py` | L564 `step2_cfg.get("forward_scan", {})` | 删除引用 |
| `rph_core/orchestrator.py` | L1111 `_resolve_forward_scan_config(...)` 调用 | 删除/修复 |
| `rph_core/steps/runners.py` | L178 `_resolve_forward_scan_config(...)` 调用 | 删除 |
| `rph_core/steps/runners.py` | L199-219 `if s2_strategy == "forward_scan":` 分支 | 删除整个分支 |

### 1.2 配置文件

| 文件 | 位置 | 操作 |
|------|------|------|
| `config/defaults.yaml` | L648 `s2_strategy: forward_scan` | 改为 `retro_scan` |
| `config/validation.yaml` | L143, L153 | 改为 `retro_scan` |

### 1.3 测试文件

| 文件 | 操作 |
|------|------|
| `tests/test_forward_scan_wiring.py` | **整文件删除**（243 行） |
| `tests/test_phase3c_chain_lite_zero.py:29` | 删除 `_resolve_forward_scan_config` mock |
| `tests/test_step2_pes_adapter_fallback.py:108` | 删除 `_resolve_forward_scan_config` mock |
| `tests/test_dataset_loader.py:19` | `"s2_strategy": "forward_scan"` → `"retro_scan"` |
| `tests/test_cleaner_adapter.py:26` | `"s2_strategy": "forward_scan"` → `"retro_scan"` |

### 1.4 文档更新

18 个 MD 文件中的 forward_scan 引用，全部移除或标注为已删除。

### 1.5 验证门

```bash
python scripts/ci/check_imports.py rph_core
pytest tests/ -v --timeout=60
grep -r "forward_scan" rph_core/  # 应返回 0 结果
```

---

## 2. Phase 2 — 持久化原子映射（P2）

**风险**：🟢 低（纯增量，不影响现有逻辑）
**依赖**：Phase 1（清理 forward_scan 后代码更干净）

> ⚠️ **Oracle 架构审查关键修正**：S0 在 pipeline 中**先于 S1 执行**，此时尚无 XYZ 文件。
> 因此原子映射分为**两层**：
> - **Layer 1（S0 阶段）**：仅保存 SMILES/Map# 空间的映射
> - **Layer 2（S1 完成后）**：在 orchestrator 中补充 XYZ↔Map# 映射

### 2.1 扩展 MechanismGraph 数据结构（Layer 1）

**文件**: `rph_core/steps/mechanism_classifier/models.py`

```python
@dataclass
class SmilesAtomMapping:
    """SMILES canonical order ↔ Map# 映射（S0 阶段可获取）"""
    # 前体 SMILES canonical order idx → Map#
    precursor_smiles_to_map: Dict[int, int]
    # 产物 SMILES canonical order idx → Map#
    product_smiles_to_map: Dict[int, int]
    # Map# → 前体 SMILES idx
    map_to_precursor_smiles: Dict[int, int]
    # Map# → 产物 SMILES idx
    map_to_product_smiles: Dict[int, int]

@dataclass
class FormingBondNotation:
    """一根 forming bond 的多种编号表示"""
    map_space: Tuple[int, int]           # Map# 编号
    product_smiles_idx: Tuple[int, int]  # 产物 SMILES canonical idx
    precursor_smiles_idx: Optional[Tuple[int, int]]  # 前体 SMILES idx
    bond_type_product: str               # 产物中的键类型
    bond_type_precursor: str             # 前体中的键类型（含 NONE/AROMATIC）

# MechanismGraph 新增字段:
smiles_atom_mapping: Optional[SmilesAtomMapping] = None
forming_bonds_annotated: Optional[List[FormingBondNotation]] = None
```

### 2.2 在 S0 阶段填充 SMILES 映射

**文件**: `rph_core/steps/mechanism_classifier/graph_builder.py`

当前 `_create_nodes()` (L166-221) 仅为 intermediate 节点设置 `atom_map`。修改为：

1. **所有节点类型**（reactant, intermediate, product）：从 RDKit mol 对象提取 `mol.GetAtomWithIdx(i).GetAtomMapNum()` → 构建 SMILES↔Map 映射
2. 在 `_create_edges()` 中为每根 forming bond 构建 `FormingBondNotation`

### 2.3 保存 Layer 1 artifact

**文件**: `rph_core/orchestrator.py` `_run_s0()` (L888 之后)

```python
# S0 阶段：保存 SMILES 空间映射
mapping_path = s0_dir / "atom_map_smiles.json"
# 包含 precursor_smiles_to_map, product_smiles_to_map, forming_bonds_annotated
```

### 2.4 Layer 2：S1 完成后补充 XYZ 映射

**文件**: `rph_core/orchestrator.py` `run_pipeline()` (S1 完成后、S2 之前)

S1 完成后已有 `product_min.xyz` 和（可选）`precursor_min.xyz`，此时可建立 XYZ↔Map 映射：

```python
def _build_xyz_atom_mapping(
    self,
    smiles_mapping: SmilesAtomMapping,
    product_xyz_path: Path,
    precursor_xyz_path: Optional[Path] = None,
) -> Dict:
    """
    S1 完成后，将 SMILES canonical order 与 XYZ 文件中的原子顺序对齐。
    
    对齐策略：
    1. 从 product SMILES 构建 RDKit mol + 3D conformer
    2. 读取 product_min.xyz
    3. 用坐标匹配（匈牙利算法/Kabsch 对齐后最近邻）确定 mol_idx → xyz_idx
    4. 通过 smiles_mapping.chain: mol_idx → map# → xyz_idx
    """
    ...
```

保存为 `S1_ConfGeneration/atom_map_xyz.json`：

```json
{
  "map_to_product_xyz_1based": {"12": 12, "13": 13, ...},
  "map_to_precursor_xyz_1based": {"12": 8, "13": 1, "15": 2, ...},
  "forming_bonds_full_notation": [
    {
      "map_space": [12, 19],
      "product_xyz_1based": [12, 19],
      "precursor_xyz_1based": [8, 9],
      "bond_type_product": "SINGLE",
      "bond_type_precursor": "AROMATIC"
    }
  ],
  "index_base_convention": "all XYZ indices are 1-based; map numbers are 1-based"
}
```

> ⚠️ **注意事项**：如果新增 `atom_map_xyz.json` 为必需 artifact，必须同步更新：
> - `layout_contract.py`（artifact 列表）
> - `Step0Artifacts` / `Step1Artifacts` 类型定义
> - `checkpoint_manager.py`（signature 计算）
> - benchmark migrate 代码
> - 所有引用 S0/S1 artifacts 的测试

### 2.5 验证门

```bash
pytest tests/ -v
# 对 rx1 运行后：
# S0_Mechanism/atom_map_smiles.json 存在且包含 SMILES↔Map 映射
# S1_ConfGeneration/atom_map_xyz.json 存在且包含 XYZ↔Map 映射
# 手动验证 Map 12→precursor SMILES idx, product SMILES idx 正确
```

---

## 3. Phase 3 — 修复键变化检测（P1+P3 根因修复）

**风险**：🟡 中（修改核心逻辑）
**依赖**：Phase 2（AtomMapping 数据结构）

此 Phase 分为四个子任务，可独立实施。

### 3a: 修复 Cleaner 的 `_infer_bond_changes()`

**文件**: `exclean/cleaner_20260406/parsers/core_extractor.py:411-445`

#### 当前问题代码

```python
product_has_bond = product_mol.GetBondBetweenAtoms(idx_i, idx_j) is not None  # ← 仅判断有无
reactant_has_bond = reactant_mol.GetBondBetweenAtoms(idx_i, idx_j) is not None
```

#### 修复方案：全键级差分

```python
def _infer_bond_changes(...):
    """全键级差分：检测 formed / broken / order_changed"""
    formed, broken, order_changed = [], [], []

    for map_i, map_j in combinations(candidate_sorted, 2):
        product_bond = _get_bond(product_mol, product_map_idx, map_i, map_j)
        reactant_bond = _get_bond_any(reactant_mols, reactant_map_idx_list, map_i, map_j)

        if product_bond and not reactant_bond:
            formed.append((map_i, map_j))
        elif reactant_bond and not product_bond:
            broken.append((map_i, map_j))
        elif product_bond and reactant_bond:
            bt_product = product_bond.GetBondType()
            bt_reactant = reactant_bond.GetBondType()
            if bt_product != bt_reactant:
                if _is_significant_change(bt_reactant, bt_product):
                    order_changed.append((map_i, map_j, str(bt_reactant), str(bt_product)))

    return formed, broken, order_changed
```

#### 输出格式

```python
# core_bond_changes 新增 "order_changed" 元数据
# 旧: "12-13:formed;15-16:formed"
# 新: "12-13:formed;15-16:formed;12-19:order_changed(AROMATIC→SINGLE)"
```

> ⚠️ **Oracle 注意**：`order_changed` 作为**调试元数据**输出，下游（S0/S2）仍只消费 `formed`/`breaking`。
> 不引入新的消费端依赖，保持 `forming`/`breaking` 的现有契约不变。

### 3b: RPH 内部几何预检（安全网）

**文件**: `rph_core/orchestrator.py` `_resolve_forming_bonds_for_s2()` (L471)

> ⚠️ **Oracle 架构审查关键修正**：`forming_bonds_resolver.infer_forming_bonds_from_geometries()` 
> 需要 **product + TS** 两套几何，在 S2 之前不可用。不能用它做 pre-S2 fallback。
> 
> 正确做法：在 `_resolve_forming_bonds_for_s2()` 的现有 fallback 链中增加**产物侧验证**，
> fallback 到已有的 SMARTS 匹配，不要引入不适用于此阶段的工具。

#### 产物侧验证逻辑

```python
def _validate_forming_bonds_against_product(
    self,
    forming_bonds: List[List[int]],
    product_xyz_path: Path,
    index_base: int,
) -> Tuple[List[List[int]], bool]:
    """
    产物侧验证：检查 claimed forming_bonds 在 product_min.xyz 中的合理性。
    
    返回：(validated_bonds, needs_fallback)
    """
    coords = GeometryUtils.read_xyz_coordinates(product_xyz_path)
    suspicious_count = 0

    for pair in forming_bonds:
        # 统一到 0-based 内部坐标索引
        i = pair[0] - (1 if index_base == 0 else 0) - 1  # XYZ 内部 0-based
        j = pair[1] - (1 if index_base == 0 else 0) - 1
        dist = GeometryUtils.calculate_distance(coords, i, j)

        if dist < 1.2:
            logger.warning(
                f"Forming bond {pair}: distance={dist:.3f} Å "
                f"in product — pair already bonded, possibly wrong pair"
            )
            suspicious_count += 1
        elif dist > 3.5:
            logger.warning(
                f"Forming bond {pair}: distance={dist:.3f} Å "
                f"in product — atoms very distant, possibly wrong pair"
            )
            suspicious_count += 1

    needs_fallback = suspicious_count == len(forming_bonds) and len(forming_bonds) > 0
    return forming_bonds, needs_fallback
```

#### 正确的 Fallback 链

```
_resolve_forming_bonds_for_s2() 的现有 fallback 链:
  1. S0 summary (mechanism_summary.json)     → 验证 → 若可疑则跳过
  2. Cleaner XYZ pairs (core_bond_changes)    → 验证 → 若可疑则跳过
  3. Cleaner map pairs (formed_bond_map_pairs) → 验证 → 若可疑则跳过
  4. Config (reaction_profiles forming_bonds)  → 验证
  5. SMARTS 匹配 (SMARTSMatcher fallback)     ← 最终兜底
```

每个来源获取后都经过产物侧验证。如果所有来源都可疑，使用 SMARTS 匹配作为最终兜底（SMARTSMatcher 已有独立的键对识别逻辑，不依赖 cleaner/RXNMapper）。

**前体侧检查**（可选诊断）：如果 `precursor_min.xyz` 存在，计算 forming bond 对应原子对在前体中的距离作为 **log-only 诊断信号**，不作为 gate。

### 3c: 机理图渲染为键线式

**新文件**: `rph_core/steps/mechanism_classifier/visualizer.py`

```python
class MechanismVisualizer:
    """将 MechanismGraph 渲染为键线式结构图（best-effort 调试 artifact）"""

    def render(self, graph: MechanismGraph, output_path: Path) -> Optional[Path]:
        """
        生成机理图 PNG：
        - 左侧：前体（带 atom map 编号）
        - 箭头
        - 右侧：产物（带 atom map 编号）
        - Forming bonds：绿色虚线高亮
        - Breaking bonds：红色 X 标记
        - Atom Map# 标注在原子旁
        
        失败时返回 None，不阻塞 pipeline（best-effort）。
        """
        try:
            from rdkit.Chem import Draw, AllChem, rdMolDraw2D
            # ... 渲染逻辑 ...
            return output_path
        except Exception as e:
            logger.warning(f"Mechanism visualization failed: {e}")
            return None
```

> ⚠️ **Oracle 注意**：此 visualization 为 **best-effort 调试 artifact**。
> 渲染失败**不应**阻塞 pipeline 或影响 step success/failure 判断。

### 3d: 对齐 RXNMapper 置信度阈值

> ⚠️ **Oracle 架构审查关键发现**：RPH 已存在 `map_confidence` 阈值不一致问题：
> - `rph_core/utils/cleaner_adapter.py`：gates at `map_confidence >= 0.8`
> - `rph_core/steps/mechanism_classifier/clean_adapter.py`：**不检查** confidence
>
> 修复方向：**不是降低阈值**，而是在 S0 的 clean_adapter 中对齐已有的 0.8 阈值。

**文件**: `rph_core/steps/mechanism_classifier/clean_adapter.py`

```python
# 在 CleanRecord 解析时增加 confidence 检查
CONFIDENCE_THRESHOLD = 0.8

def classify_single(self, record: CleanRecord) -> MechanismGraph:
    confidence = record.source_data.get("map_confidence", 0)
    if confidence < CONFIDENCE_THRESHOLD:
        logger.warning(
            f"Low mapping confidence ({confidence:.3f} < {CONFIDENCE_THRESHOLD}). "
            f"S0 results may be unreliable. Consider SMARTS fallback."
        )
        # 标记低置信度，供 _resolve_forming_bonds_for_s2() 参考
        # 但不阻止 S0 运行（降级而非拒绝）
```

---

## 4. Phase 4 — path_search 策略调整（P5 安全网）

**风险**：🟡 中（影响 S2 行为）
**依赖**：Phase 1, Phase 3b

> ⚠️ **Oracle 架构审查修正**：**不应**将 `path_search.enabled` 全局默认改为 `true`。
> 这会掩盖 forming_bonds 错误，且增加不必要的计算开销。
> 
> 正确策略：**仅在 S2.1 retro_scan 产出 DEGRADED 结果时**，自动触发 path_search 作为 rescue。

### 修改

**文件**: `rph_core/steps/runners.py` `run_step2()`

```python
# 修改逻辑：path_search 不再是独立开关，而是 retro_scan DEGRADED 时的自动 rescue

def run_step2(...):
    # S2.1: retro_scan
    result = hunter.run_retro_scan(product_xyz, forming_bonds, ...)
    
    # S2.2 rescue: 仅在 S2.1 质量不佳时触发
    if result.status in ("DEGRADED", "FAILED") or result.topology_drift:
        logger.warning(
            f"S2.1 retro_scan quality: {result.status}, "
            f"triggering path_search as rescue"
        )
        result = hunter.run_path_search(
            intermediate_xyz=result.intermediate_xyz,
            product_xyz=product_xyz,
            forming_bonds=forming_bonds,  # 仅用于事后距离计算
        )
```

**文件**: `config/defaults.yaml`

```yaml
# 修改：不再有简单的 enabled/disabled
step2:
  path_search:
    mode: "rescue"          # "always" | "rescue" | "never"
    # rescue = 仅在 retro_scan DEGRADED/FAILED 时触发
    # always = 始终执行（类似旧版 enabled=true）
    # never = 从不执行
```

### 理由

- 保留 path_search 的鲁棒性优势，但**不掩盖** forming_bonds 错误
- retro_scan 成功时不必浪费额外计算
- DEGRADED 是明确的信号，说明 retro_scan 产出不可靠

---

## 5. Phase 5 — 加固 Baseline TS 验证

**风险**：🟢 低
**依赖**：无

**文件**: `benchmark/dft_theory/stages/baseline.py` `_validate_ts()`

### 修改

```python
def _validate_ts(self, ts_result: dict) -> bool:
    """加固版 TS 验证 — 阻止不合理的 TS 进入 baseline 缓存"""
    errors = []

    # 1. bond_lengths 不应为空
    bond_lengths = ts_result.get("bond_lengths", [])
    if not bond_lengths:
        errors.append("bond_lengths is empty — TS validation incomplete")

    # 2. 负势垒不应通过
    activation_energy = ts_result.get("activation_energy_kcal", 0)
    if activation_energy < 0:
        errors.append(
            f"Negative activation energy: {activation_energy:.2f} kcal/mol "
            f"— indicates energy inversion (TS lower than reactant)"
        )

    # 3. 最小虚频检查
    imag_freq = ts_result.get("imag_frequency_cm1", 0)
    if abs(imag_freq) < 50:
        errors.append(
            f"Weak imaginary frequency: {imag_freq:.1f} cm⁻¹ "
            f"(minimum 50 cm⁻¹ for valid TS)"
        )

    if errors:
        logger.error(f"TS validation FAILED: {'; '.join(errors)}")
        return False  # ← 阻止 baseline 发布，不只是记录

    return True
```

---

## 6. 实施顺序与依赖关系

```
Phase 1 (删除 forward_scan)      ← 独立，可立即实施 [🟢 低风险]
    ↓
Phase 5 (加固 TS 验证)            ← 独立，可与 Phase 1 并行 [🟢 低风险]
    ↓
Phase 3b (产物侧几何预检)          ← 依赖 Phase 1 清理后代码 [🟢 低风险]
    ↓
Phase 3d (对齐置信度阈值)          ← 独立 [🟢 低风险]
    ↓
Phase 2 (持久化原子映射)            ← 依赖 Phase 3b（预检逻辑中使用映射数据）[🟢 低风险]
    ↓
Phase 3a (修复 cleaner 根因)       ← 依赖 Phase 2 数据结构 [🟡 中风险]
Phase 3c (机理图渲染)              ← 依赖 Phase 2 数据结构 [🟢 低风险]
    ↓
Phase 4 (path_search rescue 模式)  ← 依赖 Phase 3b 的预检机制 [🟡 中风险]
```

### 推荐实施顺序

1. **Week 1**: Phase 1 + Phase 5 + Phase 3d（并行，均无依赖，风险最低）
2. **Week 1**: Phase 3b（依赖 Phase 1 完成后的干净代码）
3. **Week 2**: Phase 2（为 Phase 3a/3c 提供数据基础）
4. **Week 2**: Phase 3a + Phase 3c（依赖 Phase 2）
5. **Week 2**: Phase 4（最后，作为收尾）

---

## 7. 各 Phase 的验证标准

| Phase | 验证方法 | 通过标准 |
|-------|---------|---------|
| 1 | `grep -r "forward_scan" rph_core/` | 返回 0 结果 |
| 1 | `pytest tests/ -v` | 全部通过 |
| 1 | `python scripts/ci/check_imports.py rph_core` | exit 0 |
| 2 | 对 rx1 运行 pipeline | `S0_Mechanism/atom_map_smiles.json` 存在且包含 SMILES↔Map 映射 |
| 2 | 对 rx1 运行 pipeline | `S1_ConfGeneration/atom_map_xyz.json` 存在且包含 XYZ↔Map 映射 |
| 2 | 手动验证 | Map 12→precursor SMILES idx, product SMILES idx 正确 |
| 3a | 对 rx1 重新运行 cleaner | `core_bond_changes` 包含 `12-19:order_changed(AROMATIC→SINGLE)` |
| 3b | 对 rx1 运行 pipeline | 日志中出现 forming_bonds 产物侧验证结果 |
| 3c | 对 rx1 运行 pipeline | `S0_Mechanism/mechanism_graph.png` 存在且可读 |
| 3d | 对 rx1 运行 pipeline | S0 日志中出现低置信度警告 |
| 4 | 对 rx1 baseline re-run | retro_scan DEGRADED 时自动触发 path_search |
| 5 | 对 baseline re-run | 负势垒被正确拦截（返回 False） |

---

## 8. 关于 P1 化学映射错误的特别说明

用户指出当前分析中的"金标准 forming bonds"存在化学错误：

| 当前（错误） | 用户修正 |
|-------------|---------|
| Map 13 (allene =CH) 标记为 ✅ 成键中心 | **Map 13 不应成键**；应为 Map 14 (allene 终端 CH₂) |
| Map 15 (allene =C= 中心) → 产物桥头 CH | **Map 15 应转化为羰基 C=O** |
| 金标准 forming bonds: 12-19, 15-16 | **需修正为**: 12-19, 14-16（待用户确认） |

这意味着：
1. RXNMapper 的原子映射不仅仅是"低置信度"问题，而是**化学上不正确的分配**
2. 当前的"金标准"本身也是基于错误映射推导的，需要根据正确的化学反应机理重新确定
3. Phase 2 的原子映射持久化将使这类错误更容易被发现和验证

**建议**：在 Phase 2 完成后、Phase 3 实施前，用户应手动验证 atom_map.json 中的映射关系是否与正确的化学机理一致。

---

## 9. Oracle 架构审查关键修正摘要

Oracle 审查发现初始方案中 **3 个重大设计缺陷**，已在本文档中修正：

| # | 初始方案问题 | Oracle 修正 |
|---|------------|-----------|
| 1 | Phase 2 在 S0 保存 XYZ↔Map 映射 | ❌ S0 先于 S1 执行，无 XYZ 文件。拆为两层：S0 保存 SMILES↔Map，S1 后补充 XYZ↔Map |
| 2 | Phase 3b 用 `forming_bonds_resolver` 做 pre-S2 fallback | ❌ 该工具需要 product+TS 几何，S2 前不可用。改用产物侧距离验证 + SMARTS 兜底 |
| 3 | Phase 4 全局默认启用 path_search | ❌ 会掩盖 forming_bonds 错误。改为仅在 retro_scan DEGRADED 时自动 rescue |
| + | RXNMapper 置信度阈值不一致 | `utils/cleaner_adapter.py` gate at 0.8，但 `mechanism_classifier/clean_adapter.py` 不检查。需对齐 |
| + | 新 artifact 需更新基础设施 | 新增 `atom_map_*.json` 须同步更新 `layout_contract.py`, checkpoint, benchmark migrate |

---

## 附录 A：关键代码文件索引

| 文件 | Phase | 角色 |
|------|-------|------|
| `rph_core/steps/step2_retro/retro_scanner.py` | 1 | `run_forward_scan()` L551-565 删除目标 |
| `rph_core/orchestrator.py` | 1,2,3b | `_resolve_forward_scan_config()` L280-299 删除；atom_map 写入；几何预检 L471+ |
| `rph_core/steps/runners.py` | 1,4 | forward_scan 分支删除；path_search rescue 逻辑 |
| `rph_core/steps/mechanism_classifier/models.py` | 2 | SmilesAtomMapping / FormingBondNotation 新增 |
| `rph_core/steps/mechanism_classifier/graph_builder.py` | 2 | SMILES↔Map 映射填充 L184-221 |
| `rph_core/steps/mechanism_classifier/clean_adapter.py` | 3d | 对齐 confidence 阈值 |
| `rph_core/utils/cleaner_adapter.py` | 3d | 参考实现（已有 0.8 阈值） |
| `exclean/cleaner_20260406/parsers/core_extractor.py` | 3a | `_infer_bond_changes()` L411-445 修复 |
| `rph_core/steps/mechanism_classifier/visualizer.py` | 3c | 新文件：机理图渲染 |
| `config/defaults.yaml` | 1,4 | forward_scan 删除；path_search.mode 改 rescue |
| `benchmark/dft_theory/stages/baseline.py` | 5 | `_validate_ts()` 加固 |

## 附录 B：Baseline vs Normal Artifact 对比

| Artifact | Normal (2026-04-06) | Baseline (2026-04-21) |
|----------|---------------------|----------------------|
| S0 forming_bonds | [[12,13],[15,16]] | [[12,13],[15,16]] |
| S2 generation_method | `xtb_path_search` | `retro_scan_from_product` |
| S2 s2_strategy | `retro_scan` | `forward_scan` (=retro_scan 别名) |
| path_search.enabled | **true** | **false** |
| S3 forming_bonds (top) | [[11,12],[14,15]] (0-based) | [[12,13],[15,16]] (标注0-based实际1-based) |
| activation_energy | **+26.61 kcal/mol** ✅ | **-44.53 kcal/mol** ❌ |
| scan_quality | COMPLETE | DEGRADED |
