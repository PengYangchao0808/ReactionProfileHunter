# S1步骤日志与续算功能问题报告

## 报告日期
2026-03-13

## 版本信息
ReactionProfileHunter v6.2

---

## 1. 执行摘要

经过对S1步骤（构象搜索与优化）代码的深入分析，发现当前实现存在两个关键问题：

1. **日志粒度不足**: 无法反映每个分子各构象的具体处理状态（OPT/SP阶段细节缺失）
2. **续算功能不完善**: 无法从上次中断的具体构象/阶段继续，导致重复计算

---

## 2. 问题详细分析

### 2.1 问题一：日志无法完整反映任务状态

#### 2.1.1 当前状态记录机制

| 层级 | 当前记录方式 | 存在的问题 |
|------|-------------|-----------|
| **Pipeline级** | `pipeline.state` JSON文件 | 只记录步骤完成状态，不记录子步骤 |
| **分子级** | `AnchorPhase` 循环中更新 `.rph_step_status.json` | 仅记录分子索引，无构象详情 |
| **构象级** | `conformer_thermo.csv` | 仅记录成功的构象，失败/中断的不记录 |
| **计算级** | 日志输出到控制台/文件 | 无结构化状态记录，无法恢复 |

#### 2.1.2 关键代码位置

**文件**: `rph_core/steps/anchor/handler.py:151-161`
```python
# 当前仅记录分子级别的状态
status_file = self.base_work_dir / ".rph_step_status.json"
status_data = {
    "step": "s1",
    "molecule": name,
    "index": idx + 1,
    "total": total_mols,
    "smiles": smiles,
    "status": "running"
}
```

**问题**: 该状态文件仅记录到分子级别，无法追踪：
- CREST搜索是否完成
- 发现了多少构象
- 各构象的OPT/SP状态
- 哪些构象失败/成功

#### 2.1.3 DFT OPT-SP循环中的日志缺陷

**文件**: `rph_core/steps/conformer_search/engine.py:690-893`

当前 `_step_dft_opt_sp_coupled` 方法：
- 仅通过 `logger.info` 输出进度
- 没有持久化状态记录
- 失败后无法知道哪些构象已完成

```python
# 当前实现：仅记录成功的records
def _step_dft_opt_sp_coupled(self, candidates: List[Path]) -> Tuple[Path, float]:
    for idx, xyz_file in enumerate(candidates):
        conf_name = f"conf_{idx:03d}"
        # ... OPT尝试循环
        for attempt in range(2):
            # ... 执行Gaussian OPT
            opt_converged, next_xyz = self._run_gaussian_opt(...)
            if not opt_converged:
                # 失败时仅打印日志，不记录状态
                continue
        # ... 执行ORCA SP
        # 仅在成功时添加到records
        records.append({...})
    
    # 仅保存成功的记录到conformer_thermo.csv
    output_file = self.dft_dir / "conformer_thermo.csv"
```

#### 2.1.4 缺失的状态信息

对于每个构象，当前系统**不记录**：
1. OPT尝试次数（0次/1次/2次Rescue）
2. OPT是否收敛
3. SP计算是否完成
4. 失败原因（超时/不收敛/其他错误）
5. 计算耗时

---

### 2.2 问题二：S1续算不能正确实现

#### 2.2.1 当前续算逻辑

**文件**: `rph_core/orchestrator.py:505-530`

```python
# Step 1 续算检查
if resume_enabled and 's1' not in skip_steps and checkpoint_mgr.is_step_completed('s1'):
    product_xyz = checkpoint_mgr.get_step_output('s1', 'product_xyz')
    if product_xyz and Path(product_xyz).exists():
        # ... 恢复S1输出
        skip_steps.append('s1')
```

**文件**: `rph_core/steps/conformer_search/engine.py:220-234`

```python
# ConformerEngine中的缓存检查
global_min_path = self.molecule_dir / f"{self.molecule_name}_global_min.xyz"
if self._is_valid_xyz(global_min_path):
    self.logger.info(f"[S1] ⏭️ Found {global_min_path.name}, skipping conformer search")
    energy = self._extract_energy_from_xyz(global_min_path)
    return global_min_path, energy
```

#### 2.2.2 续算缺陷分析

| 场景 | 期望行为 | 实际行为 |
|------|---------|---------|
| CREST已完成，DFT中断 | 跳过CREST，继续未完成构象 | 仅检查`global_min.xyz`，如存在则跳过全部，否则重新开始 |
| 部分构象OPT成功 | 跳过成功构象，继续剩余 | 无记录，全部重新计算 |
| 某构象OPT失败 | 可重试失败构象或跳过 | 失败信息丢失，下次可能重复失败 |
| 两阶段CREST中断 | 从断点继续 | 仅检查`ensemble.xyz`，中间状态丢失 |

#### 2.2.3 关键缺陷代码

**两阶段CREST状态丢失** (`engine.py:365-455`):
```python
def _step_two_stage_crest(self, input_xyz: Path) -> Path:
    # Stage 1: GFN0 - 无状态记录
    gfn0_ensemble = self._run_crest_stage(...)
    
    # Stage 1: ISOSTAT聚类 - 无状态记录
    gfn0_clustered = self._run_isostat_clustering(...)
    
    # Stage 2: GFN2 - 无状态记录
    gfn2_ensemble = self._run_crest_batch_optimization(...)
    
    # Stage 2: ISOSTAT聚类 - 无状态记录
    gfn2_clustered = self._run_isostat_clustering(...)
    
    # 仅复制最终ensemble.xyz
    shutil.copy(gfn2_clustered, final_ensemble)
```

**问题**: 如果Stage 2中断，下次运行时会：
1. 发现`ensemble.xyz`不存在
2. 重新执行完整两阶段流程
3. 即使Stage 1已完成也会重复

---

## 3. 问题影响评估

### 3.1 计算资源浪费

| 阶段 | 平均耗时 | 重复计算影响 |
|------|---------|-------------|
| CREST GFN0 | 5-30分钟 | 高（广泛采样阶段） |
| CREST GFN2 | 10-60分钟 | 高（精细优化阶段） |
| DFT OPT (每构象) | 30-120分钟 | 极高（最昂贵阶段） |
| DFT SP (每构象) | 20-60分钟 | 极高（高精度计算） |

**场景示例**: 
- 分子有10个构象需要DFT计算
- 计算到第8个构象时中断
- 当前实现：再次运行时需重新计算全部10个构象
- 理想实现：仅需计算剩余2个构象

### 3.2 用户体验问题

1. **无法监控进度**: 用户无法知道当前处理到哪个构象
2. **故障排查困难**: 失败后无法定位具体失败的构象和原因
3. **长时间任务风险**: 大分子可能需要数天，中断后全部重来无法接受

---

## 4. 建议解决方案

### 4.1 方案概述

引入**细粒度状态追踪系统**，实现：
1. 每个构象的详细状态记录
2. 每个计算阶段（OPT/SP）的完成标记
3. 从中断点精确恢复的能力

### 4.2 新增状态文件结构

建议在 `S1_ConfGeneration/<molecule>/` 下创建 `conformer_state.json`:

```json
{
  "version": "1.0",
  "molecule_name": "product",
  "smiles": "C1=CCCCC1",
  "created_at": "2026-03-13T10:00:00",
  "updated_at": "2026-03-13T10:30:00",
  
  "crest": {
    "stage1_gfn0": {
      "status": "completed",
      "output": "xtb2/stage1_gfn0/crest_conformers.xyz",
      "completed_at": "2026-03-13T10:05:00"
    },
    "stage1_clustering": {
      "status": "completed",
      "output": "xtb2/stage1_gfn0/cluster/cluster.xyz",
      "n_conformers": 15
    },
    "stage2_gfn2": {
      "status": "completed", 
      "output": "xtb2/stage2_gfn2/crest_ensemble.xyz",
      "completed_at": "2026-03-13T10:15:00"
    },
    "stage2_clustering": {
      "status": "completed",
      "output": "xtb2/stage2_gfn2/cluster/cluster.xyz",
      "n_conformers": 8
    }
  },
  
  "conformers": [
    {
      "id": "conf_000",
      "source": "cluster/cluster.xyz#1",
      "status": "completed",
      "attempts": [
        {
          "type": "opt",
          "status": "converged",
          "log_file": "dft/conf_000.log",
          "started_at": "2026-03-13T10:16:00",
          "completed_at": "2026-03-13T10:26:00"
        },
        {
          "type": "sp", 
          "status": "completed",
          "energy_hartree": -234.56789012,
          "output_file": "dft/conf_000_SP.out",
          "completed_at": "2026-03-13T10:30:00"
        }
      ]
    },
    {
      "id": "conf_001",
      "source": "cluster/cluster.xyz#2", 
      "status": "failed",
      "attempts": [
        {
          "type": "opt",
          "status": "failed",
          "attempt": 0,
          "error": "SCF not converged",
          "log_file": "dft/conf_001.log"
        },
        {
          "type": "opt_rescue",
          "status": "failed", 
          "attempt": 1,
          "error": "SCF not converged",
          "log_file": "dft/conf_001_Res.log"
        }
      ]
    }
  ],
  
  "summary": {
    "total_conformers": 8,
    "completed": 6,
    "failed": 2,
    "best_conformer": "conf_000",
    "global_min_energy": -234.56789012,
    "global_min_xyz": "product_global_min.xyz"
  }
}
```

### 4.3 代码修改建议

#### 4.3.1 新增 `ConformerStateManager` 类

```python
# rph_core/steps/conformer_search/state_manager.py

class ConformerStateManager:
    """管理单个分子的构象计算状态"""
    
    def __init__(self, molecule_dir: Path, molecule_name: str):
        self.state_file = molecule_dir / "conformer_state.json"
        self.state = self._load_or_create()
    
    def mark_crest_stage_complete(self, stage: str, output_file: Path, metadata: dict):
        """标记CREST阶段完成"""
        pass
    
    def start_conformer_opt(self, conf_id: str, source_xyz: Path):
        """开始构象OPT计算"""
        pass
    
    def mark_opt_complete(self, conf_id: str, attempt: int, 
                         converged: bool, log_file: Path, error: str = None):
        """标记OPT尝试完成"""
        pass
    
    def mark_sp_complete(self, conf_id: str, energy: float, 
                        output_file: Path, error: str = None):
        """标记SP计算完成"""
        pass
    
    def get_pending_conformers(self) -> List[str]:
        """获取待处理的构象列表"""
        pass
    
    def get_completed_conformers(self) -> List[dict]:
        """获取已完成构象的列表"""
        pass
```

#### 4.3.2 修改 `ConformerEngine`

在关键位置插入状态记录：

```python
def run(self, smiles: str) -> Tuple[Path, float]:
    # ... 初始化代码
    state_mgr = ConformerStateManager(self.molecule_dir, self.molecule_name)
    
    # 检查是否有可恢复的状态
    if state_mgr.can_resume():
        pending_conformers = state_mgr.get_pending_conformers()
        if not pending_conformers:
            # 全部完成，直接返回
            return state_mgr.get_global_min()
    
    # CREST阶段 - 检查状态后决定是否跳过
    if not state_mgr.is_crest_complete():
        # 执行CREST...
        state_mgr.mark_crest_stage_complete("stage2_clustering", ...)
    
    # DFT阶段 - 只处理未完成的构象
    pending = state_mgr.get_pending_conformers()
    for conf_id in pending:
        state_mgr.start_conformer_opt(conf_id, ...)
        # ... 执行OPT/SP
        state_mgr.mark_opt_complete(conf_id, ...)
```

### 4.4 日志改进建议

#### 4.4.1 结构化日志输出

```python
# 在 _step_dft_opt_sp_coupled 中
for idx, xyz_file in enumerate(candidates):
    conf_id = f"conf_{idx:03d}"
    
    # 检查是否已完成
    if state_mgr.is_conformer_complete(conf_id):
        self.logger.info(f"[S1] ⏭️ {conf_id} already processed, skipping")
        continue
    
    self.logger.info(f"[S1] 🔄 Processing {conf_id} ({idx+1}/{len(candidates)})")
    
    for attempt in range(2):
        self.logger.info(f"[S1]    Attempt {attempt}: OPT starting...")
        # ... 执行OPT
        if opt_converged:
            self.logger.info(f"[S1]    ✓ OPT converged")
            break
        else:
            self.logger.warning(f"[S1]    ✗ OPT failed: {error_msg}")
    
    if opt_converged:
        self.logger.info(f"[S1]    🔄 SP starting...")
        # ... 执行SP
        if sp_energy:
            self.logger.info(f"[S1]    ✓ SP completed: E={sp_energy:.6f} Ha")
        else:
            self.logger.error(f"[S1]    ✗ SP failed")
```

#### 4.4.2 实时状态报告

```python
# 新增定期状态报告
def _log_progress_summary(self, state_mgr: ConformerStateManager):
    """输出当前进度摘要"""
    summary = state_mgr.get_summary()
    self.logger.info(
        f"[S1] 📊 Progress: {summary['completed']}/{summary['total']} conformers, "
        f"{summary['failed']} failed, "
        f"best E={summary.get('best_energy', 'N/A')}"
    )
```

---

## 5. 实施优先级

| 优先级 | 改进项 | 预估工作量 | 影响程度 |
|-------|--------|-----------|---------|
| **P0** | 构象级状态追踪系统 | 3-4天 | 极高（解决续算问题） |
| **P1** | CREST阶段状态记录 | 1-2天 | 高（避免重复采样） |
| **P2** | 结构化日志输出 | 1天 | 中（改善可观测性） |
| **P3** | 实时进度报告 | 0.5天 | 低（用户体验提升） |

---

## 6. 附录：相关文件清单

### 6.1 核心文件

| 文件路径 | 行数 | 职责 |
|---------|------|------|
| `rph_core/steps/conformer_search/engine.py` | 1434 | 构象搜索引擎主实现 |
| `rph_core/steps/anchor/handler.py` | 397 | S1锚定阶段协调器 |
| `rph_core/orchestrator.py` | 1451 | Pipeline主控制器 |
| `rph_core/utils/checkpoint_manager.py` | 604 | Pipeline级断点续传 |

### 6.2 关键函数

| 函数 | 文件 | 行号 | 职责 |
|-----|------|-----|------|
| `ConformerEngine.run()` | engine.py | 197 | 主入口 |
| `ConformerEngine._step_two_stage_crest()` | engine.py | 365 | 两阶段CREST |
| `ConformerEngine._step_dft_opt_sp_coupled()` | engine.py | 690 | DFT OPT-SP循环 |
| `ConformerEngine._run_gaussian_opt()` | engine.py | 984 | Gaussian OPT执行 |
| `AnchorPhase.run()` | handler.py | 98 | 分子锚定协调 |
| `CheckpointManager.is_step_completed()` | checkpoint_manager.py | 134 | 步骤完成检查 |

---

## 7. 结论

当前S1步骤的日志和续算功能存在根本性缺陷：

1. **日志问题**: 状态记录仅到分子级别，缺少构象级别的OPT/SP详细状态
2. **续算问题**: 无法从中间状态恢复，导致昂贵计算（CREST/DFT）的重复执行

建议优先实施**构象级状态追踪系统**，这将：
- 消除重复计算，节省大量计算资源
- 提高长时间任务的可靠性
- 改善用户监控和故障排查体验

---

*报告生成者: Sisyphus AI*
*基于代码版本: ReactionProfileHunter v6.2*
