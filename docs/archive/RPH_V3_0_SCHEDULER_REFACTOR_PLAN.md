# RPH V3.0 调度架构重构方案

**Project**: ReactionProfileHunter  
**Document**: V3.0 Scheduler Refactor Plan  
**Status**: Proposed implementation plan  
**Scope**: Dataset/S0/Small-molecule cache/Precursor/DR branch/Condition scheduling  
**Primary goal**: 最小化重复 QC 计算，最大化复用已有计算结果，同时保持 S1-S4 引擎主体稳定。

---

## 1. 背景与问题定义

当前 RPH 已经完成了 S0 DR 机理补全、dataset 按 `reaction_id` 分组、condition 隔离，以及 DR branch 调度。但当前 DR branch 调度仍存在一个核心问题：

```text
BR_MAJOR 跑完整 S1-S4
BR_DR_001 再跑完整 S1-S4，skip_steps 只跳过 s0
```

这会导致相同 precursor、DMDO、acetone、小分子热力学项在不同反应/分支目录中重复计算。对于 20 条数据集、12 个 unique reaction、默认存在 DR branch 的测试模式，这种重复会迅速放大。

V3.0 的重构目标不是重写 S1-S4 计算引擎，而是重构 **任务调度层、目录归属、artifact 引用关系**。

---

## 2. V3.0 核心计算语义

V3.0 的计算粒度定义如下：

```text
小分子：全局一次
S0：每个 reaction 一次
precursor S1：每个 reaction 一次
product S1/S2/S3/S4：每个 DR branch 一次
condition thermo/merge/DR：每个 condition 一次
```

也就是说：

1. **小分子缓存全局化**：DMDO、acetone、AcOH 等由 scheduler 在所有 reaction 前预计算。
2. **S0 前置**：S0 不再只是 `run_pipeline()` 内部步骤，而是 scheduler 的任务展开源头。
3. **reaction 与 condition 继续隔离**：当前按 `reaction_id` 去重 condition 的设计保留。
4. **DR branch 改为 precursor + branches 模型**：相同 precursor 只算一次，每个 branch 只算 product 端路径。
5. **condition 后处理 branch-aware**：每个 condition 对每个 branch 单独做温度相关 thermo，随后生成 condition-level DR。

---

## 3. 目标目录契约

V3.0 目标目录不采用过度臃肿的 `shared/molecules/pathways/ts` 层，而采用更符合化学语义的结构：

```text
rph_output/
├── small_molecules/
│   ├── <DMDO_hash>/
│   │   ├── molecule_min.xyz
│   │   ├── finalDFT/
│   │   ├── thermo.json
│   │   └── cache_meta.json
│   └── <acetone_hash>/
│       ├── molecule_min.xyz
│       ├── finalDFT/
│       ├── thermo.json
│       └── cache_meta.json
│
└── RXN_xxxxxxxx/
    ├── reaction_manifest.json
    │
    ├── S0_Mechanism/
    │   ├── mechanism_graph.json
    │   ├── mechanism_summary.json
    │   ├── atom_map_smiles.json
    │   ├── dr_branch_plan.json
    │   └── s0_status.json
    │
    ├── precursor/
    │   ├── precursor_manifest.json
    │   └── S1_ConfGeneration/
    │       └── precursor/
    │
    ├── branches/
    │   ├── BR_MAJOR/
    │   │   ├── branch_manifest.json
    │   │   ├── S1_ConfGeneration/
    │   │   │   ├── product/
    │   │   │   ├── precursor -> ../../../precursor/S1_ConfGeneration/precursor
    │   │   │   └── small_molecules/
    │   │   ├── S2_Retro/
    │   │   ├── S3_TS/
    │   │   └── S4_Data/
    │   │
    │   └── BR_DR_001/
    │       ├── branch_manifest.json
    │       ├── S1_ConfGeneration/
    │       │   ├── product/
    │       │   ├── precursor -> ../../../precursor/S1_ConfGeneration/precursor
    │       │   └── small_molecules/
    │       ├── S2_Retro/
    │       ├── S3_TS/
    │       └── S4_Data/
    │
    ├── reaction_features/
    │   └── dr_prediction_default.json
    │
    └── conditions/
        └── COND_xxx/
            ├── condition_manifest.json
            ├── branches/
            │   ├── BR_MAJOR/
            │   │   ├── thermo_calculation.json
            │   │   └── merged_features.csv
            │   └── BR_DR_001/
            │       ├── thermo_calculation.json
            │       └── merged_features.csv
            └── dr_prediction.json
```

### 3.1 Artifact 归属表

| Artifact | 所属层级 | 生产者 | 消费者 | 重算策略 |
|---|---|---|---|---|
| `small_molecules/<hash>/molecule_min.xyz` | 全局 | small molecule precompute | S1/S4/condition thermo | 同 SMILES + theory 只算一次 |
| `small_molecules/<hash>/thermo.json` | 全局 | small molecule precompute | `ConditionThermoCalculator` | 同 SMILES + theory + thermo config 只算一次 |
| `S0_Mechanism/dr_branch_plan.json` | reaction | S0 | V3 scheduler | 每个 reaction 一次 |
| `precursor/S1_ConfGeneration/precursor/` | reaction | precursor S1 | all branches/S4 | 每个 reaction 一次 |
| `branches/BR_*/S1_ConfGeneration/product/` | branch | product S1 | S2/S3/S4 | 每个 branch 一次 |
| `branches/BR_*/S2_Retro/` | branch | S2 | S3/S4 | 每个 branch 一次 |
| `branches/BR_*/S3_TS/` | branch | S3 | S4/condition thermo | 每个 branch 一次 |
| `branches/BR_*/S4_Data/` | branch | S4 | merger/DR | 每个 branch 一次 |
| `conditions/COND_*/branches/BR_*/thermo_calculation.json` | condition-branch | condition thermo | feature merger/DR | 每个 condition × branch 一次 |
| `conditions/COND_*/dr_prediction.json` | condition | DR aggregator | user/ML dataset | 每个 condition 一次 |

---

## 4. 新增模块设计

### 4.1 `rph_core/scheduling/models.py`

新增 V3.0 调度数据模型。

```python
@dataclass
class ReactionJob:
    reaction_id: str
    reaction_root: Path
    representative: TaskSpec
    conditions: list[TaskSpec]
    dr_plan_path: Path
    branches: list[BranchJob]

@dataclass
class BranchJob:
    parent_reaction_id: str
    branch_id: str
    pathway_id: str
    product_smiles: str
    branch_root: Path
    generation_policy: str
    flipped_map_numbers: list[int]
    fixed_stereocenters: list[int]
    notes: list[str]

@dataclass
class ConditionJob:
    reaction_id: str
    condition_id: str
    condition_root: Path
    task: TaskSpec
    temperature_k: float
```

设计要求：

- 所有路径使用 `Path`。
- 不把 S1/S2/S3/S4 结果直接塞进 dataclass，只记录任务身份和根目录。
- `BranchJob` 必须来自 S0 的 `dr_branch_plan.json`，不能由下游自行推断。

---

### 4.2 `rph_core/scheduling/v3_scheduler.py`

新增 V3.0 主调度器。

核心职责：

1. 调用 `build_tasks_from_run_config()`。
2. 按 `reaction_id` 分组。
3. 为每个 reaction 创建 `reaction_root`。
4. 前置运行 S0。
5. 读取 `dr_branch_plan.json`。
6. 展开 `ReactionJob / BranchJob / ConditionJob`。
7. 调用 small molecule precompute。
8. 按 reaction 顺序执行：
   - precursor S1
   - branch product S1-S4
   - condition branch-aware thermo/merge/DR

建议类结构：

```python
class V3Scheduler:
    def __init__(self, hunter: ReactionProfileHunter, run_cfg: dict[str, Any]): ...

    def plan(self) -> list[ReactionJob]: ...

    def precompute_small_molecules(self, jobs: list[ReactionJob]) -> None: ...

    def run(self) -> list[PipelineResult]: ...

    def _run_reaction(self, job: ReactionJob) -> list[PipelineResult]: ...

    def _run_precursor_s1(self, job: ReactionJob) -> Path: ...

    def _run_branch(self, job: ReactionJob, branch: BranchJob) -> PipelineResult: ...

    def _run_conditions(self, job: ReactionJob) -> None: ...
```

约束：

- 不要在这里实现 QC 细节。
- 不要绕过 `ReactionProfileHunter` 已有 S1/S2/S3/S4 engine。
- 不要在 scheduler 中硬编码 Gaussian/xTB/ORCA 行为。

---

### 4.3 `rph_core/scheduling/small_molecule_precompute.py`

新增小分子预计算模块。

职责：

1. 从 dataset 的 `small_molecular_keys` 收集小分子。
2. 从 `reaction_reference_terms` 收集 DMDO/acetone/AcOH 等热力学引用物种。
3. 通过 `SmallMoleculeCatalog` 解析 SMILES。
4. 调用现有 `AnchorPhase` 或 molecule-level S1 runner 计算小分子。
5. 写入全局 cache。
6. 保证 `thermo.json` 存在。

建议接口：

```python
class SmallMoleculePrecomputer:
    def __init__(self, config: dict[str, Any], cache_root: Path): ...

    def collect_keys(self, reaction_jobs: list[ReactionJob]) -> set[str]: ...

    def ensure_all(self, keys: set[str]) -> dict[str, Path]: ...

    def ensure_one(self, key: str) -> Path: ...

    def ensure_thermo_json(self, cache_entry: Path) -> Path: ...
```

实现要求：

- 复用 `SmallMoleculeCache.acquire_compute_lock()`，避免并发 cache stampede。
- `cache_meta.json` 必须包含 theory signature。
- `thermo.json` 缺失时从 `finalDFT/*Shermo*.sum` 派生。
- 如果 Shermo 产物不可用，应写 warning 到 cache metadata，而不是静默通过。

---

### 4.4 `rph_core/scheduling/artifact_refs.py`

新增 artifact 引用工具，避免各处手写相对路径。

职责：

```python
def link_or_copy_precursor_into_branch(reaction_root: Path, branch_root: Path) -> None: ...

def materialize_small_molecule_refs(branch_root: Path, cache_root: Path, keys: Iterable[str]) -> None: ...

def read_branch_manifest(branch_root: Path) -> dict[str, Any]: ...

def write_branch_manifest(branch_root: Path, payload: dict[str, Any]) -> Path: ...
```

实现策略：

- Linux 下优先 symlink。
- Windows/受限环境 fallback 到 copy。
- 所有链接都要写入 manifest，便于 debug。

---

## 5. 现有文件逐文件修改方案

### 5.1 `config/defaults.yaml`

新增配置：

```yaml
scheduler:
  mode: legacy          # legacy | v3
  v3:
    precompute_small_molecules: true
    precursor_s1_once_per_reaction: true
    branch_product_only_s1: true
    condition_branch_thermo: true
```

保留 legacy 默认值，避免直接影响当前可运行 pipeline。

验收点：

- `scheduler.mode: legacy` 时行为与当前一致。
- `scheduler.mode: v3` 时进入 `_run_tasks_v3()`。

---

### 5.2 `rph_core/orchestrator.py`

#### Step 1：保留 `_run_tasks()`，新增 `_run_tasks_v3()`

```python
def _run_tasks_v3(hunter: ReactionProfileHunter, run_cfg: dict[str, Any]) -> list[PipelineResult]:
    from rph_core.scheduling.v3_scheduler import V3Scheduler
    return V3Scheduler(hunter=hunter, run_cfg=run_cfg).run()
```

#### Step 2：在 `main()` 或任务入口处选择 scheduler

伪代码：

```python
scheduler_cfg = config.get("scheduler", {}) or {}
mode = str(scheduler_cfg.get("mode", "legacy")).lower()
if mode == "v3":
    results = _run_tasks_v3(hunter, run_cfg)
else:
    results = _run_tasks(hunter, run_cfg)
```

#### Step 3：抽出 S0 单独运行 helper

当前 `_run_s0()` 是 `ReactionProfileHunter` 方法，可直接复用，但需要一个调度层 helper 包装：

```python
def run_s0_for_planning(...):
    checkpoint_mgr = CheckpointManager(reaction_root)
    return hunter._run_s0(...)
```

注意：

- 不要复制 S0 逻辑。
- 只包装 checkpoint/progress 上下文。
- 输出仍写入 `reaction_root/S0_Mechanism/`。

验收点：

- V3 plan 阶段能生成 `dr_branch_plan.json`。
- S0 disabled/degraded 时仍能生成有效 plan 或明确 warning。

---

### 5.3 `rph_core/steps/anchor/handler.py`

当前 `AnchorPhase.run()` 已经按 molecule 字典运行，并且已有小分子 cache 被动复用逻辑。

需要改造：

#### Step 1：新增单分子公开入口

```python
def run_single_molecule(self, name: str, smiles: str, base_work_dir: Optional[Path] = None) -> AnchorPhaseResult:
    return self.run({name: smiles}, base_work_dir=base_work_dir)
```

目的：

- precursor S1 可以单独运行。
- product S1 可以在 branch 中单独运行。
- small molecule precomputer 可以复用相同入口。

#### Step 2：小分子 cache 写入补齐 `thermo.json`

当前 `_store_small_molecule_cache()` 写入：

```text
molecule_min.xyz
finalDFT/
cache_meta.json
```

V3.0 需要补充：

```text
thermo.json
```

建议逻辑：

```python
sum_file = find best *_Shermo.sum under finalDFT/
if sum_file:
    derive_hoac_thermo_from_sum(sum_file, cache_dir / "thermo.json")
else:
    cache_meta["warnings"].append("W_MISSING_SHERMO_SUM_FOR_THERMO_JSON")
```

验收点：

- DMDO/acetone cache entry 有 `thermo.json`。
- `ConditionThermoCalculator._add_reference_species()` 能读到对应物种。

---

### 5.4 `rph_core/utils/small_molecule_cache.py`

保留现有类，增加少量辅助方法即可。

新增方法建议：

```python
def is_complete(self, smiles: str, theory_signature: Optional[dict[str, Any]] = None, require_thermo: bool = True) -> bool: ...

def get_entry_manifest(self, smiles: str) -> dict[str, Any]: ...
```

`exists()` 当前允许 `molecule_min.xyz` 或 `thermo.json` 任一存在，V3.0 预计算需要更严格的 completeness 判断。

验收点：

- 小分子预热使用 `is_complete(require_thermo=True)`。
- 旧 S1 被动 cache 命中仍可使用 `exists()`，保持兼容。

---

### 5.5 `rph_core/utils/path_manager.py`

新增 V3.0 路径 helper：

```python
def get_global_small_molecules_root(output_root: Path) -> Path:
    return Path(output_root) / "small_molecules"

def get_precursor_root(reaction_root: Path) -> Path:
    return Path(reaction_root) / "precursor"

def get_precursor_s1_dir(reaction_root: Path) -> Path:
    return get_precursor_root(reaction_root) / "S1_ConfGeneration"

def get_condition_branch_root(condition_root: Path, branch_id: str) -> Path:
    return Path(condition_root) / "branches" / str(branch_id)
```

验收点：

- 不在 scheduler 中手写路径。
- 所有新路径有单元测试。

---

### 5.6 `rph_core/steps/condition_thermo.py`

当前类签名：

```python
ConditionThermoCalculator(config, reaction_root, condition_root)
```

内部硬编码：

```python
s3_dir = self.reaction_root / "S3_TS"
```

V3.0 需要改为 branch-aware：

```python
@dataclass(frozen=True)
class ConditionThermoCalculator:
    config: Dict[str, Any]
    reaction_root: Path
    condition_root: Path
    artifact_root: Optional[Path] = None
```

内部：

```python
root = self.artifact_root or self.reaction_root
s3_dir = root / "S3_TS"
```

这样 legacy 仍可用，V3 传 `artifact_root=branch_root`。

验收点：

- legacy tests 不变。
- V3 condition thermo 可从 `branches/BR_*/S3_TS/` 读取。

---

### 5.7 `rph_core/steps/condition_feature_merger.py`

当前类签名：

```python
ConditionFeatureMerger(reaction_root, condition_root)
```

内部固定读：

```python
reaction_features_dir = self.reaction_root / "reaction_features"
```

V3.0 需要支持 branch feature source：

```python
@dataclass(frozen=True)
class ConditionFeatureMerger:
    reaction_root: Path
    condition_root: Path
    feature_root: Optional[Path] = None
```

读取策略：

1. 如果 `feature_root` 提供，则优先从：
   - `feature_root / "S4_Data" / "features_raw.csv"`
   - 或 `feature_root / "reaction_features" / "geo_electronic_features.json"`
2. 否则保持 legacy：
   - `reaction_root / "reaction_features"`

验收点：

- 每个 `conditions/COND_*/branches/BR_*` 都有自己的 `merged_features.csv`。
- merged features 中包含 `branch_id` 字段。

---

### 5.8 `rph_core/steps/dr_aggregator.py`

当前 aggregator 从：

```text
reaction_root/S4_Data/features_raw.csv
reaction_root/branches/*/S4_Data/features_raw.csv
```

读取 branch barrier。

V3.0 需要增加 condition-aware 模式：

```python
def run_for_condition(
    self,
    reaction_root: Path,
    condition_root: Path,
    temperature_k: float,
) -> dict[str, Any]: ...
```

读取：

```text
conditions/COND_xxx/branches/BR_*/merged_features.csv
```

优先列：

```text
thermo.dG_activation
```

fallback：

```text
branches/BR_*/S4_Data/features_raw.csv
```

验收点：

- 同一 reaction 不同 temperature 的 `dr_prediction.json` 可不同。
- 缺失某 branch 时 degraded，不中断全局任务。

---

### 5.9 `rph_core/utils/task_builder.py`

当前已有 `TaskSpec` 和 `BranchTaskSpec`，但 `BranchTaskSpec` 未被 DR 调度使用。

改造建议：

- 保留 `TaskSpec`。
- 可将 `BranchTaskSpec` 迁移或复用于 `scheduling.models.BranchJob`。
- 不建议在 `task_builder.py` 中展开 branch，因为 branch 必须来自 S0 `dr_branch_plan.json`，而 `task_builder.py` 不应运行 S0。

验收点：

- `task_builder.py` 只负责 dataset → TaskSpec。
- branch 展开归 V3 scheduler。

---

## 6. 分阶段实施计划

### Phase 0：配置与兼容入口

**目标**：引入 V3 模式，但不改变 legacy 行为。

修改文件：

1. `config/defaults.yaml`
   - 新增 `scheduler.mode: legacy`。
2. `rph_core/orchestrator.py`
   - 新增 `_run_tasks_v3()` 空壳。
   - 主入口根据 `scheduler.mode` 选择 legacy/v3。

测试：

```bash
python scripts/ci/check_imports.py rph_core
pytest tests/test_tsv_loader.py -v
```

验收：

- 默认仍走 legacy。
- 现有 dataset run 行为不变。

---

### Phase 1：V3 Planner 与 S0 前置

**目标**：先只生成计划，不跑 S1-S4。

新增文件：

1. `rph_core/scheduling/__init__.py`
2. `rph_core/scheduling/models.py`
3. `rph_core/scheduling/v3_scheduler.py`

修改文件：

1. `rph_core/orchestrator.py`
   - `_run_tasks_v3()` 调用 `V3Scheduler.plan()`。

测试建议新增：

```text
tests/test_v3_scheduler_plan.py
```

测试点：

- 20 条 dataset → 12 个 `ReactionJob`。
- 每个 `ReactionJob` 有 `conditions`。
- 每个 reaction 生成一个 `S0_Mechanism/dr_branch_plan.json`。
- branch 数量来自 `dr_branch_plan.json`。

---

### Phase 2：小分子全局预计算

**目标**：所有 reaction 前先完成 DMDO/acetone 等小分子 cache。

新增文件：

1. `rph_core/scheduling/small_molecule_precompute.py`

修改文件：

1. `rph_core/utils/small_molecule_cache.py`
   - 新增 `is_complete()`。
2. `rph_core/steps/anchor/handler.py`
   - 小分子 cache 补齐 `thermo.json`。
3. `rph_core/scheduling/v3_scheduler.py`
   - `precompute_small_molecules()`。

测试建议新增：

```text
tests/test_v3_small_molecule_precompute.py
```

测试点：

- DMDO/acetone 在多个 reaction 中只计算一次。
- cache entry 包含 `molecule_min.xyz`、`finalDFT/`、`thermo.json`、`cache_meta.json`。
- cache theory signature 不匹配时触发 recompute。

---

### Phase 3：precursor S1 reaction-level 拆分

**目标**：每个 reaction 的 precursor 只跑一次。

修改文件：

1. `rph_core/steps/anchor/handler.py`
   - 新增 `run_single_molecule()`。
2. `rph_core/utils/path_manager.py`
   - 新增 precursor path helpers。
3. `rph_core/scheduling/v3_scheduler.py`
   - 新增 `_run_precursor_s1()`。
4. `rph_core/scheduling/artifact_refs.py`
   - 新增 precursor link/copy helper。

测试建议新增：

```text
tests/test_v3_precursor_reuse.py
```

测试点：

- 一个 reaction 多个 branch 时，precursor S1 调用次数 = 1。
- 每个 branch S1 目录能看到 precursor artifact ref。
- S4 `_resolve_s1_artifacts()` 能解析 precursor。

---

### Phase 4：branch product-only pipeline

**目标**：每个 branch 只跑 product S1 + S2 + S3 + S4，不再重复 precursor 和小分子。

修改文件：

1. `rph_core/scheduling/v3_scheduler.py`
   - `_run_branch()`。
2. `rph_core/orchestrator.py`
   - 可新增轻量 branch runner，或扩展 `run_pipeline()` 支持 `s1_molecule_scope`。

推荐最小侵入实现：

```python
hunter.run_pipeline(
    product_smiles=branch.product_smiles,
    work_dir=branch.branch_root,
    precursor_smiles=None,
    small_molecular_keys=[],
    skip_steps=[...],
)
```

然后在 branch `S1_ConfGeneration/` 下 materialize precursor/small molecule refs，以兼容 S4。

测试建议新增：

```text
tests/test_v3_branch_product_only.py
```

测试点：

- `branches/BR_MAJOR/S1_ConfGeneration/product/` 存在。
- `branches/BR_MAJOR/S1_ConfGeneration/precursor` 是引用或复制。
- `branches/BR_DR_001/S2_Retro/S3_TS/S4_Data` 独立存在。
- branch 之间不复用 product TS。

---

### Phase 5：condition branch-aware 后处理

**目标**：每个 condition 对每个 branch 生成 thermo/merged features，再做 condition-level DR。

修改文件：

1. `rph_core/steps/condition_thermo.py`
   - 增加 `artifact_root`。
2. `rph_core/steps/condition_feature_merger.py`
   - 增加 `feature_root` 和 `branch_id`。
3. `rph_core/steps/dr_aggregator.py`
   - 增加 `run_for_condition()`。
4. `rph_core/scheduling/v3_scheduler.py`
   - `_run_conditions()`。

测试建议新增：

```text
tests/test_v3_condition_branch_outputs.py
tests/test_v3_condition_dr_temperature.py
```

测试点：

- `conditions/COND_xxx/branches/BR_MAJOR/merged_features.csv` 存在。
- `conditions/COND_xxx/branches/BR_DR_001/merged_features.csv` 存在。
- `conditions/COND_xxx/dr_prediction.json` 使用该 condition 的 temperature。
- 两个不同温度 condition 的 predicted DR 可不同。

---

### Phase 6：V3 全链路 dry-run 与 mock-QC 验证

**目标**：在不调用真实 QC 的情况下验证目录、manifest、调度次数。

测试建议新增：

```text
tests/test_v3_scheduler_end_to_end_mock.py
```

测试点：

- dataset 20 rows → 12 reactions。
- S0 调用次数 = 12。
- small molecule precompute 中 DMDO/acetone 各调用一次。
- precursor S1 调用次数 = reaction 数。
- product branch pipeline 调用次数 = branch 总数。
- condition thermo 调用次数 = condition 数 × branch 数。

---

## 7. 关键风险与解决策略

### 风险 1：S4 仍假设 precursor 在 branch 的 S1 目录下

解决：

- branch `S1_ConfGeneration/precursor` 使用 symlink/copy 指向 `reaction_root/precursor/S1_ConfGeneration/precursor`。
- 不立即重写所有 S4 extractor。

### 风险 2：小分子 cache 有几何但无 `thermo.json`

解决：

- `SmallMoleculePrecomputer.ensure_one()` 必须验证 `thermo.json`。
- `_store_small_molecule_cache()` 写入时尝试从 Shermo `.sum` 派生。

### 风险 3：DR branch 错误复用 TS

解决：

- V3.0 禁止跨不同 product branch 复用 S2/S3。
- 每个 branch 独立运行 product S1/S2/S3/S4。
- 只复用 precursor 和小分子。

### 风险 4：condition DR 仍使用默认 298.15 K

解决：

- `DRAggregator.run_for_condition()` 必须从 `condition_manifest.json` 或传入参数读取 temperature。
- `reaction_features/dr_prediction_default.json` 只作为默认视图，不作为 condition 数据源。

### 风险 5：大改造影响现有可运行 pipeline

解决：

- `scheduler.mode: legacy` 保持默认。
- V3.0 所有入口通过开关启用。
- 每个 phase 都保持 tests passing 后再进入下一步。

---

## 8. 推荐提交拆分

建议不要一次性提交整个 V3.0，而是拆成以下 commit/PR：

1. `add v3 scheduler mode and planning models`
2. `frontload s0 planning for v3 scheduler`
3. `precompute global small molecule cache`
4. `split precursor s1 from branch product pipeline`
5. `make condition thermo branch-aware`
6. `add condition-level dr aggregation`
7. `enable v3 scheduler end-to-end with mock tests`

每个提交都应满足：

```bash
python scripts/ci/check_imports.py rph_core
pytest <相关新增测试> -v
```

---

## 9. 最终验收标准

V3.0 完成后，对 20 条数据集应满足：

1. 20 条输入按 `reaction_id` 去重为 12 个 reaction。
2. 每个 reaction 只生成一次 `S0_Mechanism/dr_branch_plan.json`。
3. DMDO/acetone 小分子全局 cache 只计算一次。
4. 每个 reaction 的 precursor S1 只计算一次。
5. 每个 DR branch 只计算自己的 product S1/S2/S3/S4。
6. condition 不触发 S1/S2/S3 重算。
7. 每个 condition 输出 branch-aware merged features。
8. 每个 condition 输出 temperature-aware `dr_prediction.json`。
9. legacy scheduler 仍可运行，作为回滚路径。

---

## 10. 推荐执行顺序摘要

```text
Phase 0  配置开关 + legacy 保护
Phase 1  V3 planner + S0 前置
Phase 2  小分子全局预计算
Phase 3  precursor S1 reaction-level 拆分
Phase 4  branch product-only S1-S4
Phase 5  condition branch-aware thermo/merge/DR
Phase 6  mock-QC 全链路验证
```

这套方案的核心价值是：

```text
把 RPH 从“每个目录独立跑完整 pipeline”
改造成“由 S0 驱动的 reaction/branch/condition 任务图”。
```

这样可以保持现有 S1-S4 引擎稳定，同时从调度层消除最大量的重复计算。
