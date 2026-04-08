# RPH Git 工作流规范 v1.0

> 针对 ReactionProfileHunter 四层协议栈 v2.10 及后续开发

---

## 1. 分支策略：Git Flow 科研增强版

采用 **Git Flow + Feature Flags** 混合策略，适应 RPH "长期开发 + 版本基准"的特点。

### 1.1 核心分支

| 分支 | 用途 | 保护级别 | 合并要求 |
|------|------|----------|----------|
| `main` | 稳定发布分支，对应 published paper | 🔒 严格保护 | PR + 2 reviewers + CI全绿 + Gate A通过 |
| `develop` | 集成开发分支，接受功能合并 | 🔒 中等保护 | PR + 1 reviewer + CI全绿 |
| `release/*` | 发布候选（如 `release/v2.10`） | 🔒 严格保护 | 仅接受 bugfix，禁止新功能 |
| `feature/*` | 功能开发（如 `feature/pipeline-v210`） | 📝 常规保护 | 合并到 develop |
| `hotfix/*` | 紧急修复（如 `hotfix/orca-timeout`） | 🔒 严格保护 | 直接合 main + cherry-pick develop |
| `benchmark/*` | 基准测试专用（如 `benchmark/v210-vs-v290`） | 📝 常规保护 | 仅用于验证，不合并 |

### 1.2 v2.10 具体分支规划

```text
main
 ├── develop
 │    ├── feature/pipeline-skeleton        # Phase 1: Pipeline骨架 + 工具函数
 │    ├── feature/lite-zero-protocols      # Phase 2: lite/zero协议实现
 │    ├── feature/full-protocol            # Phase 3: full协议（可选）
 │    ├── feature/s2-pes-adapter           # Phase 4: S2 PES适配器
 │    ├── feature/provenance-output        # Phase 5: provenance.json规范
 │    └── feature/s4-degradation           # Phase 6: S4容错降级
 │
 ├── release/v2.10                         # 发布候选，通过三门验收后合并main
 │    ├── hotfix/pipeline-timeout          # 仅bugfix
 │    └── hotfix/provenance-missing        # 仅bugfix
 │
 └── benchmark/v210-vs-v290                # 基准验证分支（保留不删）
```

### 1.3 分支命名规范

```
feature/<kebab-case-description>[-<issue-id>]
hotfix/<kebab-case-description>
benchmark/<comparison-description>
release/v<major>.<minor>[.<patch>]

# 示例
feature/pipeline-skeleton-42          # 带issue编号
feature/lite-zero-protocols
hotfix/orca-sp-mpi-crash
benchmark/v210-ext-vs-v290-ext
release/v2.10
```

---

## 2. Commit Message 规范（Conventional Commits + RPH扩展）

### 2.1 格式

```
<type>(<scope>): <subject>

<body>

<footer>
```

### 2.2 Type 定义

| Type | 含义 | 使用场景 |
|------|------|----------|
| `feat` | 新功能 | 新增 stage、新协议 |
| `fix` | 修复 | bug修复、配置修正 |
| `refactor` | 重构 | 代码结构调整，无功能变更 |
| `perf` | 性能优化 | 缓存、并行优化 |
| `test` | 测试 | 新增/修改测试 |
| `docs` | 文档 | V2.10_Plan.md、README、注释 |
| `config` | 配置 | defaults.yaml、模板修改 |
| `benchmark` | 基准测试 | 新增 benchmark 脚本、结果 |
| `gate` | 验收门 | Gate A/B/C 相关脚本或修复 |

### 2.3 Scope 定义

```
s1              # Step 1 相关
s2              # Step 2 相关  
s3              # Step 3 相关
s4              # Step 4 相关
pipeline        # Pipeline框架
protocols       # 协议规格（ProtocolSpec）
stages/<name>   # 具体 stage
config          # 配置文件
qc/<engine>     # QC接口（orca/gaussian/xtb/crest）
docs            # 文档
```

### 2.4 示例

```
feat(pipeline): add ProtocolSpec and ConformerPipeline skeleton

- Implement ProtocolSpec dataclass with 4 predefined specs
- Add ConformerPipeline executor with stage orchestration
- Add BaseStage ABC and PipelineContext
- Include stage auto-discovery via @register_stage decorator

Relates to #42
```

```
feat(stages/lite): implement low_level_rerank with mRRHO support

- Add r2SCAN-3c SP calculation via ORCAInterface
- Add GFN2 mRRHO correction via xtb --hess
- Compute G = E_SP + G_solv(CPCM) + G_mRRHO
- Support two-step optimization: coarse SP filter → mRRHO on subset

Gate B: verified on 2 representative molecules
```

```
fix(config): increase full protocol prescreening window from 4.0 to 8.0 kcal/mol

Per review feedback: 4.0 window too aggressive, may discard true global min.
Official CENSO recommends wider initial window.

BREAKING CHANGE: full protocol now uses 8.0 kcal/mol prescreening threshold
```

```
gate(a): add ext equivalence regression test suite

- Compare product_min.xyz coordinates (RMSD < 1e-6 Å)
- Compare final SP energy (< 1e-6 Hartree)
- Compare conformer_thermo.csv content
- Compare conformer_energies.json structure

Fails CI if ext output differs from baseline.
```

---

## 3. Pull Request 规范

### 3.1 PR 模板

```markdown
## 变更类型
- [ ] feat: 新功能
- [ ] fix: 修复
- [ ] refactor: 重构
- [ ] perf: 性能优化
- [ ] test: 测试
- [ ] docs: 文档
- [ ] config: 配置
- [ ] gate: 验收门

## 关联 Issue
Closes #<issue-number>

## 变更摘要
<!-- 一句话描述 -->

## 详细变更
<!-- 分点列出关键变更 -->

## 验收门状态
- [ ] Gate A: ext 等价性测试通过
- [ ] Gate B: 新协议输出契约测试通过  
- [ ] Gate C: 下游稳定性测试通过

## 基准测试影响
- [ ] 不影响现有基准
- [ ] 新增基准测试（附结果）
- [ ] 修改基准参数（说明原因）

## 文档更新
- [ ] V2.10_Plan.md 已更新
- [ ] defaults.yaml 已更新
- [ ] 代码注释已更新

## 破坏性变更
<!-- 如有，列出并说明迁移路径 -->
```

### 3.2 PR 分类标签

| 标签 | 用途 | 颜色 |
|------|------|------|
| `protocol/ext` | 涉及 ext 协议 | 🔵 blue |
| `protocol/full` | 涉及 full 协议 | 🟣 purple |
| `protocol/lite` | 涉及 lite 协议 | 🟢 green |
| `protocol/zero` | 涉及 zero 协议 | 🟡 yellow |
| `stage/*` | 具体 stage | ⚪ gray |
| `gate-a` | Gate A 相关 | 🔴 red |
| `gate-b` | Gate B 相关 | 🔴 red |
| `gate-c` | Gate C 相关 | 🔴 red |
| `benchmark` | 基准测试 | 🟠 orange |
| `breaking` | 破坏性变更 | 🔴 red |
| `docs-only` | 仅文档 | ⚫ black |

### 3.3 Review 规则

| 目标分支 | 必需 Reviewers | 特殊要求 |
|----------|----------------|----------|
| `main` | 2人（含代码owner） | 必须通过 Gate A/B/C |
| `develop` | 1人 | CI全绿 |
| `release/*` | 2人 | 仅 bugfix，禁止新功能 |
| `hotfix/*` | 1人（紧急情况） | 事后补PR到develop |

---

## 4. 版本管理

### 4.1 版本号规范（SemVer）

```
v<major>.<minor>.<patch>[-<prerelease>][+<build>]

# 示例
v2.10.0           # 正式发布
v2.10.0-rc1       # 发布候选1
v2.10.0-beta.2    # Beta 2
v2.10.0+build.123 # 构建号
```

### 4.2 版本标签管理

```bash
# 发布标签
git tag -a v2.10.0 -m "Release v2.10.0: Four-layer protocol stack

- Add ProtocolSpec + ConformerPipeline architecture
- Add lite protocol (r2SCAN-3c + mRRHO + final SP)
- Add zero protocol (GFN2+mRRHO + narrow window + final SP)
- Add full protocol implementation (code ready, not exposed)
- Add S2 PES adapter for xTB/DFT mismatch
- Add provenance.json output contract
- Three-methodology firewalls enforced
- Three-gate acceptance criteria

See V2.10_Plan.md for details."

# 推送标签
git push origin v2.10.0
```

### 4.3 发布检查清单

- [ ] 所有 Gate A/B/C 通过
- [ ] `CHANGELOG.md` 已更新
- [ ] `V2.10_Plan.md` 状态标记为 "Released"
- [ ] 版本标签已打
- [ ] GitHub Release 已创建（附二进制/文档）
- [ ] Docker 镜像已构建推送
- [ ] Paper 附录已同步

---

## 5. 基准测试分支管理

### 5.1 专用 benchmark 分支

```bash
# 创建基准验证分支
git checkout -b benchmark/v210-lite-vs-ext main

# 运行基准测试（结果提交到分支，不合并）
python scripts/benchmark/run_suite.py --protocols ext,lite --molecules test_set_20.json

# 提交结果
 git add benchmark/results/v210/
 git commit -m "benchmark(v210): add lite vs ext comparison on 20 molecules
 
 Results:
 - ΔG‡ MAE: 0.8 kcal/mol (< 1.0 target ✅)
 - S2 success rate: 95% (ext: 96%, not degraded ✅)
 - Wall time reduction: 70% vs ext"

# 推送保留（不创建 PR）
git push origin benchmark/v210-lite-vs-ext
```

### 5.2 基准结果归档

```
benchmark/
├── results/
│   ├── v210/
│   │   ├── lite-vs-ext-20mol/
│   │   ├── zero-vs-ext-20mol/
│   │   └── full-vs-ext-10mol/     # full 内部验证
│   └── v290/                       # 下一版本基准
├── scripts/
│   ├── run_suite.py
│   └── analyze_results.py
└── README.md                       # 如何运行基准测试
```

---

## 6. 大型重构的 Git 策略（v2.10 专用）

### 6.1 增量合并策略

避免 5000+ 行的巨型 PR，采用 "骨架先行，功能填充"：

```bash
# PR 1: 纯骨架（可 review 的结构）
feature/pipeline-skeleton
├── pipeline/__init__.py        # 空壳
├── pipeline/spec.py            # ProtocolSpec 定义（无实现）
├── pipeline/pipeline.py        # ConformerPipeline 空壳
├── pipeline/base_stage.py      # BaseStage ABC
└── tests/test_pipeline_skeleton.py

# PR 2: 工具函数（独立可测）
feature/pipeline-utils
├── pipeline/utils/solvent_resolver.py
├── pipeline/utils/boltzmann.py
└── tests/test_pipeline_utils.py

# PR 3-5: stages 分批实现
feature/stages-generation       # rdkit_embed, crest_sampling
feature/stages-filtering        # isostat_clustering, ensemble_processing
feature/stages-selection        # boltzmann_cutoff, energy_window_selection

# PR 6: lite 协议集成
feature/lite-protocol-integration

# PR 7: zero 协议集成
feature/zero-protocol-integration

# PR 8: S2 PES 适配器（跨 step 边界，需额外 review）
feature/s2-pes-adapter
```

### 6.2 功能开关（Feature Flags）

在 `defaults.yaml` 中使用功能开关，允许半成品的 feature 合并到 develop：

```yaml
# config/defaults.yaml
experimental:
  protocol_stack:
    enabled: false          # 全局开关
    expose_full: false      # full 协议开关（v2.10 保持 false）
    expose_lite: true       # lite 协议开关
    expose_zero: true       # zero 协议开关
```

代码中检查：

```python
def build_s1_engine(protocol: str, config: Dict, ...):
    exp = config.get("experimental", {}).get("protocol_stack", {})
    
    if protocol == "full" and not exp.get("expose_full", False):
        raise ValueError("full protocol is experimental and not yet exposed")
    
    ...
```

优点：
- PR 可以小步快跑合并到 develop
- 未完成的功能默认关闭，不影响他人
- 可通过配置随时开启测试

---

## 7. CI/CD 集成

### 7.1 必需 CI Checks

| Check | 触发条件 | 失败后果 |
|-------|----------|----------|
| `lint-ruff` | 所有 PR | ❌ 阻止合并 |
| `type-check-mypy` | 所有 PR | ❌ 阻止合并 |
| `import-style-check` | 所有 PR | ❌ 阻止合并 |
| `unit-tests` | 所有 PR | ❌ 阻止合并 |
| `gate-a-ext-equivalence` | 涉及 s1/conformer_search | ❌ 阻止合并 |
| `gate-b-contract` | 涉及 pipeline/ | ❌ 阻止合并 |
| `gate-c-downstream` | 涉及 step2/ | ⚠️ 警告（允许合并但需记录） |
| `docs-build` | 涉及 docs/ | ❌ 阻止合并 |

### 7.2 GitHub Actions 示例结构

```yaml
# .github/workflows/ci.yml
name: CI

on:
  push:
    branches: [main, develop, release/*]
  pull_request:
    branches: [main, develop, release/*]

jobs:
  lint:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - run: pip install ruff
      - run: ruff check rph_core/

  gate-a:
    runs-on: self-hosted  # 需要 ORCA/Gaussian 的 runner
    if: contains(github.event.pull_request.changed_files, 'rph_core/steps/conformer_search/')
    steps:
      - uses: actions/checkout@v4
      - run: pytest tests/gate_a/test_ext_equivalence.py -v
```

---

## 8. 文档与配置版本同步

### 8.1 V2.10_Plan.md 版本标记

在文档头部添加版本和状态：

```markdown
# RPH v2.10 S1 四层协议栈修改方案

> **版本**: v2.10.0-rc1  
> **状态**: 🚧 In Progress / ✅ Released / 📋 Planned  
> **最后更新**: 2026-04-08  
> **对应分支**: `release/v2.10`
```

### 8.2 配置版本校验

在代码中添加配置版本检查，防止新旧代码混用：

```python
# rph_core/utils/config_loader.py
CONFIG_VERSION = "2.10.0"

def load_config(path: Path) -> Dict:
    config = yaml.safe_load(path.read_text())
    
    # 检查最小版本要求
    min_version = config.get("_meta", {}).get("min_version", "2.0.0")
    if parse_version(min_version) > parse_version(CONFIG_VERSION):
        raise ConfigError(
            f"Config requires RPH >= {min_version}, "
            f"but current version is {CONFIG_VERSION}"
        )
    
    return config
```

### 8.3 配置变更日志

```yaml
# config/defaults.yaml 顶部注释
# RPH Configuration Schema v2.10.0
# 
# Changelog:
#   v2.10.0:
#     - Add step1.protocol and step1.protocol_stack sections
#     - Add experimental.protocol_stack.enabled flag
#     - Add step2.pes_adapter configuration
#   v2.9.0:
#     - ...
```

---

## 9. 常用 Git 工作流命令

### 9.1 日常开发

```bash
# 开始新功能
git checkout develop
git pull origin develop
git checkout -b feature/lite-zero-protocols

# 日常提交（小而频）
git add pipeline/stages/low_level_rerank.py
git commit -m "feat(stages/lite): implement r2SCAN-3c SP calculation

- Use ORCAInterface for single-point energy
- Support solvent model mapping via solvent_resolver
- Add error handling for ORCA failures

Part of #42"

# 推送并创建 PR
git push origin feature/lite-zero-protocols
gh pr create --title "feat: implement lite and zero protocols" \
             --body-file .github/PULL_REQUEST_TEMPLATE.md \
             --base develop

# 处理 review 意见后更新
git add .
git commit -m "refactor(stages/lite): address review comments

- Extract common SP logic into utils
- Add type hints for all public methods
- Improve error messages for missing input files"
git push origin feature/lite-zero-protocols
```

### 9.2 发布流程

```bash
# 1. 从 develop 创建 release 分支
git checkout develop
git pull origin develop
git checkout -b release/v2.10

# 2. 版本号更新（更新 __version__.py, CHANGELOG.md, V2.10_Plan.md）
vim rph_core/__version__.py
vim CHANGELOG.md
git add -A
git commit -m "chore(release): bump version to v2.10.0-rc1"

# 3. 运行完整验收门
pytest tests/gate_a/ tests/gate_b/ tests/gate_c/ -v

# 4. 合并到 main（通过 PR）
gh pr create --title "release: v2.10.0-rc1" --base main

# 5. 打标签
git checkout main
git pull origin main
git tag -a v2.10.0-rc1 -m "Release candidate 1 for v2.10.0"
git push origin v2.10.0-rc1

# 6. 合并 release 回 develop
git checkout develop
git merge release/v2.10
```

### 9.3 紧急修复

```bash
# 生产环境发现紧急 bug（main 分支）
git checkout main
git checkout -b hotfix/orca-sp-mpi-crash

# 修复并测试
vim rph_core/utils/orca_interface.py
git commit -m "fix(qc/orca): handle MPI abort in SP calculation

Orca may abort after SCF convergence on some systems.
Add retry with single-core fallback.

Fixes #911"

# 快速 review 并合并到 main
gh pr create --title "hotfix: orca SP MPI crash" --base main --label hotfix

# 同步到 develop
git checkout develop
git cherry-pick <hotfix-commit-hash>
```

---

## 10. 科研可重复性最佳实践

### 10.1 计算结果归档

每个重要版本发布后，归档计算结果：

```bash
# 在 release 标签处归档基准结果
git checkout v2.10.0
mkdir -p archive/v2.10.0/
cp -r benchmark/results/v210/* archive/v2.10.0/
git add archive/v2.10.0/
git commit -m "archive(v2.10.0): add benchmark results for v2.10.0 release"
git tag -a v2.10.0-archive -m "Archived benchmark results for v2.10.0"
```

### 10.2 环境复现

```bash
# 生成环境快照
pip freeze > requirements-v2.10.0.txt
conda env export > environment-v2.10.0.yml

# 提交到版本库
git add requirements-v2.10.0.txt environment-v2.10.0.yml
git commit -m "chore(env): add frozen dependencies for v2.10.0

- Python 3.9.18
- ORCA 5.0.4
- xTB 6.6.1
- RDKit 2023.09.1"
```

### 10.3 Paper 附录同步

```bash
# 创建 paper-appendix 分支
git checkout -b paper/nature-2026-rph-v210 main

# 仅保留计算相关文件（使用 git sparse-checkout）
git sparse-checkout init --cone
git sparse-checkout set \
  rph_core/steps/ \
  config/defaults.yaml \
  docs/V2.10_Plan.md \
  tests/gate_a/ \
  tests/gate_b/ \
  tests/gate_c/ \
  benchmark/results/v210/

# 提交到 Zenodo / Figshare
git archive --format=zip -o rph-v210-paper-appendix.zip HEAD
```

---

## 11. 快速参考卡

| 操作 | 命令 |
|------|------|
| 开始功能 | `git checkout -b feature/name develop` |
| 小步提交 | `git commit -m "type(scope): description"` |
| 推送 PR | `git push origin feature/name && gh pr create` |
| 同步 develop | `git checkout feature/name && git rebase develop` |
| 紧急修复 | `git checkout -b hotfix/name main` |
| 打标签 | `git tag -a v2.10.0 -m "Release v2.10.0"` |
| 查看历史 | `git log --oneline --graph --all` |
| 比较分支 | `git diff develop..feature/name` |

---

> **一句话总结**: 用 Git Flow 管理功能开发，用 Feature Flags 控制实验性功能，用三门验收门保护质量，用小而频的 PR 替代巨型合并，用版本标签保证科研可重复性。