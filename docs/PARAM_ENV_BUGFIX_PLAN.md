# RPH 参数传递与环境发现关键 BUG 修复计划

范围：ReactionProfileHunter v3.0.1 / `RPH_V3.0.1_pure`

目标：先修复会导致 QC 参数或可执行程序路径实际失效的关键问题，再删除/合并重复发现逻辑。S4 仍保持外部化，不在本仓库新增 S4 逻辑。

## 当前结论

RPH 的配置主链路是 `defaults.yaml -> load_config() -> resolve_resources() -> normalize_qc_config() -> step engines`。资源参数继承基本成立，但“可执行程序自动发现”目前只完成了日志报告，没有形成稳定的运行时声明；同时 xTB、ORCA、ISOSTAT、Gaussian rescue、Gau_XTB 等路径各自有重复解析逻辑，导致同一配置在不同路径上行为不一致。

优先处理的核心问题不是“发现能力不够智能”，而是发现结果没有被统一消费、统一回写、统一错误处理。

## P0 必修 BUG

### 1. discovery.known_dirs 配置失效

证据：

- `config/defaults.yaml` 声明了 `executables.discovery.known_dirs.*`。
- `rph_core/utils/resource_utils.py:316` 的 `resolve_all_executables()` 实际使用模块常量 `_KNOWN_DIRS`，在 `resource_utils.py:329` 没有读取 YAML 里的 `discovery.known_dirs`。
- `_resolve_mpi()` 同样使用 `_MPI_KNOWN_DIRS`，没有读取 `discovery.known_dirs.mpirun`。

影响：

- 用户修改 `defaults.yaml` 中的 known_dirs 不会影响自动发现。
- 配置文件看似是 single source of truth，但实际有第二份硬编码真相源。

修复计划：

1. 在 `resolve_all_executables(config)` 中读取：
   - `config["executables"]["discovery"]["known_dirs"][prog]`
   - 缺省时才 fallback 到 `_KNOWN_DIRS[prog]`
2. `_resolve_mpi(config, allow_discovery)` 同步支持 `known_dirs.mpirun`。
3. 新增测试：
   - 临时目录创建 fake executable。
   - 设置 `discovery.known_dirs.xtb` 指向该目录。
   - 断言 `resolve_all_executables()` source 为 `known_dirs` 且 path 为 fake executable。

### 2. 自动发现只记录日志，不回写/声明到 config

证据：

- `ReactionProfileHunter.__init__()` 在 `orchestrator.py:172-173` 做资源和 QC 配置规范化，然后在 `orchestrator.py:223` 调用 `_resolve_executables()`。
- `_resolve_executables()` 在 `orchestrator.py:231` 拿到 `resolve_all_executables()` 结果，只打印 `[discovery]` 日志，不把发现结果写回 `self.config["executables"][prog]["path"]`。
- 后续仍有代码直接读取 `config.executables.*.path`：
  - `gau_xtb_interface.py:87-88` 直接读 Gaussian path。
  - `gau_xtb_interface.py:297-300` 直接读 xTB path 并写入 `XTB_PATH`。
  - `conformer_search/engine.py:249-251` 直接读 ISOSTAT/Shermo path。
  - `small_molecule_precompute.py:152-156` 直接读 Shermo path。
  - `qc_interface.py:1158` 的 `XTBInterface.enso_thermo()` 直接读 xTB path。

影响：

- 启动期发现了 `$PATH` 或 env 中的工具，后续直读配置的路径仍可能是 `/opt/software/...` 旧值。
- 日志显示 found，但实际执行路径可能继续失败，属于“发现链”和“执行链”断裂。

修复计划：

1. 新增一个单一入口，例如 `declare_resolved_executables(config, results)`：
   - 仅当 `found=True` 时写回 `executables.<prog>.path = str(path)`。
   - 保留原始显式配置可选记录到 `executables.<prog>.declared_from` 或仅在结果中返回，不建议污染 YAML schema 过多。
   - 对 ORCA MPI，把 `_resolve_mpi()` 已推导的 `mpi_bin_dir/mpi_lib_dir` 保留。
2. `_resolve_executables()` 调用后回写 `self.config`。
3. 新增测试：
   - mock `resolve_all_executables()` 返回 xTB/Shermo found。
   - 初始化或直接调用 `_resolve_executables()` 后断言 config path 已更新。

### 3. 环境变量语义不统一：目录变量被当成可执行文件

证据：

- `resource_utils.py:303-310` 声明了 `ORCA_PATH`, `MOLCLUS_HOME`, `MPI_HOME` 等 env vars。
- `find_executable()` 对 env var 的处理只接受 `Path(env_value).is_file()`。
- 但运行时声明里 `orca_interface.py:860-861` 把 `ORCA_PATH`/`ORCA_DIR` 设置为 ORCA 所在目录。
- `isostat_runner.py:78-81` 把 `MOLCLUS_HOME` 当目录，并拼接 `isostat`。
- `scripts/Gau_XTB/xtb.sh` 支持 `XTB_PATH` 为可执行文件或目录。

影响：

- 用户按常见习惯设置 `ORCA_PATH=/path/to/orca_dir` 时，统一发现层可能判断无效。
- `MOLCLUS_HOME` 在统一发现层无效，但在 ISOSTAT 私有解析层有效，行为不一致。

修复计划：

1. 给 `_ENV_VARS` 增加语义元数据，或在 `find_executable()` 支持 env value 为目录时拼接 `binary_names`：
   - file: `/path/to/orca`
   - dir direct: `/path/to/orca_dir/orca`
   - dir bin: `/path/to/orca_dir/bin/orca`
2. 明确每个 env var 语义：
   - `*_BIN` 优先解释为文件。
   - `*_PATH`, `*_HOME`, `*_DIR` 可解释为目录或文件。
3. 新增测试：
   - `ORCA_PATH=<tmpdir>`，目录下有 `orca`。
   - `MOLCLUS_HOME=<tmpdir>`，目录下有 `isostat`。
   - `XTB_PATH=<tmpdir>`，目录下有 `bin/xtb`。

### 4. Gaussian 执行路径不一致

证据：

- `GaussianInterface.__init__()` 在 `qc_interface.py:1693-1706` 尊重 `use_wrapper`，direct 模式调用 `resolve_executable_config()`。
- `GaussianRunner.run()` 在 `qc_interface.py:929` 硬编码 `["g16", "input.gjf"]`。
- `qst2_rescue.py:211-219` 强制使用 `wrapper_path`，没有尊重 `use_wrapper=false`。
- `irc_driver.py:189-197` 在 direct 模式下也 fallback 到字面量 `g16`，没有使用 `executables.gaussian.path` 或统一解析结果。

影响：

- 主 TS Berny 路径可能能跑，QST2/IRC/rescue 路径却因 wrapper 缺失或 `g16` 不在 PATH 失败。
- `executables.gaussian.path` 与自动发现结果对部分 Gaussian 入口无效。

修复计划：

1. 新增统一 helper：`resolve_gaussian_command(config) -> tuple[list[str] | str, bool]` 或更简单的 `get_gaussian_command(config)`。
2. 替换以下位置：
   - `GaussianRunner.run()`
   - `QST2RescueDriver._submit_gaussian_job()`
   - `IRCDriver` 的 Gaussian submit 路径
   - Gau_XTB 的 Gaussian submit 路径
3. 行为契约：
   - `use_wrapper=true`：要求 wrapper 存在，否则 fail-fast。
   - `use_wrapper=false`：使用 declared/resolved Gaussian path；找不到才 fallback 到 `g16`，且 warning 明确。
4. 新增测试：
   - `use_wrapper=false + path=/tmp/fake/g16` 时，QST2/IRC/GaussianRunner 命令均使用 fake path。
   - `use_wrapper=true + wrapper_path missing` 时初始化或提交前明确报错。

### 5. xTB 参数传递不完整

证据：

- `XTBRunner._verify_executable()` 在 `xtb_runner.py:87-140` 自己实现查找，不使用 `resolve_executable_config()`。
- `XTBRunner.optimize()` 在 `xtb_runner.py:262` 读 `resources.nproc`。
- `XTBRunner.run_scan()` 和 path search 也读 `resources.nproc`。
- `XTBInterface.optimize()` 会把 `self.nproc` 写回临时 config 的 `resources.nproc`，但其他直接构造 `XTBRunner` 或 `XTBInterface.enso_thermo()` 路径仍可能绕过 `step2.xtb_settings.nproc`、`step1.conformer_search.common.threads`。
- `XTBInterface.enso_thermo()` 在 `qc_interface.py:1158` 直接读取 `config.executables.xtb.path`，跳过 fallback。

影响：

- S2/S1/MRRHO 中 xTB 核数可能看起来配置了 step 子树，实际执行使用全局资源或旧 path。
- xTB 自动发现行为在 optimize/scan/enso 之间不一致。

修复计划：

1. `XTBRunner._verify_executable()` 改为调用 `resolve_executable_config(config, "xtb", env_vars=["XTB_PATH", "XTB_BIN"])`。
2. 删除 `XTBRunner.FALLBACK_PATHS` 或仅作为 `_KNOWN_DIRS` 的兼容输入，不再维护第二份 fallback。
3. `XTBInterface.enso_thermo()` 改用统一解析/声明后的 path。
4. 明确资源入口：
   - XTBRunner 只消费 `resources.nproc`。
   - 各上层 interface 负责把 step-level nproc 映射到 `resources.nproc`。
5. 新增测试：
   - `executables.xtb.path` 无效且 `fail_on_invalid_explicit=false` 时，XTBRunner 能走 env/PATH/known_dirs。
   - `XTBInterface(nproc=4).enso_thermo()` 调用 `run_xtb_enso()` 时 nproc 为 4。

## P1 重构与删除冗余代码

### A. 删除/合并重复 executable lookup

建议合并到 `resource_utils.resolve_executable_config()`：

- `xtb_runner.py:52` 的 `FALLBACK_PATHS`
- `xtb_runner.py:87-140` 的 `_verify_executable()` 私有查找链
- `orca_interface.py:725-762` 的 `_find_orca_binary()` 中二次 `shutil.which("orca")` fallback
- `isostat_runner.py:64-105` 的 `_iter_isostat_candidates()` / `_resolve_isostat_path()`
- `qc_interface.py:1158` 的 xTB path 直读
- `gau_xtb_interface.py:87-88`、`297-300` 的 Gaussian/xTB path 直读

保留策略：

- 每个 runner 可以保留一个薄 wrapper，例如 `_resolve_xtb_path()`，但内部只能调用统一解析函数。
- 工具特有的运行时环境构造可以保留在各自模块，例如 ORCA 的 MPI PATH/LD_LIBRARY_PATH、Gaussian 的 GAUSS_SCRDIR。

### B. 删除或降级静态 Gaussian 模板

证据：

- `config/AGENTS.md` 声明 `templates/*.gjf/*.com` 要运行时读取。
- 现有代码实际由 `InputFactory.create()`、`GaussianInterface.write_input_file()`、`qst2_rescue.py` 等内联组装输入。
- `config/templates/` 中同时存在 Jinja 风格 `{{mem}}`、Python format 风格 `{mem}` 和静态占位 `[coordinates]`，格式不统一。

二选一决策：

1. 推荐短期方案：把模板降级为 `docs/examples/gaussian_inputs/` 示例文件，并修改 `config/AGENTS.md`，删除“运行时读取”声明。
2. 中期方案：真正引入模板渲染器，并将 Gaussian 输入生成全部收敛到模板系统。

当前建议：

- 不在本轮引入模板引擎，避免扩大行为面。
- 先删除“运行时模板”承诺，或把 `config/templates/` 移到 docs 示例目录。

### C. 配置副本清理

当前 `config/` 下存在：

- `defaults.yaml`
- `defaults.yaml.bak`
- `defaults_alcl3_test.yaml`
- `defaults_la_test.yaml`
- `defaults_la_test_dataset.yaml`
- `defaults_la_test_fast.yaml`
- `defaults_zncl2_as_mgcl2_remaining.yaml`

问题：

- `config/AGENTS.md` 明确说 `defaults.yaml` 是唯一配置源，避免 forks。
- `docs/RECOMPUTE_PLAN_V3.0.1.md` 仍引用部分 `defaults_la_test*.yaml`。

清理计划：

1. 先检查这些文件是否被测试、脚本或文档实际调用。
2. 若只是实验运行配置：
   - 移到 `docs/config_examples/` 或 `config/examples/`。
   - 主运行只保留 `defaults.yaml`。
3. 如果仍要支持实验配置，不要命名为 `defaults_*`，改为 `example_*.yaml`，并在 README 中明确“不作为行为源”。
4. 删除 `defaults.yaml.bak`，除非有明确回滚用途；回滚应靠版本控制，不靠 bak。

## P2 质量门与 schema

### 1. 最小 schema 校验

不建议立刻引入重型依赖。先在 `config_loader.py` 后增加轻量校验：

- 顶层必须是 dict。
- 必须有 `executables`, `resources`, `theory`, `step1`, `step2`, `step3`, `run`。
- `executables.discovery.known_dirs` 必须是 dict。
- `resources.nproc` 为 null 或正整数。
- `resources.mem` 为 null 或可由 `mem_to_mb()` 解析。

### 2. discovery report 持久化

建议将启动期发现结果写到运行目录：

- `pipeline.state` 或 `run_manifest.json`
- 字段：program, found, source, path, required, declared_path

这样外部排查时不用只看滚动日志。

## 推荐实现顺序

1. 修 `resource_utils.resolve_all_executables()` 读取 YAML `known_dirs`，补测试。
2. 实现发现结果回写/声明到 runtime config，补测试。
3. 修 env var 目录/文件双语义，补 ORCA_PATH/MOLCLUS_HOME/XTB_PATH 测试。
4. 收敛 xTBRunner 和 `XTBInterface.enso_thermo()` 到统一解析。
5. 收敛 Gaussian submit helper，覆盖 GaussianRunner/QST2/IRC/Gau_XTB。
6. 再处理 ISOSTAT/ORCA 私有 fallback 的删除。
7. 最后处理模板和配置副本的文档/目录清理。

## 验收测试清单

必须运行：

```bash
pytest tests/test_executable_resolution.py -v
pytest tests/test_gaussian_env.py -v
pytest tests/test_orca_interface.py -v
python scripts/ci/check_imports.py rph_core
```

建议新增：

```bash
pytest tests/test_executable_resolution.py::test_known_dirs_from_yaml_override -v
pytest tests/test_executable_resolution.py::test_discovery_declares_runtime_paths -v
pytest tests/test_executable_resolution.py::test_env_dir_semantics_for_orca_xtb_isostat -v
```

## 不纳入本轮

- 不新增 S4 feature extraction。
- 不重写 QC 执行框架。
- 不引入递归文件系统扫描、module/lmod、Slurm 自动探测。
- 不把所有 direct subprocess 一次性重构掉；先处理会影响参数/path 传递的 QC 主路径。

## 风险点

- 现有默认 `executables.*.path` 多为 `/opt/software/...`，且 `fail_on_invalid_explicit=true`。如果用户机器没有这些路径，自动发现会在显式路径无效时直接 fail-fast。需要决定默认配置是否应把 path 置空，让 env/PATH 生效。
- 在 Windows/WSL 混合路径下，`Path.exists()` 对 Linux 路径可能在 Windows Python 中失败。真正跑 QC 的环境应以 Linux/WSL 为准，Windows 侧只做代码和 mock 测试。
- ORCA 的 `ORCA_PATH` 在不同生态里既可能表示目录也可能表示 binary；本轮应兼容两者，不强行改变用户习惯。
