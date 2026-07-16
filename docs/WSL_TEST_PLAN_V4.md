# RPH 4.0 WSL 测试方案

状态：执行前基线审计  
日期：2026-07-11  
适用版本：RPH V4 S0–S4

## 1. 当前开发状态

| 模块 | 当前状态 | 测试结论 |
|---|---|---|
| V4 入口 | `rph_core.v4_orchestrator`、`bin/rph_run`、`python -m rph_core` 已切换 | 已通过编译检查 |
| S0 | 由 `SMARTSMatcher` 生成 `S0_Mechanism/mechanism.json` | 已实现，但当前运行顺序是在 S1 生成几何后执行 |
| S1 | V4 强制 `censo_lite`；CREST/GFN2、B97-3c SP、xTB mRRHO、二面角签名去重 | 已切换到独立 `CensoLiteRuntime`，需真实 QC 实测 |
| S2 | `PEBScanner` 目前是旧 `RetroScanner` 的 V4 别名 | 功能可接入，破坏性拆分尚未完成 |
| S3 | 新 `LowLevelEngine` + `StageCalculator`；ORCA B97-3c OPT/OptTS + r2SCAN-3c SP | 静态检查通过，尚无真实 ORCA 运行结果 |
| S4 | 新 `HighLevelEngine`；ORCA M062X OPT/OptTS/FREQ + ORCA wB97M-V SP | 静态检查通过，需真实 ORCA 运行验证 |
| 检查点 | `pipeline.state` 记录签名并参与 S0–S4 manifest 复用 | 已实现，需做中断/配置变更实测 |
| QCTaskRunner | V4 活动路径和旧 V3 文件已移除 | 需通过全量 import 和 WSL 运行确认无残留 |

## 2. 已完成的本地验证

在 Windows 工作区中已验证：

```text
Python 3.12.7
py_compile                  PASS
scripts/ci/check_imports.py PASS（142 个 Python 文件）
V4/检查点/路径相关测试       PASS（10 passed）
```

完整 pytest 当前不能作为通过依据，原因是当前 Windows Python 环境缺少：

- `rdkit`
- `rph_features`（外部 RPH_Postprocess 测试依赖）

WSL 当前状态为 Ubuntu/WSL2，32 核、约 54 GiB 可用内存；仓库已挂载到 `/mnt/e`，但尚未发现 `pytest`、`rdkit`、ORCA、CREST 或 xTB 可执行文件。

## 3. WSL 环境准备

### 3.1 进入仓库和创建 Python 环境

```bash
export RPH_ROOT='/mnt/e/Calculations/AI4S_ML_Studys/[4+3] Mechain learning/ReactionProfileHunter/RPH_V4.0.0'
cd "$RPH_ROOT"

sudo apt update
sudo apt install -y python3-venv python3-dev build-essential cmake gfortran \
  openmpi-bin libopenmpi-dev libgomp1

python3 -m venv ~/.venvs/rph4
source ~/.venvs/rph4/bin/activate
python -m pip install --upgrade pip wheel
python -m pip install -r requirements.txt
python -m pip install pytest
```

若 `rdkit` 的 pip wheel 在目标 Ubuntu/Python 组合上不可用，使用一个带 RDKit 的 Conda/Mamba 环境替代虚拟环境；不要在同一个环境中混装多个 RDKit 来源。

### 3.2 配置 QC 软件

V4 的默认配置期望以下路径，建议先保持路径一致：

```text
/opt/software/orca/orca
/opt/software/crest/crest
/opt/software/xtb/bin/xtb
```

ORCA、CREST 和 xTB 必须分别确认：

```bash
test -x /opt/software/orca/orca
test -x /opt/software/crest/crest
test -x /opt/software/xtb/bin/xtb
```

ORCA 需要其动态库目录可被找到。若实际安装路径不同，应复制 `config/defaults.yaml` 为本地配置并只修改 `executables.*.path`、库路径和资源配置，不要修改代码中的路径。

### 3.3 运行时依赖检查

```bash
python - <<'PY'
import importlib
for name in ('numpy', 'scipy', 'yaml', 'pytest', 'rdkit'):
    module = importlib.import_module(name)
    print(name, getattr(module, '__version__', 'ok'))
PY

for exe in orca crest xtb; do command -v "$exe" || true; done
nproc
free -h
df -h /tmp "$RPH_ROOT"
```

## 4. 分阶段测试顺序

### T0：配置和静态门禁（无 QC）

```bash
cd "$RPH_ROOT"
python - <<'PY'
from pathlib import Path
import yaml
cfg = yaml.safe_load(Path('config/defaults.yaml').read_text())
assert cfg['step1']['protocol'] == 'censo_lite'
assert cfg['step1']['allowed_protocols'] == ['censo_lite']
assert cfg['theory']['s3_low_level']['optimization']['engine'] == 'orca'
assert cfg['theory']['s3_low_level']['optimization']['method'] == 'B97-3c'
assert cfg['theory']['s3_low_level']['single_point']['method'] == 'r2SCAN-3c'
assert cfg['theory']['s4_high_precision']['single_point']['method'] == 'wB97M-V'
print('V4_CONFIG_OK')
PY

python -m py_compile \
  rph_core/v4_orchestrator.py \
  rph_core/steps/stage_calculator.py \
  rph_core/steps/conformer_search/censo_lite.py \
  rph_core/steps/conformer_search/deduplicator.py \
  rph_core/steps/conformer_search/torsion_signature.py \
  rph_core/steps/step2_retro/peb_scanner.py \
  rph_core/steps/step3_lowlevel/engine.py \
  rph_core/steps/step4_highlevel/engine.py \
  rph_core/utils/qc_jobs.py \
  rph_core/utils/qc_models.py \
  rph_core/utils/v4_checkpoint.py

python scripts/ci/check_imports.py rph_core
pytest -q \
  tests/test_v4_protocol_contract.py \
  tests/test_v4_checkpoint.py \
  tests/test_sandbox_toxic_paths.py \
  tests/test_checkpoint_partial_resume.py
```

T0 通过标准：编译成功、导入门禁成功、V4 相关测试全部通过。

### T1：QC 接口 smoke test

先分别用一个小的闭壳层分子验证 ORCA 的输入生成和执行，再进入 RPH。至少验证：

1. ORCA `B97-3c Opt`；
2. ORCA `B97-3c OptTS` 输入能够生成；
3. ORCA `r2SCAN-3c SP`；
4. ORCA `wB97M-V/def2-TZVPP` SP；
5. ORCA `M062X/def2-SVP Opt`；
6. ORCA `M062X/def2-SVP OptTS` 与 `Freq` 输入和运行；
7. ORCA 输出能被现有接口解析出能量、坐标和频率。

每个 smoke job 必须保留输入、输出、退出码和运行时间。T1 不要求找到真实过渡态，但要求 `OptTS` 作业能启动并给出可诊断结果。

### T2：S1 CENSO-LITE 单阶段测试

使用一个已知可以被当前 SMARTS 注册表识别的最小反应样本，不要先使用大规模数据集：

```bash
rm -rf /tmp/rph4_s1_smoke
bin/rph_run --smiles '<已登记的最小反应产物 SMILES>' \
  --reaction-type '[4+3]_default' \
  --output /tmp/rph4_s1_smoke
```

若完整入口会继续执行 S2–S4，可先用 Python 直接调用 `CensoLiteEngine`，或为 CLI 增加明确的 `--stop-after s1` 测试开关。验收文件：

```text
S1_ConfSearch/product/manifest.json
S1_ConfSearch/product/candidates/conf_*.xyz
S1_ConfSearch/product/crest/crest_conformers.xyz
```

验收内容：

- manifest 中 `protocol == censo_lite`；
- 没有高精度 OPT/FREQ 输出；
- 候选保留窗口和二面角签名去重生效；
- 镜像对映体和立体化学保护策略符合配置；
- B97-3c 仅作为 S1 排序 SP，不作为 S1 几何优化；
- 失败候选仍能留下原始 XYZ、错误信息和可追溯来源。

### T3：S0/S2 机制和 PEB 测试

检查：

```text
S0_Mechanism/mechanism.json
S2_PEB/manifest.json
S2_PEB/peb_profile.*
S2_PEB/ts_guess.xyz 或 peb_peak 对应文件
S2_PEB/intermediate.xyz
```

验收内容：

- `forming_bonds` 使用 0-based 索引且与 XYZ 原子数一致；
- S0 和 S2 manifest 的键元数据一致；
- PEB 失败时状态为 degraded/failed，并保留扫描曲线，而非静默丢失；
- TS seed 与 intermediate seed 可以独立被 S3 读取；
- 不再生成旧的 S2 高级别预优化产物。

当前特别需要记录：代码实际先运行 S1 生成几何，再执行 S0 SMARTS 识别。因此若项目要求严格的 S0→S1 顺序，必须把 S0 改为“SMILES/拓扑验证”，再把几何相关验证作为 S1 后的 S0 几何补充，不能仅靠测试文档解释。

### T4：S3 低级别批量测试

使用 S1 的 1 个 minimum、S2 的 intermediate、S2 的 TS seed，验证：

```text
S3_LowLevel/<id>/opt/
S3_LowLevel/<id>/sp/
S3_LowLevel/manifest.json
```

验收矩阵：

| 结构 | OPT | SP | 允许结果 |
|---|---|---|---|
| minimum | ORCA B97-3c Opt | ORCA r2SCAN-3c SP | complete 或 opt_failed_sp_complete |
| intermediate | ORCA B97-3c Opt | ORCA r2SCAN-3c SP | complete 或 opt_failed_sp_complete |
| TS seed | ORCA B97-3c OptTS | ORCA r2SCAN-3c SP | complete 或 opt_failed_sp_complete |

必须确认：

- S3 运行日志中没有 Gaussian 调用；
- S3 运行日志中没有 `QCTaskRunner`、FREQ、NBO 或旧 TS rescue；
- OPT 失败时仍尝试对原始输入进行 SP；
- OPT 和 SP 均成功时 `usable_for_ml == true`；
- 每个结构的输入、优化坐标、SP 输入和输出路径都写入 manifest。

### T5：S4 高精度测试

S4 必须接收 S3 的优化坐标；S3 OPT 失败时才允许回退到 S3 输入坐标：

```text
S4_HighLevel/<id>/opt/
S4_HighLevel/<id>/sp/
S4_HighLevel/manifest.json
```

验收矩阵：

| 结构 | OPT | SP |
|---|---|---|
| minimum/intermediate | ORCA M062X/def2-SVP Opt + Freq 验证 | ORCA wB97M-V/def2-TZVPP + def2/J SP |
| TS seed | ORCA M062X/def2-SVP OptTS + Freq 验证 | ORCA wB97M-V/def2-TZVPP + def2/J SP |

必须确认：

- S4 对每一个 S3 结构执行完整 OPT+SP，不只计算最低能结构；
- ORCA 的优化、频率和单点输出文件均被保留；
- ORCA OPT 不收敛时，manifest 仍保留 S3 坐标，将 S4 标记为 incomplete，并令该结构 `usable_for_ml == false`；
- S3 低精度结果不会被 S4 失败覆盖；
- ML 使用的最低保障字段仍可从 S3 manifest 读取。

### T6：端到端回归和失败注入

至少完成以下失败注入：

1. 删除 ORCA：S3/S4 manifest 应报告失败，不应 Python 崩溃；
2. 将 S4 engine 配成非 ORCA：配置应被明确拒绝，S3 结果仍完整；
3. 让一个结构 OPT 超时：其他结构继续运行；
4. 删除一个 XYZ：对应结构 degraded，其他结构继续运行；
5. 修改配置后重跑：签名变化，不能错误复用旧 manifest；
6. 传入旧 `pipeline.state`：V4 明确拒绝，不得静默兼容旧状态。

## 5. 最终验收命令和门槛

```bash
cd "$RPH_ROOT"
pytest -q
bin/rph_run --smiles '<已登记的最小反应产物 SMILES>' \
  --reaction-type '[4+3]_default' \
  --output /tmp/rph4_acceptance
```

V4 进入 benchmark 前必须同时满足：

1. T0、T1、T2、T3、T4、T5 全部有日志和 manifest；
2. 至少一个真实样本完成 S3 的 OPT+SP；
3. 至少一个真实样本完成 S4 的 ORCA OPT+SP；
4. S3/S4 不再调用 Gaussian、NBO 和 QCTaskRunner；
5. S4 失败时 S3 低级别坐标和能量仍可直接用于 ML；
6. 输出目录中不得混入旧版 `S1_ConfGeneration`、`S2_Retro`、`S3_TransitionAnalysis` 作为 V4 的主产物；
7. 全量测试中的外部 `rph_features` 依赖必须单独安装或从 RPH 核心测试集合中隔离，不能把 collection error 当作 V4 通过。

## 6. 测试期间需要继续修复的代码项

按优先级排列：

1. 删除 `config/defaults.yaml` 中仍残留的 `ext/full/lite/zero` 协议栈，确保只有 `censo_lite`；
2. 把 `CensoLiteEngine` 从旧 `ConformerEngine` 的高精度 DFT 代码中彻底解耦，移除旧 `QCTaskRunner` 导入链；
3. 将 `PEBScanner` 从旧 `RetroScanner` 别名改为独立 V4 PEB 实现或明确的最小适配层；
4. 让 `V4Checkpoint.reusable()` 真正参与 orchestrator 的阶段短路和 manifest 复用；
5. 更新 README/README.zh-CN 的流程和输出合同，当前文档仍描述旧 S0–S3 及外部 S4；
6. 为 `StageCalculator` 增加 mock QC 单元测试，覆盖 minimum、TS、OPT 失败但 SP 成功、OPT/SP 均失败四种路径。

完成上述六项后，再进行大规模 [4+3] benchmark；在此之前只做小样本和单结构真实 QC，避免把旧架构残留与量化化学失败混在一起。
