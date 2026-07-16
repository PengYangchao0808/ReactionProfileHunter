# RPH V4 S1 构象搜索阶段问题汇总

**日期：** 2026-07-14
**范围：** V4 S1 `CENSO-LITE` 构象搜索阶段
**依据：** 当前实现代码、`config/defaults.yaml`、`rx_id_1` precursor 实际运行日志与输出目录
**文档性质：** 问题盘点与影响分析，不包含修复方案

## 1. 执行流程概况

S1 对 precursor 及每个 product variant 分别执行：

```text
RDKit ETKDGv3 -> CREST GFN2/ALPB -> 扭转角/重原子 RMSD 去重
  -> xTB 相对能窗口预筛 -> 多样性截断 -> ORCA B97-3c/CPCM SP
  -> B97-3c 相对能窗口 -> thermo pool 截断
  -> xTB GFN1 SPH+mRRHO -> E(B97-3c)+G(RRHO) 最终排序
  -> 固定窗口/min_keep/max_keep -> manifest.json + selected.xyz
```

S1 不执行 DFT 几何优化或频率计算，符合 V4 理论合同。与完整 CENSO 分阶段优化/精修协议不同。

## 2. 本次实际运行暴露的问题

### 2.1 B97-3c SP 候选数仍然较大

`rx_id_1` precursor 实际数量：

```text
CREST 原始构象：191
B97-3c SP：60（达 prefilter.max_candidates_soft 上限）
mRRHO pool：约 30
```

B97-3c 前已执行 RMSD 去重、`xtb_window_kcal: 5.0`、`max_candidates_soft: 60` 截断。较大的 SP 数量是阈值+多样性保留策略的共同结果，非缺漏 preselection。

### 2.2 B97-3c 构象级并行默认未启用

`ThreadPoolExecutor` 已实现，但默认 `parallel_jobs: 1, cores_per_job: 16`。构象间串行，单任务 16 核。~20s/SP，60 任务约 20 min wall time。属默认资源策略问题。

### 2.3 xTB mRRHO 多线程数值 Hessian 不稳定

28 个 mRRHO 目录中：成功 23，失败 4，未执行 1。

| 构象 | 错误 |
|---|---|
| `conf_0004` | rc=-11，Biased Numerical Hessian 崩溃 |
| `conf_0018` | Fortran 153：allocatable array 未分配 |
| `conf_0025` | Fortran 174：SIGSEGV |
| `conf_0009` | Fortran 153 + Coulomb evaluator setup failure |

所有失败均发生在 `--bhess normal` 阶段。隔离复现：保持 GFN1/ALPB(acetone)/charge/UHF/`--bhess`/`--enso` 不变，线程数从 4 降为 1 后全部成功。

根因分析：
```text
xTB 6.7.1
stack limit: 8192 KB
OMP_STACKSIZE: 未设置
可用内存: ~53 GiB
```

`omp_stacksize: null` 导致代码主动移除继承的 `OMP_STACKSIZE`。设置足够的 OpenMP 线程栈后 4 线程任务成功。本次问题的直接原因是 xTB 6.7.1 并行 Hessian 线程栈需求与运行环境不匹配，非内存耗尽。

### 2.4 GFN0 fallback 因参数搜索路径缺失全部失败

mRRHO fallback（GFN1/4t → GFN0/1t）均报错：

```text
Parameter file param_gfn0-xtb.txt not found
```

文件实际存在于 `/opt/software/xtb/share/xtb/`，但子进程 `XTBPATH` 和 `XTBHOME` 均为空。属环境传播问题，非 GFN0 未安装。

### 2.5 Fallback 同时改变线程数和理论方法

```text
primary:   GFN1 / 4 threads
fallback:  GFN0 / 1 thread
```

一次重试改变两个变量，无法判断成功来自线程还是方法变化。隔离测试已证明 GFN1/1t 即可成功。若部分构象 GFN1 成功、部分 GFN0 fallback 成功，ensemble 会静默混用不同 GFN level 的 `G(RRHO)`，且无 provenance 记录。

### 2.6 重试结构无法表达"同方法降线程重试"

代码仅当 `fallback_gfn_level != primary` 时添加 fallback attempt。无法表达：

```text
GFN1 / 多线程失败 -> GFN1 / 单线程重试
```

重试机制与本次真实故障模式不匹配。

### 2.7 单 mRRHO 失败使整个 ensemble 放弃热校正

默认 `failure_policy: ensemble_electronic_only`。任一构象 mRRHO 最终失败 → 所有候选退回 B97-3c 电子能排序。已成功的 23 个 mRRHO 不参与最终排序，计算时间浪费。

## 3. 科学语义与结果一致性问题

### 3.1 电荷/自旋默认固定为中性闭壳层

S1 从局部配置读取 charge/multiplicity，未配置时默认 charge=0, mult=1。未从 S0 分子记录自动传递。对离子、开壳层体系造成直接理论输入错误。

### 3.2 mRRHO fallback 后 ranking formula 可能与实际 score 不一致

ensemble electronic-only fallback 后 score 重置为 B97-3c 电子能，但 manifest 的 `ranking_formula` 仅由 `xtb_thermo.enabled` 决定，可能仍写为 `E_B97-3c_SP + G(RRHO)_xTB`。成功构象 metadata 可能仍保留 `g_rrho_correction_hartree`，但 `s1_score_hartree` 已不包含。

### 3.3 mRRHO 方法 provenance 缺失

`run_mrrho()` 只返回浮点校正值或 `None`，未返回：GFN level、线程数、是否 fallback、primary/fallback failure 原因、xTB 版本。manifest 无法判断热化学方法来源。

### 3.4 缺少 Boltzmann population

最终筛选仅使用固定 `final_window_kcal`、`min_keep`、`max_keep`。未计算 Boltzmann weight、normalized population、cumulative population。无法按累计 population 筛选。

### 3.5 RPH `censo_lite` 与论文 `CENSO-light` 名称相近但协议不等同

2024 年 RTCONF55-16K 论文正式定义了 `CENSO-zero` / `CENSO-light`。RPH `censo_lite` 只实现了部分低层级 ensemble ranking 思想，未实现论文定义的完整输出和 representative refinement。详见第 11 节。

### 3.6 固定候选数上限可能覆盖能量分布

`max_candidates_soft` 和 `thermo_pool_max` 是固定数量上限，非纯能量窗口。多样性截断优先覆盖不同 torsion signature，cluster 数超上限时删除部分 cluster。

## 4. 可恢复性与缓存问题

### 4.1 SP 内存缓存对 S1 实际无效

`ORCAInterface` SP cache 是实例级、进程内缓存。每个构象新建实例 → 日志 `cache hit rate: 0/1`。不能跨构象、跨进程、跨重启复用。

### 4.2 缺少 B97-3c 逐构象恢复

S1 仅整分子级 checkpoint。B97-3c 阶段中断后：CREST 可复用，但已有 ranking 目录不会被可靠解析，重启可能重算。

### 4.3 缺少 mRRHO 逐构象恢复

同理，已生成的 `xtb_enso.json` 无独立 checkpoint，已完成/失败构象无稳定状态索引。重启可能重算。

### 4.4 Manifest 只在末尾写入

候选状态只在整个 S1 结束后写入。中途无结构化候选文件，中断时无法区分成功/失败/未运行/已过滤。

## 5. 日志与审计问题

### 5.1 缺少实时相对能与排序日志

B97 日志输出绝对 Hartree 能量，无完整 `conf_id / E_B97-3c / relative energy / G(RRHO) / final score / rank` 摘要。

### 5.2 中间筛选数量没有及时输出

`pipeline_log` 内部维护候选中继数量，但运行中无稳定漏斗统计。附件日志无法直接看出 191→60 的筛选过程。

### 5.3 被过滤/失败候选缺少完整历史

manifest 只保留最终候选。未保留：去重/窗口/截断排除的构象、全部 B97-3c 能量、mRRHO 状态、淘汰步骤和原因。

### 5.4 mRRHO 错误日志冗长重复

primary 崩溃 + fallback 崩溃输出两套长日志，关键错误信息淹没在重复信息中。

### 5.5 日志截面易导致状态误判

附件在 B97/mRRHO 中间截取，曾出现判断："只执行了 45 个 B97"、"mRRHO 未实现"、"缺少 xTB preselection"。实际分别为 60 个 SP 任务和已启动 mRRHO。

## 6. 配置契约问题

### 6.1 未实际消费的配置字段

```yaml
prefilter.preserve_torsion_clusters
prefilter.min_per_cluster
deduplication.preserve_stereochemistry
deduplication.preserve_mirror_pairs
retention.energy_window_kcal
retention.always_keep_lowest
```

修改它们不会改变当前行为。

### 6.2 `torsion_bin_deg` 与实际去重判据不同

实际 duplicate 判定使用 `torsion_rmsd_deg` + 重原子 RMSD。`torsion_bin_deg` 影响 signature key 和 diversity grouping，非直接 duplicate 判据。

### 6.3 mRRHO 环境依赖未纳入配置合同

executable 路径已配置，但启动前未验证：GFN parameter file 可发现性、`XTBPATH` 有效性、OpenMP stack 是否满足、xTB build 与 fallback level 兼容性。

## 7. 实现与仓库约束问题

### 7.1 Stage 模块直接调用 `subprocess.run`

`xtb_thermo.py` 直接使用 `subprocess.run`，不符合"所有外部 QC job 通过 `qc_jobs.py`"的编码规则。

### 7.2 外部任务路由不统一

B97-3c（ORCAInterface）与 mRRHO（XTBInterface + `run_xtb_enso()` + 子进程）在缓存、重试、环境变量、checkpoint、失败记录、输出复用上使用不同实现。

### 7.3 Precursor + product variants 各自重复完整 S1

V4 orchestrator 先处理 precursor，再遍历所有 product variants。每个 variant 独立运行 CREST + B97-3c + mRRHO 全流程。

## 8. 测试覆盖问题

### 8.1 缺少真实 xTB mRRHO 集成测试

现有测试覆盖：`xtb_enso.json` 解析、命令参数、模拟 crash fallback、部分能量解析。未覆盖：多线程 Hessian 栈需求、GFN 参数搜索路径、GFN1/GFN0 一致性、部分 mRRHO 失败的 score、fallback 后 manifest formula 真实性。

### 8.2 缺少自由能组合回归测试

无回归断言验证 `s1_score = b973c_sp + g_rrho_correction`。无 ensemble electronic-only fallback 后 score/metadata/formula 同步变化测试。

### 8.3 缺少中断恢复测试

未覆盖：CREST 完成后、部分 B97-3c 后、部分 mRRHO 后、primary 失败 fallback 成功后、manifest 写入前的中断恢复。

### 8.4 缺少参考 CENSO 数值 benchmark

无固定版本/分子集的对照：最低构象一致性、前若干构象 RMSD、相对能排序、population、筛选召回率、QC 任务数和 wall time。

## 9. 问题优先级总览

优先级仅表示影响严重程度，不代表修复顺序。

### P0：可能影响科学正确性或最终结果语义

1. 电荷/自旋未按分子自动传递，默认中性闭壳层。
2. mRRHO fallback 可能混用 GFN1/GFN0 校正。
3. Ensemble fallback 后 manifest ranking formula 与实际 score 不一致。
4. mRRHO 方法和 fallback provenance 未写入结果。

### P1：已造成任务失败、计算浪费或不可恢复

1. xTB 6.7.1 多线程 `--bhess` 在当前 OpenMP 栈环境下不稳定。
2. GFN0 fallback 因 `XTBPATH` 缺失全部失败。
3. 重试结构不能表达同 GFN level 降线程重试。
4. 单 mRRHO 失败使整个 ensemble 放弃已成功热校正。
5. B97-3c 和 mRRHO 缺少逐构象 checkpoint 与磁盘复用。
6. ORCA SP instance cache 对 S1 实际无效。

### P2：效率、可观测性和协议清晰度问题

1. B97-3c 默认构象级串行。
2. 默认 SP pool 上限仍可产生 60 个 B97-3c 任务。
3. 缺少实时相对能、最终 rank 和 Boltzmann population。
4. 被筛除/失败候选缺少完整历史。
5. 多项配置字段未被实际使用。
6. `censo_lite` 名称易与 RTCONF55-16K 论文 `CENSO-light` 混淆。
7. Precursor + product variants 重复完整 S1 成本。
8. mRRHO 错误日志冗长，缺少聚合摘要。

### P3：实现一致性与维护问题

1. Stage 模块直接 `subprocess.run`，违反 QC job 路由规则。
2. ORCA SP 与 xTB mRRHO 使用不同任务生命周期/缓存/恢复语义。
3. 真实二进制、失败回退、分数一致性、中断恢复测试不足。

## 10. 总结

S1 已具备完整基本漏斗：CREST/GFN2、xTB 预筛、B97-3c 排序、xTB mRRHO 自由能校正。问题集中在四类：

1. **运行稳定性：** xTB 多线程 Hessian 栈配置不当 + fallback 环境不可用。
2. **结果一致性：** charge/spin、混合 GFN fallback、ensemble fallback、manifest formula 之间的科学语义风险。
3. **恢复与审计：** 缺少逐构象 checkpoint、磁盘复用和完整候选历史。
4. **效率与可观测性：** B97 候选池偏大、默认串行、缺乏实时漏斗和排序摘要。

这些共同导致本次运行 B97-3c 数量大、mRRHO 大量报错、fallback 全部失效、已成功热化学计算可能被整体弃用。

## 11. RPH `censo_lite` 与 RTCONF55-16K 论文 `CENSO-light` 的核心差异

### 11.1 对比依据

本节"论文 CENSO-light"指以下文章定义的协议：

> B. B. Mészáros et al., *J. Chem. Theory Comput.* **2024**, 20, 7385–7392. DOI: [10.1021/acs.jctc.4c00565](https://doi.org/10.1021/acs.jctc.4c00565).

论文定义四种协议：

| 协议 | Ensemble optimization | Ensemble ranking | Representative refinement |
|---|---|---|---|
| CENSO-zero | xTB | xTB | RSH//GGA |
| CENSO-light | xTB | GGA | RSH//GGA |
| CENSO-default | narrowed GGA ensemble | narrowed RSH ensemble | RSH//GGA |
| CENSO-brute-force | GGA | RSH | RSH//GGA |

论文 CENSO-light = 低成本 ensemble 热力学 + GGA 单点排序 + 单代表构象 RSH//GGA refinement。RPH `censo_lite` 与论文名称相近但非同一协议。

### 11.2 论文 CENSO-light 的核心自由能定义

论文将一个状态的构象 ensemble Gibbs free energy 分解为：

$$G_{\mathrm{ens}} = G_0 + G_{\mathrm{conf}}^{\mathrm{rel}}$$

其中 $G_0$ 为最低自由能构象的精修自由能，$G_{\mathrm{conf}}^{\mathrm{rel}} = -RT\ln\sum_i g_i\exp[-(G_i-G_0)/RT]$ 为其余构象的构象自由能稳定化。

低层级排序自由能：$G_i^{\mathrm{light}} = E_i^{\mathrm{B97-3c//xTB}} + \Delta G_{i,\mathrm{therm}}^{\mathrm{xTB//xTB}}$（热统计校正，非 xTB 总自由能）。

最终语义：
$$G_{\mathrm{ens}}^{\mathrm{CENSO-light}} = G_0^{\mathrm{RSH//GGA}} + G_{\mathrm{conf}}^{\mathrm{rel,GGA//xTB}}$$

### 11.3 论文实际计算层级

| 环节 | 方法 |
|---|---|
| 构象搜索 | CREST GFN2-xTB, iMTD-sMTD, 6.0 kcal/mol window |
| xTB 溶剂 | ALPB(DCM) |
| Ensemble geometry | xTB optimized |
| GGA ranking | B97-3c/CPCM(DCM) SP |
| Thermal contribution | xTB//xTB, 全构象 |
| Repr. geometry refinement | B97-3c/CPCM(DCM) OPT |
| Repr. final SP | ωB97M-V/def2-TZVPP/CPCM(DCM) |

### 11.4 RPH 重构前实际语义

RPH S1 低层级排序：$S_i = E_i^{\mathrm{B97-3c//CREST}} + G_i^{\mathrm{RRHO,GFN1}}$（`G(T)` 作为 correction 使用，避免了重复 xTB electronic energy）。

但 RPH 只输出 rank-1 构象及其排序分数，不组装 $G_{\mathrm{conf}}^{\mathrm{rel}}$，也不将 conformational correction 传递到后续阶段。

### 11.5 核心差异对照表

| 比较项 | 论文 CENSO-light | RPH censo_lite | 差异性质 |
|---|---|---|---|
| 目标输出 | 状态的 ensemble Gibbs free energy | S1 rank-1 构象、有限候选及排序分数 | 输出物理量不同 |
| Ensemble 保留到最终热力学 | 是，通过 $G_{\mathrm{conf}}^{\mathrm{rel}}$ | 否，主要用于选 rank-1 | 核心语义缺失 |
| Partition function | 计算 $Z_{\mathrm{rel}}$ | 未计算 | 核心差异 |
| Boltzmann population | 参与 conformational correction | 未输出 | 核心差异 |
| CREST workflow | GFN2-xTB iMTD-sMTD | GFN2 `imtd_gc` | 搜索协议不同 |
| CREST window | 6.0 kcal/mol | 6.0 kcal/mol | 一致 |
| 溶剂 | ALPB/CPCM(DCM) | ALPB/CPCM(acetone) | benchmark 条件不同 |
| B97 SP 覆盖范围 | 全 xTB ensemble | 去重 + 5 kcal + cap ≤ 60 | RPH 提前截断 |
| Thermal correction 范围 | 全 partition function 构象 | B97 window + cap ≤ 30 | 范围不同 |
| Thermal level | GFN2-xTB//xTB | 默认 GFN1，失败可能 GFN0 | 方法层级不同 |
| 热校正前筛选 | ensemble free-energy 语义 | B97 electronic window 后再 mRRHO | 排序顺序不同 |
| 代表构象优化 | B97-3c OPT | S1 不做，S3 做 B97-3c OPT/OptTS | 阶段与对象不同 |
| 代表构象最终 SP | ωB97M-V on B97-3c geom | ωB97M-V on M062X geom (S4) | geometry level 不同 |
| 最终组装 | $G_0^{\mathrm{RSH//GGA}} + G_{\mathrm{conf}}^{\mathrm{rel}}$ | 未加 conformational correction | 核心语义缺失 |
| mRRHO 失败 | 无此降级（全 thermal） | 任一失败 → pool 退回电子能 | RPH 特有降级 |
| 适用对象 | 55 reaction states benchmark | precursor/product/intermediate/TS | 外推范围不同 |

### 11.6 提前漏斗与 ensemble 语义不一致

RPH 的两轮固定数量截断 + "先 electronic 筛选后 thermal correction"的顺序改变了进入 partition function 的构象集合。即使保留构象的单点 score 与论文相近，也不等于 ensemble free-energy 语义。

### 11.7 Representative refinement 与 S3/S4 不能等同

论文代表构象精修链：low-level rank-1 → B97-3c OPT → ωB97M-V SP。
RPH：S1 rank-1 → S2 PEB（引入 TS/intermediate）→ S3 B97-3c OPT + r2SCAN-3c SP → S4 M062X OPT + ωB97M-V SP。S3/S4 计算对象不再只是 S1 rank-1 平衡态构象，S4 ωB97M-V 在 M062X 而非 B97-3c 优化几何上计算。

### 11.8 TS 管线适用边界

RTCONF55-16K 验证的是反应物/产物等反应状态，TS 构象评估仍在进行。论文的 MAE ≈ 0.6 (CENSO-light) / 0.9 (CENSO-zero) kcal/mol 及 "~10x faster" 不能直接作为 RPH 过渡态和反应剖面的误差或加速比。

### 11.9 定性结论

RPH 已实现论文 CENSO-light 的部分计算思想——CREST/GFN2 ensemble、B97-3c ranking、xTB thermal correction、rank-1 进入后续管线、不做全 ensemble DFT optimization——但尚未实现决定其最终物理语义的核心部分：partition function、$G_{\mathrm{conf}}^{\mathrm{rel}}$、Boltzmann population、rank-1 refined energy 与 conformational correction 的组装、GFN2 thermal level、RSH//GGA representative refinement。

当前更准确的定义为：**CENSO-light-inspired rank-1 conformer selection stage**，非完整的论文 `CENSO-light ensemble free-energy protocol`。

## 12. 2026-07-14 ensemble thermodynamics 重构状态

本节记录第 11 节审计之后已经落地的代码变化。第 11 节保留为“重构前差异基线”，不应再被解释为当前实现状态。

### 12.1 已解决的核心语义

S1 manifest 已升级为 `s1_censo_lite_v2`，当前对通过 CREST 和几何去重的完整唯一构象集合计算：

$$
G_i = E_i^{\mathrm{B97-3c}} + G_i^{\mathrm{RRHO,xTB}}
$$

$$
Z_{\mathrm{rel}} = \sum_i d_i\exp\left(-\frac{\Delta G_i}{RT}\right)
$$

$$
p_i = \frac{d_i\exp(-\Delta G_i/RT)}{Z_{\mathrm{rel}}}, \qquad
G_{\mathrm{conf}}^{\mathrm{rel}}=-RT\ln Z_{\mathrm{rel}}
$$

具体变化：

- 默认不再使用额外的 xTB 5 kcal/mol window、60 构象 cap、B97 electronic window 和 30 构象 thermo cap；CREST 自身 6 kcal/mol window 内的唯一构象保留到 partition function。
- B97-3c 与 mRRHO 使用同一 ensemble；不再允许在热校正之前按 B97 电子能删除成员。
- 每个候选输出 `relative_free_energy_kcal`、`relative_partition_weight`、`boltzmann_population` 和 `degeneracy`。
- ensemble 汇总输出 `partition_function_relative`、`conformational_free_energy_correction_kcal`、`population_sum`、reference free energy 和 ensemble free-energy estimate。
- `candidates` 保存完整 ensemble；`representative_candidates` 只表示便于后续选种和查看的低能子集，不再决定 partition function 成员。
- 默认 B97 和 mRRHO failure policy 均为 `strict`，避免以缺失成员或混合能量层级静默构造 partition function。

当前 `degeneracy` 默认值为 1。CREST 重复采样次数不会被误当作物理简并度；只有后续存在可靠的对称/简并度来源时才应改变该值。

### 12.2 与 S3/S4 最终能量的连接

对于直接来源于 S1 的 precursor/product minimum，S1 同时输出：

```text
ensemble_thermochemistry_correction_hartree
  = G_RRHO(reference)_S1 + G_conf_rel_S1
```

该 correction 随结构记录传入 S3 和 S4。若相应阶段单点能成功，则 manifest 增加：

```text
ensemble_corrected_sp_free_energy_hartree
  = E_stage_SP + G_RRHO(reference)_S1 + G_conf_rel_S1
```

因此 ensemble 不再只用于选 rank-1；其热统计贡献可追踪地进入后续精修能量。S2 生成的 intermediate/TS 没有对应的平衡态构象 ensemble，目前不会借用 product correction。

### 12.3 并行资源模型

B97-3c 和 mRRHO 均按 `resources.nproc` 分配总核心预算：

```text
parallel_jobs = min(floor(total_cores / cores_per_job), max_parallel_jobs)
```

默认 16 核配置下：

| 环节 | 并行任务 | 每任务核心 | 最大占用核心 |
|---|---:|---:|---:|
| ORCA B97-3c SP（主体） | 16 | 1 | 16 |
| ORCA B97-3c SP（自适应尾部） | ≤8 | ≥2 | ≤16 |
| xTB mRRHO | 8 | 1 | 8 |

mRRHO 不通过 CREST 执行；这里采用的是与 CREST 类似的受控 fan-out 思路。单核 mRRHO 同时规避了已定位的 xTB 6.7.1 parallel BHESS 线程栈崩溃，默认 `OMP_STACKSIZE=1G`，崩溃重试保持同一 GFN2 层级并降低到 1 核，不再默认切换 GFN0。

### 12.4 仍与论文 CENSO-light 不同的边界

- RPH S1 仍遵守 V4 theory contract，不执行 B97-3c geometry optimization；论文 representative refinement 的阶段定义不同。
- CREST 搜索模式仍为 `imtd_gc`，并非论文写明的 iMTD-sMTD。
- RPH 使用 acetone，论文 benchmark 使用 DCM；这是目标体系条件差异，不应为复刻 benchmark 而盲目改动。
- S4 ωB97M-V 单点位于 M062X 几何，而非论文 RSH//GGA 的 B97-3c 几何。
- 当前只为 S1 直接产生的平衡态 minimum 组装 conformational correction；TS/intermediate 需要各自独立的 ensemble 定义与验证。
- 本轮只完成静态检查和 mock/unit tests，尚未替代 `docs/WSL_TEST_PLAN_V4.md` 的 T0–T6 实机门禁，也没有据此宣称真实加速倍数。

因此，重构后的准确表述是：**RPH 已实现 CENSO-light-style ensemble partition thermodynamics，并将 correction 接入 S3/S4；但整体几何精修链和适用对象仍不是论文协议的逐项复刻。**

## 13. 首次 2×8 实机运行暴露的 mRRHO 坐标单位错误

首次 `rx-id=1 --stop-after s1` 实测确认 ORCA B97-3c 的 `2 jobs × 8 cores` fan-out 正常工作：precursor 191/191 个 SP 全部成功，耗时约 35.5 分钟。

运行报告最初将后续失败描述为“GFN1 SCC 不收敛”，但原始命令日志表明实际使用的是：

```text
--gfn 2 --bhess normal --enso --parallel 1 --alpb acetone
```

失败构象为 `conf_0155`、`0158`、`0160`、`0164`、`0183`、`0190` 和 `0191`，均在默认 250 次 iteration 后生成 `.sccnotconverged`。

进一步定位发现根因不在 `failure_policy=strict`，而在 `_xyz_to_coord()`：XYZ 坐标单位为 Å，RPH 却将数值原样写入普通 Turbomole `$coord`；普通 `$coord` 默认单位为 Bohr，因此 xTB 实际读取到缩短为 0.529 倍的畸变几何。该错误会影响全部 mRRHO，而不仅是 7 个显式失败构象。

修复内容：

- 写入 `$coord` 前使用 `BOHR_TO_ANGSTROM` 执行 Å→Bohr 转换。
- 保持 GFN2、ALPB、SCC accuracy 和 electronic temperature 不变。
- 保留 `strict`，避免以缺失构象构造 partition function。
- 增加同一 GFN 层级的 SCC iteration retry（默认 1000）作为真正的收敛后备，而非切换 GFN 方法。

隔离复测结果：上述 7 个原失败构象在修正坐标单位后，全部在默认 250 次 SCC 上限内成功生成 `xtb_enso.json`，无需放宽收敛精度或提高 electronic temperature。

由于此前成功的 product mRRHO 同样使用了错误单位，首次运行得到的 `G_conf=-0.66/-0.68 kcal/mol` 不应继续使用；必须以修复后的坐标重新计算。新增配置项会改变 S1 checkpoint signature，从而阻止错误的 product S1 manifest 被静默复用。

## 14. B97-3c SP 自适应尾部资源调度

固定 `16 jobs × 1 core` 在主体阶段具有更高 ensemble throughput，但最后少量任务会导致空闲核心。当前调度器将最后8个 SP 标记为 adaptive tail：

- bulk 阶段持续维持最多16个单核 ORCA 任务。
- 最后8个任务不会提前以单核全部启动。
- bulk 任务每释放至少2个核心，调度器立即启动一个更宽的 tail SP。
- 所有运行中任务的 `nprocs` 总和始终不超过 `resources.nproc`。
- 不要求整个 bulk 批次先完成，因此避免简单 wave/batch 调度的完整栅栏等待。

默认16核示意：

```text
bulk: 16 × 1
tail=8: 8 × 2
tail=4（小 ensemble）: 4 × 4
tail=2（小 ensemble）: 2 × 8
tail=1（小 ensemble）: 1 × 16
```

当 bulk 与 tail 暂时重叠时，实际布局可能为 `12×1 + 2×2`、`8×1 + 4×2` 等，但 active core sum 始终不超过16。该策略优化的是尾部资源利用率；实际 wall-time 收益仍取决于 ORCA 对单任务多核的缩放效率。
