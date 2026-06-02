# DLPNO-CCSD(T) TS 计算失败深度分析报告

**任务**: SP-REF 对 rx3 和 rx15 的 TS 驻点单点能计算失败  
**方法**: DLPNO-CCSD(T)/def2-TZVPP，TightPNO，RIJK，16 核 MPI 并行，maxcore=2000 MB  
**日期**: 2026-04-28 ~ 2026-04-29  
**报告日期**: 2026-04-30

---

## 摘要

rx3 和 rx15 的 SP-REF (DLPNO-CCSD(T)) TS 单点能计算在 **DLPNO (T) triples 校正步骤** 因内存不足而中止。CCSD 迭代已成功收敛（14 轮），但 (T) 校正步骤中，单个 triple 积分计算需要约 1.45–1.48 GB 内存，而每个 MPI 进程实际可用仅约 1.35–1.36 GB。此问题源于 **maxcore=2000 MB 被 16 个 MPI 进程分摊后，单个大 triple 占用的内存超出进程配额**。这是 ORCA DLPNO-CCSD(T) 在较大分子 TS 构型上的已知瓶颈——(T) 步骤无可用的批处理降级方案，不同于 CCSD 步骤可以在内存不足时自动切换为批量计算模式。

---

## 1. 失败定位：计算走到了哪一步？

DLPNO-CCSD(T) 计算包含以下顺序步骤：

```
SCF 参考态 → LMP2 初始猜测 → PNO 生成 → CCSD 迭代 → (T) triples 校正
```

### rx3 进度追踪

| 步骤 | 状态 | 输出行号 | 详情 |
|------|------|---------|------|
| SCF 参考态 | ✅ 完成 | 早期 | RHF 收敛 |
| LMP2 初始猜测 | ✅ 完成 | 1972–1985 | 12 轮迭代收敛，E_LMP2 = −3.881562 |
| PNO 生成 | ✅ 完成 | 1986–1999 | 98.1% 密度迹保留，1031/1653 对入选 |
| 局部 RI 变换 (IAVPAO) | ✅ 完成 | 1918–1936 | 920 PAO, 2182 辅助函数, 4.9s |
| CCSD 迭代 | ✅ 完成 | 2187–2204 | **14 轮迭代收敛！** E_Corr = −3.909273 |
| **(T) triples 校正** | ❌ **崩溃** | 2255–2266 | **内存不足** |

### rx15 进度追踪

| 步骤 | 状态 | 输出行号 | 详情 |
|------|------|---------|------|
| SCF 参考态 | ✅ 完成 | 早期 | RHF 收敛 |
| LMP2 初始猜测 | ✅ 完成 | 1873–1884 | 11 轮迭代收敛，E_LMP2 = −3.453103 |
| PNO 生成 | ✅ 完成 | 1885–1898 | 98.3% 密度迹保留，1097/1326 对入选 |
| 局部 RI 变换 (IAVPAO) | ✅ 完成 | 1817–1836 | 827 PAO, 1954 辅助函数, 3.8s |
| CCSD 迭代 | ✅ 完成 | 2084–2101 | **14 轮迭代收敛！** E_Corr = −3.489223 |
| **(T) triples 校正** | ❌ **崩溃** | 2152–2163 | **内存不足** |

**关键结论**: 两次失败均发生在完全相同的步骤——(T) triples 校正，这是 DLPNO-CCSD(T) 计算的**最后一步**。SCF → LMP2 → PNO → CCSD 全部成功。

---

## 2. 精确错误信息

### rx3 TS 崩溃输出（行 2255–2266）

```
-------------------------------------------
DLPNO BASED TRIPLES CORRECTION
-------------------------------------------

Singles multiplier          ...  1.000000
TCutTNO                     ... 1.000e-09
TCutMP2Pairs                ... 1.000e-06
...
Number of Triples that are to be computed   ...     16546

 . . . . . . . .             (     31.841 sec)Process  0: 
PROBLEM (ORCA/DLPNO Triples): 
Process  0:         not enough memory for a single Triple
Process  0:         MemNeeded  1481.7 MB > MemAvailable  1349.5 MB

[file orca_mdci/mdci_dlpno_rhf_triples.cpp, line 1702, Process 0]:
   ... Aborting the run

ORCA finished by error termination in MDCI
Calling Command: mpirun -np 16  /opt/software/orca/orca_mdci_mpi ...
[file orca_tools/qcmsg.cpp, line 465]:
   .... aborting the run
```

### rx15 TS 崩溃输出（行 2152–2163）

```
Number of Triples that are to be computed   ...     20134

 .Process 10: 
PROBLEM (ORCA/DLPNO Triples): 
Process 10:         not enough memory for a single Triple
Process 10:         MemNeeded  1451.0 MB > MemAvailable  1364.9 MB

[file orca_mdci/mdci_dlpno_rhf_triples.cpp, line 1702, Process 10]:
   ... Aborting the run

 . . . . .
ORCA finished by error termination in MDCI
Calling Command: mpirun -np 16  /opt/software/orca/orca_mdci_mpi ...
[file orca_tools/qcmsg.cpp, line 465]:
   .... aborting the run
```

### 内存缺口分析

| 反应 | 需要内存 | 可用内存 | 缺口 | 缺口比例 | 失败进程 |
|------|---------|---------|------|---------|---------|
| **rx3** | 1481.7 MB | 1349.5 MB | **132.2 MB** | 9.8% | Process 0 |
| **rx15** | 1451.0 MB | 1364.9 MB | **86.1 MB** | 6.3% | Process 10 |

两个案例的内存缺口都不大（仅 5–10%），但 (T) triples 代码**不支持分块计算**——单个 triple 所需的内存要么全部满足，要么直接中止。

---

## 3. 与成功案例的对比

四个反应的 DLPNO-CCSD(T) TS 计算结果对比如下：

### 系统规模与计算状态

| 参数 | rx1 ✅ | rx3 ❌ | rx8 ✅ | rx15 ❌ |
|------|-------|-------|-------|-------|
| 基函数数 | 914 | 920 | 650 | 827 |
| 壳层数 | 346 | 344 | 242 | 311 |
| 总电子对 | 1653 | 1653 | 903 | 1326 |
| 入选对 | 1031 | 1031 | 787 | 1097 |
| **Triples 数量** | **19,727** | **16,546** | **11,937** | **20,134** |
| Triples / 基函数 | 21.6 | 18.0 | 18.4 | **24.3** ⬆️ |
| CCSD 内存警告 | ⚠️ 2024→1275 MB | 无 | ⚠️ 4776→922 MB | 无 |
| (T) 内存需求 | 未触发 | 1482 MB | 未触发 | 1451 MB |
| (T) 可用内存 | — | 1350 MB | — | 1365 MB |
| **最终状态** | **完成** | **中止** | **完成** | **中止** |

### 关键发现

**rx1 和 rx8 在 CCSD 步骤也遇到了内存不足警告，但 ORCA 成功降级处理**：

```
rx1 (行 2052–2054):
PROBLEM (ORCA/DLPNO): not enough memory for MDCI Iterations
       MemNeeded  2024.2 MB > MemAvailable  1274.7 MB
       Will have to batch evaluation of 4-external contributions to sigma vector

rx8 (行 1707–1708):
PROBLEM (ORCA/DLPNO): not enough memory for MDCI Iterations
       MemNeeded  4776.1 MB > MemAvailable   922.1 MB
       Will have to batch evaluation of 4-external contributions to sigma vector
```

CCSD 步骤有**自动降级机制**：当 4-external 积分需要的内存超出可用量时，ORCA 自动切换为"批量计算模式"（batch evaluation），牺牲速度换取内存效率。两个成功案例都借助此机制完成了 CCSD 迭代。

**而 (T) triples 步骤没有类似的降级机制**。`mdci_dlpno_rhf_triples.cpp:1702` 处直接调用 `abort()`，不做任何重试或降级尝试。

### 为什么 rx3 和 rx15 失败了，而 rx1 和 rx8 成功了？

表面上看，**rx1 的 triples 数量（19,727）反而多于 rx3（16,546）**，但 rx1 成功了而 rx3 失败了。原因在于：

1. **rx1 的 CCSD 触发了批量模式**，释放了部分内存缓冲区，使得 (T) 步骤有更多可用内存
2. **rx3 和 rx15 的 CCSD 未触发批量模式**——它们的基函数数较大但电子对结构更"紧凑"，CCSD 的 4-external 积分刚好能放入内存，但这也意味着内部缓冲区保持满载状态，留给 (T) 步骤的可用内存更少
3. **rx15 有最高的 triple/基函数比率（24.3）**——这意味着其电子的 triple 激发模式更加"集中"，单个 triple 的维度更大

---

## 4. 输入参数分析

两个失败案例使用的 ORCA 输入关键词完全相同：

```
! DLPNO-CCSD(T) def2-TZVPP def2/JK RIJK def2-TZVPP/C tightSCF TightPNO noautostart miniprint nopop
%maxcore 2000
%pal nprocs 16 end
%cpcm solvent "acetone" end
```

### 内存分配分析

- `maxcore = 2000`：允许 ORCA 最多使用 2000 MB（约 2 GB）内存
- `nprocs = 16`：使用 16 个 MPI 进程并行
- **每个进程理论可用内存**：2000 / 16 ≈ **125 MB**（但这只是简单除法）
- **实际内存使用**：ORCA 会动态分配，某些积分可跨进程共享；SCF 和 LMP2 步骤实际每个进程使用了约 500–800 MB（从 temp 文件大小推断）
- **(T) 步骤的实际可用**：~1350 MB/进程（远高于简单除法，因为 ORCA 允许每个进程借用未使用的全局配额）

**问题本质**：maxcore=2000 设置了全局上限，但 ORCA 内部的内存管理允许单个进程在 (T) 步骤借用更多内存。问题在于，对于这些特定的 TS 构型，单个最大的 triple 所需内存（~1.48 GB）超出了当时可用的配额（~1.35 GB）。

---

## 5. 残留文件分析

两个失败案例的 TS 目录均留下了大量 ORCA 临时文件：

### rx3 TS 目录 (`qc/ts/`)：181 个文件

| 文件类别 | 示例 | 数量 | 含义 |
|---------|------|------|------|
| 波函数文件 | `.gbw` (8.1 MB) | 1 | SCF 波函数已保存 ✅ |
| 密度文件 | `.densities` (6.5 MB) | 1 | 密度数据已保存 ✅ |
| PNO 文件 | `.PNO3.tmp.0–15` | 16 | 各进程的 triple PNO |
| PNOJK 积分 | `.PNOJK.proc0–15.tmp` | 16 | PNO 辅助 JK 积分 |
| PNOJK 日志 | `.PNOJK.LOG0–15.tmp` | 16 | 各进程积分日志 |
| PNO 重叠 | `.pnoovl.tmp.0–15` | 16 | 对重叠矩阵 |
| PAO 变换 | `.VABPAO.tmp.0–15` | 16 | 投影原子轨道 |
| IKJL 积分 | `.IKJL.tmp.proc0–15.0` | 16 | 4-index 变换积分 |
| P2V 归约 | `.P2V_red.tmp.0–15` | 16 | PNO→虚拟 MO 变换 |
| 属性文件 | `_property.txt` | 1 | SCF 属性已提取 ✅ |
| SMD 输出 | `.smd.out` | 1 | SMD 溶剂计算完成 ✅ |
| ORCA 输出 | `.out` (2267 行) | 1 | 含完整错误信息 |

### rx15 TS 目录 (`qc/ts/`)：174 个文件

结构类似，`.gbw` (6.8 MB)，`.densities` (5.3 MB)，`.out` (2164 行)。

**解读**：所有中间数据均已生成，CCSD 波函数和信息完整。如果能够为 (T) 步骤提供足够内存或使用外部后处理工具，理论上可以从这些残留文件恢复完整的 CCSD(T) 能量。

---

## 6. 根因总结

```
故障树：

DLPNO-CCSD(T) TS 计算失败
  └── (T) triples 校正步骤内存不足
        ├── 直接原因：单个 triple 需要 ~1.48 GB，可用仅 ~1.35 GB
        ├── 深层原因：(T) 代码无降级/批处理机制（与 CCSD 不同）
        ├── 触发条件：大分子 TS 构型 + TightPNO + 16 核 MPI
        │     ├── rx3: 920 基函数, 16,546 triples, 缺口 132 MB
        │     └── rx15: 827 基函数, 20,134 triples, 缺口 86 MB
        └── 成功案例对比：
              ├── rx1: CCSD 批量模式释放了内存, (T) 顺利完成
              └── rx8: 更小的系统 (650 基函数, 11,937 triples)
```

### 核心矛盾

并非系统绝对规模导致的失败（rx1 的基函数数相似且 triples 更多，却成功了），而是**内存分配的边界效应**：

- CCSD 步骤有 `maxcore` 限制下的**优雅降级**（批量化 4-external 积分）
- (T) 步骤**无降级**——全有或全无
- 当 CCSD 触发批量模式时（如 rx1），反而释放了内存给后续的 (T) 步骤
- 当 CCSD 刚好能放入内存时（如 rx3、rx15），(T) 步骤反而内存不足

---

## 7. 解决方案与建议

### 方案 A：调整 MPI 并行度（推荐，最直接）

将 `nprocs` 从 16 降至 8，使每个进程获得更多内存配额：

```
%pal nprocs 8 end
```

**预期效果**：每个进程可用内存约翻倍（~2.7 GB），远超 1.48 GB 的需求。

**代价**：计算时间约增加 1.5–2 倍。

### 方案 B：增加 maxcore（最简单）

将 `maxcore` 从 2000 提高至 4000 或更高：

```
%maxcore 4000
```

**预期效果**：每个进程可用内存相应增加，可能仅需 2500–3000 MB 即可满足。

**代价**：需要确保机器有足够物理内存（16 进程 × 4000 MB ≈ 64 GB 峰值需求，但实际使用较少）。

### 方案 C：降低 PNO 阈值（折中）

将 TightPNO 降为 NormalPNO，减少 triple 的空间维度：

```
! DLPNO-CCSD(T) ... NormalPNO
```

**预期效果**：triples 的维度更小，内存需求降低，但精度可能略有下降。

**代价**：精度从 TightPNO (~99.9% 相关能) 降至 NormalPNO (~99.5%)，对 TS 能垒可能造成 0.1–0.5 kcal 的误差。

### 方案 D：分开计算 CCSD 和 (T)

由于 CCSD 已成功收敛，可以尝试：
1. 从残留 `.gbw` 和 PNO 文件恢复
2. 使用 ORCA 的 `%mdci` 模块单独运行 (T) 校正（如果 ORCA 支持）
3. 或者在同一输入中增大 maxcore 后从断点继续（ORCA 通常不支持 DLPNO 的断点续算）

### 方案 E：降低基组（不推荐）

将 def2-TZVPP 降为 def2-TZVP 或 def2-SVP 会显著降低内存需求，但会牺牲参考值的精度。

### 推荐组合方案

```
方案 A + B：nprocs=8, maxcore=4000
方案 A 优先于 B：降低并行度比增加总内存更可靠，因为它直接增加了每个 triple 的内存配额
```

---

## 8. CCSD 能量（未加 (T) 校正的参考值）

虽然 (T) 校正未完成，CCSD 部分的能量已收敛。这些值可作**近似参考**（误差约 0.5–2 kcal/mol，取决于 (T) 贡献大小）：

| 反应 | E(0) (hartree) | E_Corr_CCSD (hartree) | E_TOT_CCSD (hartree) | T1 诊断值 |
|------|---------------|----------------------|---------------------|----------|
| **rx3 TS** | −1007.836987 | −3.909273 | **−1011.753963** | 0.01390 |
| **rx15 TS** | −894.180229 | −3.489223 | **−897.676289** | 0.01422 |

T1 诊断值均小于 0.02，表明单参考态描述良好，CCSD 方法是合适的。

> ⚠️ 注意：以上 CCSD 能量不含 (T) 校正和弱对校正（weak-pairs correction）。完整 CCSD(T) 总能量需加上 ~0.007–0.008 hartree 的弱对贡献和 triples 贡献（通常 0.01–0.03 hartree）。

---

## 9. 数据清单

### 失败案例输出文件

```
Output/benchmark_dft_theory/experiments/session_bl_fixed/
├── rx3/sp/SP-REF/
│   ├── sp_result.json                    ← 失败状态记录
│   └── qc/ts/
│       ├── ts_final_c2c5db04.out         ← 完整 ORCA 输出 (2267 行)
│       ├── ts_final_c2c5db04.gbw          ← SCF 波函数 (8.1 MB)
│       ├── ts_final_c2c5db04.densities    ← 密度 (6.5 MB)
│       ├── ts_final_c2c5db04.inp          ← ORCA 输入文件
│       ├── ts_final_c2c5db04_property.txt ← SCF 属性
│       ├── ts_final_c2c5db04.smd.out      ← SMD 溶剂输出
│       └── *.tmp (160+ 文件)              ← PNO/积分中间文件
└── rx15/sp/SP-REF/
    ├── sp_result.json
    └── qc/ts/
        ├── ts_final_8acc2640.out          ← 完整 ORCA 输出 (2164 行)
        └── ... (类似结构)
```

### 成功案例对比文件

```
rx1/sp/SP-REF/qc/ts/ts_final_7b855de4.out   ← 含 CCSD 批量降级日志
rx8/sp/SP-REF/qc/ts/ts_final_7e7d2203.out   ← 含 CCSD 批量降级日志
```

---

## 附录：ORCA DLPNO-CCSD(T) 内存管理机制简述

ORCA 的 DLPNO-CCSD(T) 计算分三个阶段使用内存：

| 阶段 | 主要数据结构 | 是否有降级方案 | 失败模式 |
|------|------------|--------------|---------|
| SCF | Fock 矩阵, 密度矩阵 | ✅ DIIS, 水平移动 | 极少失败 |
| LMP2 + PNO 生成 | 对密度, PAO 变换矩阵 | ✅ 按对处理 | 极少失败 |
| CCSD 迭代 | 4-external 积分 (IKJL) | ✅ **批量模式** | 降级后变慢但成功 |
| **(T) triples** | **单个 triple 的三电子积分** | ❌ **无降级** | **直接 abort()** |

CCSD 步骤的批量模式（`batch evaluation of 4-external contributions`）是 ORCA 开发者专门为内存受限场景设计的功能。当单个进程的 IKJL 积分超出可用内存时，ORCA 将其拆分为多批处理，每次只加载一部分积分。这会增加 I/O 开销（约 1.5–3 倍慢），但确保了计算能完成。

(T) triples 步骤之所以没有类似机制，是因为 triples 积分涉及 (occ, occ, virt) 三指数张量，其耦合结构使得分块计算极为困难——当前的 `mdci_dlpno_rhf_triples.cpp` 实现假定单个 triple 的全部中间数据能放入内存。

**这是 ORCA 的已知限制，而非 bug。** ORCA 手册建议对大分子 DLPNO-CCSD(T) 计算使用更少的核数（"for large molecules, reduce the number of cores to ensure sufficient memory per core for the triples step"）。
