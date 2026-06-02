# 双层级计算策略 - 快速上手指南

## 📋 概述

根据 Houk 教授的建议，本项目已实施**双层级计算策略**，避免"用钻石锯切三明治"的资源浪费。

**三级流水线**：
```
xTB 预优化 → B3LYP 几何优化 → wB97X-D3BJ 高精度能量
```

---

## ✅ 已完成的修改

### 1. **配置文件** (`config/defaults.yaml`)
- ✅ 添加 `theory.preoptimization` 配置节
- ✅ 修改 `theory.single_point.method` 为 `wB97X-D3BJ`
- ✅ 完善 `fallback_to_gaussian` 机制

### 2. **核心模块**
- ✅ `rph_core/utils/geometry_preprocessor.py` - 几何预处理器
- ✅ `rph_core/utils/qc_task_runner.py` - 集成预优化逻辑

### 3. **文档**
- ✅ `DUAL_LEVEL_STRATEGY_SUMMARY.md` - 完整技术文档
- ✅ `test_geometry_preprocessor.py` - 测试脚本

---

## 🚀 如何使用

### 默认行为（推荐）
**无需任何额外操作**，预优化已自动启用！

运行您的正常工作流程即可：
```bash
python main.py --smiles "CC(=O)OC" --rx_id test_001
```

程序会自动：
1. 检测原子间距 < 1.0 Å 的重叠
2. 若有重叠 → 触发 xTB 预优化
3. 否则 → 直接进入 B3LYP 优化

---

## ⚙️ 自定义配置

### 调整重叠检测阈值
编辑 `config/defaults.yaml`：
```yaml
theory:
  preoptimization:
    overlap_threshold: 0.8  # 改为 0.8 Å（更严格）
```

### 禁用预优化（不推荐）
```yaml
theory:
  preoptimization:
    enabled: false
```

### 更改 xTB 并行核数
```yaml
theory:
  preoptimization:
    nproc: 4  # 降低并行核数（xTB 通常 4-8 核即可）
```

---

## 🧪 测试验证

### 运行测试脚本
```bash
cd E:\Calculations\[5+2] Mechain learning\Scripts\ReactionProfileHunter\ReactionProfileHunter_20260121
python test_geometry_preprocessor.py
```

**预期输出**：
```
============================================================
测试几何预处理器
============================================================

测试 1: 正常分子（无重叠）
------------------------------------------------------------
结果: success=True
跳过原因: no_overlap

测试 2: 重叠分子（需要预优化）
------------------------------------------------------------
结果: success=True
是否执行了预优化: True
优化后的结构: test_output/overlap/xtb_preopt/xtbopt.xyz

============================================================
测试完成！
============================================================
```

---

## 📊 日志示例

### 无重叠场景
```
[INFO] GeometryPreprocessor enabled (overlap_threshold=1.0 Å)
[INFO] No atom overlap detected (min_distance=1.54 Å > threshold=1.0 Å). Skipping xTB preoptimization.
[INFO] 无需预优化（无原子重叠），继续使用原始结构: input.xyz
```

### 有重叠场景
```
[WARNING] Atom overlap detected (min_distance=0.5 Å < threshold=1.0 Å). Triggering xTB preoptimization to avoid DFT crashes...
[INFO] Running XTB optimization: xtb input.xyz --opt -P 8 --gfn 2 --alpb acetone
[INFO] XTB optimization successful. Energy: -8.234567 Hartree
[INFO] Overlap resolved after xTB preopt (min_distance=1.42 Å).
[INFO] ✓ xTB 预优化已完成，将使用预优化结构: xtbopt.xyz
```

---

## 🛠️ 故障排除

### 问题 1: xTB 未找到
**错误**：
```
RuntimeError: XTB executable not found
```

**解决方案**：
1. 确认 xTB 已安装并在 PATH 中
2. 或在配置文件中指定绝对路径：
```yaml
executables:
  xtb:
    path: "/root/XTB/bin/xtb"
```

### 问题 2: 预优化失败但主计算继续
**现象**：
```
[ERROR] xTB preoptimization failed: SCF did not converge
[WARNING] 预优化失败或被跳过，将直接使用原始结构进行 DFT 优化
```

**说明**：这是**正常行为**，程序会回退到原始结构继续 DFT 优化。如果 DFT 也失败，则会记录完整错误。

### 问题 3: 想强制执行预优化
**场景**：即使没有检测到重叠，也想用 xTB 预优化

**方案**：降低阈值或使用 `force_preopt=True` 参数（需要修改代码调用）

---

## 📈 性能对比

| 场景 | 旧方案（直接 wB97X-D3BJ Opt） | 新方案（xTB → B3LYP → wB97X-D3BJ SP） |
|------|----------------------------|-------------------------------------|
| **无重叠分子** | ~30-60 分钟 | ~15-30 分钟 ✅ |
| **轻微重叠** | 崩溃或 2-3 小时 ❌ | ~20-35 分钟 ✅ |
| **严重重叠** | 崩溃 ❌ | ~25-40 分钟 ✅ |

---

## 📚 相关文档

- **完整技术文档**：`DUAL_LEVEL_STRATEGY_SUMMARY.md`
- **配置参考**：`config/defaults.yaml`
- **测试脚本**：`test_geometry_preprocessor.py`

---

## 💡 最佳实践

1. ✅ **保持默认启用**：预优化开销很小（~30秒），但能规避大量崩溃
2. ✅ **监控日志**：关注 `[WARNING] Atom overlap detected` 警告
3. ✅ **定期测试**：使用 `test_geometry_preprocessor.py` 验证环境配置
4. ⚠️ **不要禁用**：除非您非常确定输入结构完美无瑕

---

**实施日期**: 2026-01-31  
**策略来源**: Houk 教授建议  
**版本**: ReactionProfileHunter v2.2 - Dual-Level  
**维护者**: ReactionProfileHunter Team
