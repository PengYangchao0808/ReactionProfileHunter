# ReactionProfileHunter 修改日志 / Modification Log

## v2.0.0 版本更新 (Current / 当前版本)

### 发布日期 / Release Date: 2025-03-18

---

## 新增功能 / New Features

### 1. 前向扫描过渡态搜索 (Forward-Scan TS Search) **[新架构]**
- **描述 / Description**: 实现基于 xTB 原生 `$scan` 的任意环加成反应过渡态搜索
- **支持的反应类型 / Supported Reaction Types**:
  - `[4+3]` - 4+3 环加成 (前向扫描)
  - `[4+2]` - Diels-Alder 反应 (前向扫描)
  - `[3+2]` - 1,3-偶极环加成 (前向扫描)
  - `[5+2]` - 5+2 环加成 (逆向扫描)
- **相关文件 / Related Files**:
  - `rph_core/steps/step2_retro/retro_scanner.py`
  - `rph_core/utils/xtb_runner.py`

### 2. 反应配置文件 (Reaction Profiles)
- **描述 / Description**: 按反应类型配置扫描参数
- **配置示例 / Config Example**:
  ```yaml
  reaction_profiles:
    "[4+3]_default":
      forming_bond_count: 2
      s2_strategy: forward_scan
      scan:
        scan_start_distance: 1.8
        scan_end_distance: 3.2
        scan_steps: 20
  ```
- **相关文件 / Related Files**: `config/defaults.yaml`

### 3. S1 活化特征 (Step1 Activation Features)
- Boltzmann/Gibbs 加权能量
- 构象熵与柔性
- 离去基团几何
- α-H 门控机制
- **相关文件 / Related Files**:
  - `rph_core/steps/step4_features/extractors/step1_activation.py`

### 4. S2 环化特征 (Step2 Cyclization Features)
- 动力学/热化学
- 过渡态几何 (成键)
- CDFT 指标 (Fukui/Dual Descriptor/QTAIM)
- **相关文件 / Related Files**:
  - `rph_core/steps/step4_features/extractors/step2_cyclization.py`

### 5. Multiwfn 集成 (Multiwfn Integration)
- Tier-1 支持 (Fukui/Dual Descriptor/QTAIM)
- 非交互式批处理
- 故障容错
- 缓存机制
- **相关文件 / Related Files**:
  - `rph_core/utils/multiwfn_runner.py`
  - `rph_core/steps/step4_features/extractors/multiwfn_features.py`

### 6. 小分子缓存系统 (Small Molecule Cache)
- 全局小分子复用
- 缓存键生成
- 检测与匹配
- **相关文件 / Related Files**:
  - `rph_core/utils/small_molecule_cache.py`
  - `rph_core/utils/small_molecule_catalog.py`

---

## 重构与优化 / Refactoring & Optimizations

### 1. 目录结构重构 (Directory Structure Refactor)
- **变更 / Change**: `S1_Product` 重命名为 `S1_ConfGeneration`
- **原因 / Reason**: 更清晰的语义表达
- **相关提交 / Related Commits**:
  - `b37a2d1` - refactor(s1): rename S1_Product to S1_ConfGeneration

### 2. S2 路径更新 (S2 Path Updates)
- 更新 S2 相关路径引用
- 确保向后兼容性
- **相关提交 / Related Commits**:
  - `b37a2d1` - refactor(s1): rename S1_Product to S1_ConfGeneration and update S2 paths

### 3. 项目结构清理 (Project Structure Cleanup)
- 清理过期备份目录 (`rph_core_backup_20260115/`)
- 更新 `.gitignore`
- 优化测试文件组织
- **相关提交 / Related Commits**:
  - `ca5e2f9` - chore: cleanup project structure and update .gitignore

---

## 测试增强 / Testing Enhancements

### 新增测试文件 / New Test Files
| 文件 / File | 描述 / Description |
|------------|-------------------|
| `tests/test_s1_progress_parser.py` | S1 进度解析器测试 |
| `tests/test_forward_scan_wiring.py` | 前向扫描接线测试 |
| `tests/test_geometry_guard.py` | 几何守护测试 |
| `tests/test_s2_boundary_degrade.py` | S2 边界降级测试 |
| `tests/test_s3_checkpoint.py` | S3 检查点测试 |
| `tests/test_s4_extractor_degrade_behavior.py` | S4 提取器降级行为测试 |
| `tests/test_s4_meta_warnings_and_weights.py` | S4 元数据警告与权重测试 |
| `tests/test_s4_v62_final_verification.py` | S4 V6.2 最终验证测试 |
| `tests/test_qc_interface_gaussian_failures.py` | Gaussian 失败处理测试 |
| `tests/test_xtb_path_integration.py` | xTB 路径集成测试 |
| `tests/test_xtb_scan_input.py` | xTB 扫描输入测试 |
| `tests/test_xtb_scan_params.py` | xTB 扫描参数测试 |
| `tests/test_cleaner_adapter.py` | 清理适配器测试 |
| `tests/test_dataset_loader.py` | 数据集加载器测试 |
| `tests/test_fchk_reader_multiline.py` | FCHK 多行读取器测试 |
| `tests/test_orca_interface.py` | ORCA 接口测试 |
| `tests/test_gau_xtb_interface.py` | Gaussian-xTB 接口测试 |

### 测试优化 / Test Optimizations
- 归档重复测试到 `tests/deprecated/`
- 新增快速 CI gate 测试
- 改进 Mock QC 集成

---

## 配置更新 / Configuration Updates

### 新增配置项 / New Config Options

```yaml
# 反应配置文件 (Reaction Profiles)
reaction_profiles:
  "[4+3]_default":
    forming_bond_count: 2
    s2_strategy: forward_scan
    scan:
      scan_start_distance: 1.8
      scan_end_distance: 3.2
      scan_steps: 20
      scan_mode: concerted
      scan_force_constant: 0.5

# S1 特征提取 (S1 Feature Extraction)
step1:
  features:
    activation:
      boltzmann_weighted_energy: true
      conformer_entropy: true
      leaving_group_geometry: true
      alpha_h_gating: true

# Multiwfn 配置 (Multiwfn Configuration)
step4:
  multiwfn:
    enabled: true
    tiers:
      - fukui
      - dual_descriptor
      - qtaim
```

---

## 文档更新 / Documentation Updates

### 新增文档 / New Documents
| 文档 / Document | 描述 / Description |
|---------------|-------------------|
| `README.zh-CN.md` | 中文版 README |
| `docs/S2_S1_S2_2_Plan_20260315.md` | S2 规划文档 |
| `docs/S2_XTB_Path_Integration_Report.md` | xTB 路径集成报告 |
| `docs/S4_FEATURES_SUMMARY.md` | S4 特征总结 |
| `docs/S4_REPORT_COMPLETE.md` | S4 完成报告 |
| `docs/TESTING_GUIDE.md` | 测试指南 |
| `docs/QC_VALIDATION_GUIDE.md` | QC 验证指南 |
| `docs/DEV_STATUS_REPORT_20250311.md` | 开发状态报告 |
| `docs/S2_S3_Architecture_Report.md` | S2/S3 架构报告 |
| `docs/MODIFICATION_REPORT_20250311.md` | 修改报告 |
| `docs/DUAL_LEVEL_STRATEGY_SUMMARY.md` | 双层策略总结 |

---

## 依赖更新 / Dependency Updates

### Python 依赖 / Python Dependencies
- 保持 Python 3.8+ 兼容性
- 优化包管理配置

### QC 工具支持 / QC Tool Support
- Gaussian 16
- ORCA
- xTB
- CREST
- Multiwfn

---

## 已知问题 / Known Issues

1. **NBO 分析**: 默认禁用以减少运行时间，需要时可通过配置启用
2. **大分子构象搜索**: 可能需要较长计算时间
3. **并行化**: 当前版本主要支持顺序执行

---

## 贡献者 / Contributors

- QCcalc Team

---

**最后更新 / Last Updated**: 2025-03-18

**版本 / Version**: 2.0.0
