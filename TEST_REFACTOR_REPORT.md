# 测试脚本修整报告

## 概述

**当前状态**: 75个测试文件（含conftest.py）  
**目标状态**: 约20个核心测试文件  
**删除/合并**: 约55个文件

---

## 一、核心测试文件（保留）

### 1. CI快速门控（必须保留）

| 文件名 | 说明 | 优先级 |
|--------|------|--------|
| `conftest.py` | pytest配置 | P0 |
| `test_imports_step4_features.py` | 快速导入检查（55行） | P0 |
| `test_s4_no_qc_execution.py` | 验证S4不执行QC（6行） | P0 |

**推荐命令**:
```bash
pytest tests/test_imports_step4_features.py tests/test_s4_no_qc_execution.py -v
```

### 2. S4特性提取测试（核心功能）

| 文件名 | 说明 | 行数 | 优先级 |
|--------|------|------|--------|
| `test_s4_extractor_degrade_behavior.py` | 提取器降级行为测试 | 616 | P0 |
| `test_s4_v62_final_verification.py` | v6.2完整验证 | 393 | P0 |
| `test_s4_meta_warnings_and_weights.py` | 警告和权重策略 | 549 | P0 |
| `test_degradation_final.py` | 降级基础测试 | 76 | P0 |

**说明**: 这些测试覆盖了S4的核心功能：
- 降级行为（NaN填充而非跳过）
- 警告代码格式（W_前缀）
- 权重策略（qc.sample_weight）
- feature_meta.json结构

### 3. S4合同/契约测试

| 文件名 | 说明 | 优先级 |
|--------|------|--------|
| `test_s4_mech_packager.py` | 机制打包器 | P1 |
| `test_s4_artifact_integration.py` | 产物集成 | P1 |
| `test_s4_gedt_labeling.py` | GEDT标签 | P1 |
| `test_s4_path_compat.py` | 路径兼容 | P1 |
| `test_step4_cache_key.py` | 缓存键 | P1 |

### 4. M4机制模块测试

| 文件名 | 说明 | 优先级 |
|--------|------|--------|
| `test_m4_mech_index_contract.py` | mech_index契约 | P1 |
| `test_m4_mech_context_resolver.py` | 上下文解析 | P1 |
| `test_m4_qc_artifacts_structure.py` | QC产物结构 | P1 |
| `test_m4_qc_artifacts_mech_index.py` | QC产物索引 | P1 |
| `test_m4_template_structure.py` | 模板结构 | P1 |

### 5. M2模式版本测试

| 文件名 | 说明 | 优先级 |
|--------|------|--------|
| `test_m2_schema_versioning.py` | 模式版本 | P1 |
| `test_m2_step4_resume_semantics.py` | 恢复语义 | P1 |
| `test_m2_precursor_fallback.py` | 前体回退 | P1 |

### 6. M3 QC模拟测试

| 文件名 | 说明 | 优先级 |
|--------|------|--------|
| `test_m3_qc_mock_simple.py` | 简单QC模拟 | P1 |
| `test_m3_qc_collection_mock.py` | QC集合模拟 | P1 |
| `test_m3_gaussian_templates.py` | Gaussian模板 | P1 |

### 7. 简单集成测试

| 文件名 | 说明 | 优先级 |
|--------|------|--------|
| `test_mock_integration.py` | Mock集成（34行） | P1 |
| `test_mock_qc_e2e.py` | E2E mock测试（13行） | P1 |

### 8. 工具/沙箱测试

| 文件名 | 说明 | 优先级 |
|--------|------|--------|
| `test_sandbox_toxic_paths.py` | 毒性路径测试 | P2 |
| `test_step2_path_compat.py` | S2路径兼容 | P2 |
| `test_step3_path_compat.py` | S3路径兼容 | P2 |

### 9. S1/S2/S3基础测试

| 文件名 | 说明 | 优先级 |
|--------|------|--------|
| `test_two_stage_conformer.py` | 两阶段构象搜索 | P2 |
| `test_retro_scanner_v52.py` | 逆向扫描v5.2 | P2 |
| `test_s3_checkpoint.py` | S3检查点 | P2 |

---

## 二、冗余测试文件（删除）

### 1. Formchk重复测试（9个文件 → 合并为1个）

**待删除**:
```
test_formchk_unit.py
test_formchk_unit_tests.py
test_formchk_valid.py
test_formchk_working.py
test_formchk_final.py
test_try_formchk.py
test_try_formchk_standalone.py
test_try_formchk_validation.py
test_try_formchk_simplified.py
```

**合并方案**: 创建 `test_formchk_integration.py` (约150行)
- 保留所有测试用例的核心逻辑
- 统一使用pytest格式
- 删除重复用例

### 2. 降级测试重复（2个文件 → 保留1个）

**待删除**:
```
test_degradation.py  (93行，较旧)
```

**保留**:
```
test_degradation_final.py  (76行，更简洁)
```

### 3. 集成测试重复（3个文件 → 保留1个）

**待删除**:
```
test_integration.py        (127行，示例性质)
test_integration_final.py  (如果存在重复)
```

**保留**:
```
test_mock_integration.py   (34行，简洁有效)
test_mock_qc_e2e.py        (13行，补充)
```

### 4. 旧版本/实验性测试

**待删除**（未在当前架构中使用）:
```
test_stepwise_refactor.py
test_v2_2_pipeline.py
test_s4_step1_input_loop.py
test_update_reason_constants.py
test_e2e_precursor_leaving_group.py  (如果与test_two_stage_conformer重复)
```

### 5. 已弃用的FCHK解析测试

**待删除**:
```
test_fchk_reader_multiline.py  (如果已合并到主测试)
```

### 6. 工具类重复测试

检查以下文件是否有重复:
```
test_molecular_graph.py
test_molecule_utils.py
test_fragment_manipulation.py
test_sp_report.py
test_tsv_loader.py
test_geometry_preprocessor.py
test_fmo_cdft_dipolar.py
test_reference_states_runner_mock.py
test_reference_index_schema.py
```

### 7. 其他待确认文件

```
test_s4_path_compat.py  (与test_step2_path_compat重复?)
test_forward_scan_wiring.py
test_xtb_scan_input.py
test_orca_interface.py
test_qc_interface_v52.py
test_qctaskrunner_integration.py
test_thermo_validation_v52.py
```

---

## 三、推荐保留清单（22个文件）

```
# CI门控（3个）
tests/conftest.py
tests/test_imports_step4_features.py
tests/test_s4_no_qc_execution.py

# S4核心（6个）
tests/test_s4_extractor_degrade_behavior.py
tests/test_s4_v62_final_verification.py
tests/test_s4_meta_warnings_and_weights.py
tests/test_degradation_final.py
tests/test_mock_integration.py
tests/test_mock_qc_e2e.py

# S4合同（5个）
tests/test_s4_mech_packager.py
tests/test_s4_artifact_integration.py
tests/test_s4_gedt_labeling.py
tests/test_s4_path_compat.py
tests/test_step4_cache_key.py

# M4机制（5个）
tests/test_m4_mech_index_contract.py
tests/test_m4_mech_context_resolver.py
tests/test_m4_qc_artifacts_structure.py
tests/test_m4_qc_artifacts_mech_index.py
tests/test_m4_template_structure.py

# M2版本（3个）
tests/test_m2_schema_versioning.py
tests/test_m2_step4_resume_semantics.py
tests/test_m2_precursor_fallback.py
```

---

## 四、推荐删除清单（约53个文件）

### 高优先级删除（明显重复）
```
test_formchk_unit.py
test_formchk_unit_tests.py
test_formchk_valid.py
test_formchk_working.py
test_formchk_final.py
test_try_formchk.py
test_try_formchk_standalone.py
test_try_formchk_validation.py
test_try_formchk_simplified.py
test_degradation.py
test_integration.py
test_stepwise_refactor.py
test_v2_2_pipeline.py
test_s4_step1_input_loop.py
test_update_reason_constants.py
```

### 需要确认后删除（可能仍有价值）
```
test_m3_qc_mock_simple.py  (如果与test_m3_qc_collection_mock重复)
test_small_molecule_cache.py
test_small_molecule_catalog.py
test_gaussian_route_wrapping.py
test_thermo_validation_v52.py
test_fchk_reader_multiline.py
```

---

## 五、测试运行命令更新

### 快速CI（<10秒）
```bash
pytest tests/test_imports_step4_features.py tests/test_s4_no_qc_execution.py -v
```

### S4核心测试（<60秒）
```bash
pytest tests/test_s4_*.py tests/test_degradation_final.py -v
```

### 全部核心测试（<5分钟）
```bash
pytest tests/test_*.py -v --ignore=tests/tmp_v2_2_test/
```

---

## 六、实施建议

1. **第一步**: 备份整个tests目录
2. **第二步**: 删除明显重复的formchk测试（9个文件）
3. **第三步**: 删除test_degradation.py（保留test_degradation_final.py）
4. **第四步**: 删除test_integration.py（保留mock版本）
5. **第五步**: 删除旧版本/实验性测试
6. **第六步**: 运行测试验证
7. **第七步**: 更新tests/AGENTS.md

---

## 七、预期结果

- 测试文件数: 75 → 22
- 维护成本: 大幅降低
- 运行时间: 显著缩短
- 覆盖率: 保持核心功能覆盖
- 清晰度: 大幅提升

---

报告生成时间: 2026-03-11  
分析人: Claude Code
