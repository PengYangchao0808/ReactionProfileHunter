# ReactionProfileHunter V4.0.0

这是 RPH V4.0.0 的核心源码仓库。V4 只支持以下运行链：

```text
S0 mechanism -> S1 CENSO-LITE -> S2 PEB -> S3 low-level -> S4 high-level
```

项目根目录迁移后名称为 `RPH_V4.0.0/`。源码、配置、测试和文档均以 V4 合同为准；V3 归档位于仓库外部，不属于当前运行面。

## 目录结构

```text
RPH_V4.0.0/
├── rph_core/       # V4 编排器、S0-S4 阶段和 QC 接口
├── config/         # config/defaults.yaml：唯一配置源
├── tests/          # V4 协议、检查点、阶段和架构守卫测试
├── scripts/ci/     # 导入风格检查
├── docs/            # V4 设计、WSL 测试与清理记录
└── AGENTS.md       # 必须遵守的 V4 开发规范
```

## 开发验证

```bash
python -m py_compile rph_core/v4_orchestrator.py rph_core/steps/stage_calculator.py
python scripts/ci/check_imports.py rph_core
pytest -q tests/test_v4_protocol_contract.py tests/test_v4_checkpoint.py tests/test_v4_stage_calculator.py
```

WSL 测试计划见 [docs/WSL_TEST_PLAN_V4.md](docs/WSL_TEST_PLAN_V4.md)。计算输出、缓存和本地备份不会进入 Git；V3 历史说明保留在 `docs/ARCHIVE_V3.md` 等归档文档中。
