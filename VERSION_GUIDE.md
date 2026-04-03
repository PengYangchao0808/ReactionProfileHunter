# ReactionProfileHunter 版本管理指南

## Version Management Guide

**版本**: 2.0.0  
**生效日期**: 2025-03-18  
**文档语言**: 中文 (English version available: VERSION_GUIDE_EN.md)

---

## 目录 / Table of Contents

1. [版本命名规范](#1-版本命名规范)
2. [分支管理策略](#2-分支管理策略)
3. [发布流程](#3-发布流程)
4. [变更日志规范](#4-变更日志规范)
5. [版本号更新 Checklist](#5-版本号更新-checklist)
6. [Git 标签管理](#6-git-标签管理)
7. [向后兼容性策略](#7-向后兼容性策略)

---

## 1. 版本命名规范

### 1.1 版本格式

采用 **语义化版本 (Semantic Versioning)** 格式：

```
MAJOR.MINOR.PATCH[-PRERELEASE][+BUILD]
```

| 组件 | 说明 | 示例 |
|------|------|------|
| MAJOR | 主版本 - 不兼容的 API 变更 | 2.0.0 |
| MINOR | 次版本 - 向后兼容的新功能 | 2.1.0 |
| PATCH | 补丁版本 - 向后兼容的 bug 修复 | 2.1.1 |
| PRERELEASE | 预发布版本 | 2.0.0-alpha.1 |
| BUILD | 构建元数据 | 2.0.0+20250318 |

### 1.2 版本阶段定义

| 阶段 | 描述 | 命名示例 |
|------|------|----------|
| **开发中 (Development)** | 新功能开发中 | 2.0.0-dev |
| **测试版 (Beta)** | 公开测试 | 2.0.0-beta.1 |
| **候选版 (RC)** | 发布候选 | 2.0.0-rc.1 |
| **正式版 (Release)** | 稳定发布 | 2.0.0 |
| **长期支持 (LTS)** | 长期维护 | 2.0.0-lts |

### 1.3 当前版本规划

```
v2.0.0 (当前开发版本)
    │
    ├── v2.0.1-dev (后续开发)
    ├── v2.0.0-lts (长期支持)
    └── v2.1.0 (下一功能版本)
```

---

## 2. 分支管理策略

### 2.1 分支模型: Git Flow

采用 **Git Flow** 分支模型：

```
main (正式发布分支)
│
├── develop (开发主分支)
│   │
│   ├── feature/* (功能分支)
│   │   └── feature/forward-scan-optimization
│   │
│   ├── bugfix/* (热修复分支)
│   │   └── bugfix/fix-xtb-scan-error
│   │
│   └── release/* (发布分支)
│       └── release/v2.0.0
```

### 2.2 分支说明

| 分支类型 | 命名规则 | 生命周期 | 合并目标 |
|----------|----------|----------|----------|
| `main` | 固定 | 永久 | - |
| `develop` | 固定 | 永久 | main (发布时) |
| `feature/*` | feature/描述 | 数天-数周 | develop |
| `bugfix/*` | bugfix/描述 | 数小时-数天 | develop 或 release |
| `hotfix/*` | hotfix/描述 | 数小时-数天 | main + develop |
| `release/*` | release/v版本 | 数天-数周 | main + develop |

### 2.3 当前分支状态

```
main (v1.x 历史版本)
    │
    └── v2_4+3_basic (当前开发分支)
            │
            ├── feature-login
            └── fix/m4-repair
```

### 2.4 建议分支重组

建议将当前开发分支整合为标准的 Git Flow：

```bash
# 创建 develop 分支 (从当前开发分支)
git checkout -b develop

# main 保持历史版本
# develop 用于日常开发
# feature/* 用于新功能
```

---

## 3. 发布流程

### 3.1 发布流程图

```
┌─────────────────────────────────────────────────────────────┐
│                      发布流程 Release Flow                   │
└─────────────────────────────────────────────────────────────┘

1. 准备发布
   │
   ├── git checkout -b release/v2.0.0 develop
   │
   ├── 更新版本号 (version.py)
   │
   ├── 更新 CHANGELOG.md
   │
   └── 运行完整测试 ✓
           │
           ▼
2. 发布审查
   │
   ├── 代码审查 (Code Review)
   │
   ├── 文档完整性检查
   │
   └── 测试覆盖率验证
           │
           ▼
3. 正式发布
   │
   ├── git merge release/v2.0.0 main
   │
   ├── git tag -a v2.0.0 -m "Release v2.0.0"
   │
   ├── git push main --tags
   │
   └── git merge main develop
           │
           ▼
4. 发布后处理
   │
   ├── 创建 GitHub Release
   │
   ├── 更新文档网站 (如有)
   │
   └── 通知相关人员
```

### 3.2 发布 Checklist

- [ ] 所有测试通过 (`pytest -v`)
- [ ] 代码格式检查通过 (`black`, `isort`)
- [ ] 类型检查通过 (`mypy`)
- [ ] 更新 `rph_core/version.py` 版本号
- [ ] 更新 `CHANGELOG.md` 变更日志
- [ ] 更新 `README.md` 版本徽章
- [ ] 创建 git tag
- [ ] 推送 tags 到远程
- [ ] 合并到 main 和 develop
- [ ] 创建 GitHub Release (如使用 GitHub)

### 3.3 发布命令示例

```bash
# 1. 创建发布分支
git checkout -b release/v2.0.0 develop

# 2. 更新版本号
# 编辑 rph_core/version.py
vim rph_core/version.py
# __version__ = "2.0.0"

# 3. 更新 CHANGELOG
vim CHANGELOG.md
# 添加 ## v2.0.0 (2025-03-18) 章节

# 4. 提交更改
git add .
git commit -m "release: prepare v2.0.0"

# 5. 合并到 main
git checkout main
git merge release/v2.0.0

# 6. 创建标签
git tag -a v2.0.0 -m "Release v2.0.0: 2.0 major release with forward-scan"

# 7. 推送
git push main v2.0.0

# 8. 合并回 develop
git checkout develop
git merge main

# 9. 清理发布分支
git branch -d release/v2.0.0
```

---

## 4. 变更日志规范

### 4.1 变更日志格式

采用 **Keep a Changelog** 格式：

```markdown
## [2.0.0] - 2025-03-18

### Added
- 新功能描述

### Changed
- 现有功能变更

### Deprecated
- 即将弃用的功能

### Removed
- 已移除的功能

### Fixed
- Bug 修复

### Security
- 安全相关修复
```

### 4.2 提交信息规范

使用 **Conventional Commits** 规范：

```
<type>(<scope>): <description>

[optional body]

[optional footer]
```

| Type | Description |
|------|-------------|
| `feat` | 新功能 |
| `fix` | Bug 修复 |
| `docs` | 文档更新 |
| `style` | 代码格式 |
| `refactor` | 重构 |
| `test` | 测试相关 |
| `chore` | 构建/工具 |

示例：
```bash
git commit -m "feat(s2): add forward-scan for [4+3] cycloaddition"
git commit -m "fix(xtb): resolve scan parameter validation error"
git commit -m "docs: update CHANGELOG for v2.0.0"
```

---

## 5. 版本号更新 Checklist

### 5.1 主版本更新 (MAJOR - 2.0.0)

- [ ] 重大架构变更
- [ ] 不兼容的 API 变更
- [ ] 移除已弃用的功能
- [ ] 项目重命名
- [ ] 重大重构

### 5.2 次版本更新 (MINOR - x.1.0)

- [ ] 新增功能 (向后兼容)
- [ ] 新增可选参数
- [ ] 新增配置文件选项
- [ ] 性能优化
- [ ] 新增支持的反应类型

### 5.3 补丁版本更新 (PATCH - x.x.1)

- [ ] Bug 修复
- [ ] 文档修正
- [ ] 类型注解修正
- [ ] 测试补充
- [ ] 小优化

---

## 6. Git 标签管理

### 6.1 标签命名规范

| 标签类型 | 命名示例 | 说明 |
|----------|----------|------|
| 正式版 | `v2.0.0` | 正式发布版本 |
| 预发布 | `v2.0.0-beta.1` | 测试版 |
| RC | `v2.0.0-rc.1` | 发布候选 |
| 开发版 | `v2.0.0-dev` | 开发快照 |

### 6.2 标签操作命令

```bash
# 创建轻量标签
git tag v2.0.0

# 创建附注标签 (推荐)
git tag -a v2.0.0 -m "Release v2.0.0"

# 列出所有标签
git tag -l

# 查看标签信息
git show v2.0.0

# 删除本地标签
git tag -d v2.0.0

# 删除远程标签
git push origin :refs/tags/v2.0.0

# 推送单个标签
git push origin v2.0.0

# 推送所有标签
git push origin --tags
```

### 6.3 当前建议标签结构

```
v2.0.0    ← 当前开发版本 (即将发布)
v1.0.0    ← 历史版本 (如保留)
```

---

## 7. 向后兼容性策略

### 7.1 兼容性级别

| 级别 | 说明 | 变更规则 |
|------|------|----------|
| **完全兼容** | 现有代码无需修改 | MINOR/PATCH |
| **警告兼容** | 现有代码会产生警告 | MINOR (添加警告) |
| **破坏兼容** | 需要修改现有代码 | MAJOR |

### 7.2 兼容性保障措施

1. **废弃预告 (Deprecation)**
   - 至少保留一个 MINOR 版本
   - 提供替代方案文档
   - 产生弃用警告

2. **配置文件变更**
   - 新版本需兼容旧配置文件
   - 提供配置迁移指南
   - 保留默认值向后兼容

3. **输出文件格式**
   - 保持输出格式向后兼容
   - 如需变更，在 CHANGELOG 中说明

### 7.3 弃用周期示例

```
v2.0.0: 添加新参数 'new_param' (可选)
    │
    ├── v2.1.0: 标记旧参数 'old_param' 为 deprecated
    │           (产生警告，但仍可用)
    │
    └── v3.0.0: 移除 'old_param'
```

---

## 附录 A: 快速参考命令

### A.1 版本更新快速命令

```bash
# 完整发布流程
make release VERSION=2.0.0

# 或手动执行
./scripts/release.sh 2.0.0

# 仅更新版本号
bumpversion patch  # 2.0.0 -> 2.0.1
bumpversion minor  # 2.0.0 -> 2.1.0
bumpversion major  # 2.0.0 -> 3.0.0
```

### A.2 版本检查

```bash
# 查看当前版本
python -m rph_core --version

# 或
cat rph_core/version.py

# 查看最近标签
git describe --tags --abbrev=0
```

---

## 附录 B: 版本历史

| 版本 | 日期 | 说明 |
|------|------|------|
| 2.0.0 | 开发中 | 2.0 重大版本 - 前向扫描、重构 |
| 1.0.0 | 历史 | 初始版本 (v6.x 时代) |

---

## 附录 C: 相关文件

- `rph_core/version.py` - 版本号定义
- `CHANGELOG.md` - 变更日志
- `CHANGELOG_EN.md` - 英文变更日志
- `README.md` / `README.zh-CN.md` - 项目文档

---

**文档维护者**: ReactionProfileHunter Team  
**下次评审**: 2025-06-01  
**版本**: 1.0
