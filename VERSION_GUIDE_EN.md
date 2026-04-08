# ReactionProfileHunter Version Management Guide

**Version**: 2.0.0  
**Effective Date**: 2025-03-18

---

## Table of Contents

1. [Version Naming Convention](#1-version-naming-convention)
2. [Branch Management Strategy](#2-branch-management-strategy)
3. [Release Process](#3-release-process)
4. [Changelog Convention](#4-changelog-convention)
5. [Version Update Checklist](#5-version-update-checklist)
6. [Git Tag Management](#6-git-tag-management)
7. [Backward Compatibility Strategy](#7-backward-compatibility-strategy)

---

## 1. Version Naming Convention

### 1.1 Version Format

采用 **语义化版本 (Semantic Versioning)** 格式：

```
MAJOR.MINOR.PATCH[-PRERELEASE][+BUILD]
```

| Component | Description | Example |
|-----------|-------------|---------|
| MAJOR | Incompatible API changes | 2.0.0 |
| MINOR | Backward-compatible new features | 2.1.0 |
| PATCH | Backward-compatible bug fixes | 2.1.1 |
| PRERELEASE | Pre-release versions | 2.0.0-alpha.1 |
| BUILD | Build metadata | 2.0.0+20250318 |

### 1.2 Version Stages

| Stage | Description | Example |
|-------|-------------|---------|
| **Development** | Feature development in progress | 2.0.0-dev |
| **Beta** | Public testing | 2.0.0-beta.1 |
| **RC** | Release Candidate | 2.0.0-rc.1 |
| **Release** | Stable release | 2.0.0 |
| **LTS** | Long Term Support | 2.0.0-lts |

### 1.3 Current Version Planning

```
v2.0.0 (Current Development Version)
    │
    ├── v2.0.1-dev (Subsequent Development)
    ├── v2.0.0-lts (Long Term Support)
    └── v2.1.0 (Next Feature Release)
```

---

## 2. Branch Management Strategy

### 2.1 Branch Model: Git Flow

```
main (Production Release Branch)
│
├── develop (Development Main Branch)
│   │
│   ├── feature/* (Feature Branches)
│   │   └── feature/forward-scan-optimization
│   │
│   ├── bugfix/* (Bug Fix Branches)
│   │   └── bugfix/fix-xtb-scan-error
│   │
│   └── release/* (Release Branches)
│       └── release/v2.0.0
```

### 2.2 Branch Types

| Branch Type | Naming Rule | Lifecycle | Merge Target |
|-------------|-------------|----------|--------------|
| `main` | Fixed | Permanent | - |
| `develop` | Fixed | Permanent | main (on release) |
| `feature/*` | feature/description | Days-Weeks | develop |
| `bugfix/*` | bugfix/description | Hours-Days | develop or release |
| `hotfix/*` | hotfix/description | Hours-Days | main + develop |
| `release/*` | release/vversion | Days-Weeks | main + develop |

### 2.3 Current Branch Status

```
main (v1.x Historical Version)
    │
    └── v2_4+3_basic (Current Development Branch)
            │
            ├── feature-login
            └── fix/m4-repair
```

### 2.4 Recommended Branch Reorganization

```bash
# Create develop branch (from current development branch)
git checkout -b develop

# main keeps historical versions
# develop for daily development
# feature/* for new features
```

---

## 3. Release Process

### 3.1 Release Flow

```
┌─────────────────────────────────────────────────────────────┐
│                      Release Flow                            │
└─────────────────────────────────────────────────────────────┘

1. Prepare Release
   │
   ├── git checkout -b release/v2.0.0 develop
   │
   ├── Update version (version.py)
   │
   ├── Update CHANGELOG.md
   │
   └── Run full test suite ✓
           │
           ▼
2. Release Review
   │
   ├── Code Review
   │
   ├── Documentation Completeness Check
   │
   └── Test Coverage Verification
           │
           ▼
3. Official Release
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
4. Post-Release
   │
   ├── Create GitHub Release
   │
   ├── Update Documentation Site (if applicable)
   │
   └── Notify Stakeholders
```

### 3.2 Release Checklist

- [ ] All tests pass (`pytest -v`)
- [ ] Code format check passes (`black`, `isort`)
- [ ] Type check passes (`mypy`)
- [ ] Update `rph_core/version.py` version number
- [ ] Update `CHANGELOG.md` changelog
- [ ] Update `README.md` version badges
- [ ] Create git tag
- [ ] Push tags to remote
- [ ] Merge to main and develop
- [ ] Create GitHub Release (if using GitHub)

### 3.3 Release Commands

```bash
# 1. Create release branch
git checkout -b release/v2.0.0 develop

# 2. Update version number
# Edit rph_core/version.py
vim rph_core/version.py
# __version__ = "2.0.0"

# 3. Update CHANGELOG
vim CHANGELOG.md
# Add ## v2.0.0 (2025-03-18) section

# 4. Commit changes
git add .
git commit -m "release: prepare v2.0.0"

# 5. Merge to main
git checkout main
git merge release/v2.0.0

# 6. Create tag
git tag -a v2.0.0 -m "Release v2.0.0: 2.0 major release with forward-scan"

# 7. Push
git push main v2.0.0

# 8. Merge back to develop
git checkout develop
git merge main

# 9. Clean up release branch
git branch -d release/v2.0.0
```

---

## 4. Changelog Convention

### 4.1 Changelog Format

Using **Keep a Changelog** format:

```markdown
## [2.0.0] - 2025-03-18

### Added
- New feature description

### Changed
- Existing feature changes

### Deprecated
- Upcoming deprecations

### Removed
- Removed features

### Fixed
- Bug fixes

### Security
- Security fixes
```

### 4.2 Commit Message Convention

Using **Conventional Commits**:

```
<type>(<scope>): <description>

[optional body]

[optional footer]
```

| Type | Description |
|------|-------------|
| `feat` | New feature |
| `fix` | Bug fix |
| `docs` | Documentation |
| `style` | Code formatting |
| `refactor` | Refactoring |
| `test` | Test-related |
| `chore` | Build/Tools |

Examples:
```bash
git commit -m "feat(s2): add forward-scan for [4+3] cycloaddition"
git commit -m "fix(xtb): resolve scan parameter validation error"
git commit -m "docs: update CHANGELOG for v2.0.0"
```

---

## 5. Version Update Checklist

### 5.1 Major Version Update (MAJOR - 2.0.0)

- [ ] Major architectural changes
- [ ] Incompatible API changes
- [ ] Remove deprecated features
- [ ] Project rename
- [ ] Major refactoring

### 5.2 Minor Version Update (MINOR - x.1.0)

- [ ] Add new features (backward compatible)
- [ ] Add optional parameters
- [ ] Add new config file options
- [ ] Performance improvements
- [ ] New supported reaction types

### 5.3 Patch Version Update (PATCH - x.x.1)

- [ ] Bug fixes
- [ ] Documentation corrections
- [ ] Type annotation fixes
- [ ] Test additions
- [ ] Minor optimizations

---

## 6. Git Tag Management

### 6.1 Tag Naming Convention

| Tag Type | Example | Description |
|----------|---------|-------------|
| Release | `v2.0.0` | Official release version |
| Pre-release | `v2.0.0-beta.1` | Beta test version |
| RC | `v2.0.0-rc.1` | Release Candidate |
| Dev | `v2.0.0-dev` | Development snapshot |

### 6.2 Tag Commands

```bash
# Create lightweight tag
git tag v2.0.0

# Create annotated tag (recommended)
git tag -a v2.0.0 -m "Release v2.0.0"

# List all tags
git tag -l

# Show tag info
git show v2.0.0

# Delete local tag
git tag -d v2.0.0

# Delete remote tag
git push origin :refs/tags/v2.0.0

# Push single tag
git push origin v2.0.0

# Push all tags
git push origin --tags
```

### 6.3 Recommended Tag Structure

```
v2.0.0    ← Current development version (about to release)
v1.0.0    ← Historical version (if kept)
```

---

## 7. Backward Compatibility Strategy

### 7.1 Compatibility Levels

| Level | Description | Change Rule |
|-------|-------------|-------------|
| **Full Compatible** | No code changes needed | MINOR/PATCH |
| **Warning Compatible** | Warnings generated | MINOR (with warnings) |
| **Breaking Compatible** | Code changes required | MAJOR |

### 7.2 Compatibility Measures

1. **Deprecation**
   - Keep at least one MINOR version
   - Provide alternative documentation
   - Generate deprecation warnings

2. **Config File Changes**
   - New version must be compatible with old configs
   - Provide migration guide
   - Keep default values backward compatible

3. **Output File Format**
   - Keep output format backward compatible
   - Document changes in CHANGELOG if changed

### 7.3 Deprecation Cycle Example

```
v2.0.0: Add new parameter 'new_param' (optional)
    │
    ├── v2.1.0: Mark old parameter 'old_param' as deprecated
    │           (Generates warning, still usable)
    │
    └── v3.0.0: Remove 'old_param'
```

---

## Appendix A: Quick Reference Commands

### A.1 Version Update Quick Commands

```bash
# Full release process
make release VERSION=2.0.0

# Or manually
./scripts/release.sh 2.0.0

# Only update version number
bumpversion patch  # 2.0.0 -> 2.0.1
bumpversion minor  # 2.0.0 -> 2.1.0
bumpversion major  # 2.0.0 -> 3.0.0
```

### A.2 Version Check

```bash
# View current version
python -m rph_core --version

# Or
cat rph_core/version.py

# View latest tag
git describe --tags --abbrev=0
```

---

## Appendix B: Version History

| Version | Date | Description |
|---------|------|-------------|
| 2.0.0 | In Development | 2.0 Major Release - Forward-scan, Refactoring |
| 1.0.0 | Historical | Initial Version (v6.x Era) |

---

## Appendix C: Related Files

- `rph_core/version.py` - Version number definition
- `CHANGELOG.md` - Changelog
- `CHANGELOG_EN.md` - English Changelog
- `README.md` / `README.zh-CN.md` - Project documentation

---

**Document Maintainer**: ReactionProfileHunter Team  
**Next Review**: 2025-06-01  
**Version**: 1.0
