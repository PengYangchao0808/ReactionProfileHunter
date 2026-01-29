# Import Style Guidelines

## Purpose

These guidelines ensure consistent, maintainable import statements across the ReactionProfileHunter codebase. Following these conventions prevents import chain errors and improves code clarity.

## Core Principles

### 1. Cross-Package/Module Imports: Use Absolute Paths

When importing from `rph_core/utils` (the shared utilities module), **always use absolute imports**:

```python
# ✅ CORRECT - Absolute import
from rph_core.utils.file_io import read_xyz
from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.qc_interface import XTBInterface

# ❌ INCORRECT - Multi-dot relative import
from ...utils.file_io import read_xyz  # Resolves to rph_core.steps.utils (doesn't exist!)
from ....utils.file_io import read_xyz
```

**Why?**
- `rph_core/utils/` is at the top level under `rph_core/`, not under `steps/`
- From subdirectories like `steps/step3_opt/` or `steps/step4_features/`, the relative path would be:
  - 3 dots (`...`) → `rph_core/steps/` (wrong directory)
  - 4 dots (`....`) → `rph_core/` (correct, but fragile)
- Absolute imports are unambiguous and survive directory restructure

### 2. Same-Subpackage Imports: Use Short Relative Imports

For imports within the same subpackage (e.g., extractors importing from each other), use short relative imports:

```python
# ✅ CORRECT - Same subpackage (extractors/)
from .base import BaseExtractor, register_extractor
from . import thermo, geometry, qc_checks

# ❌ INCORRECT - Absolute import for same subpackage
from rph_core.steps.step4_features.extractors.base import BaseExtractor
```

**Why?**
- Single-dot imports are clear and concise
- Maintains encapsulation within the subpackage
- Moves the subpackage as a unit

### 3. Import Order (PEP 8 Compliant)

Order imports as follows:

```python
# 1. Standard library imports
import logging
from pathlib import Path
from typing import Optional, Dict, Any

# 2. Third-party imports
import numpy as np
from rdkit import Chem

# 3. Local imports (absolute preferred)
from rph_core.utils.log_manager import LoggerMixin
from rph_core.steps.step4_features.context import FeatureContext
from .base import BaseExtractor  # Same-subpackage relative
```

## Common Pitfalls

### Pitfall 1: Multi-Dot Relative Imports to Wrong Directory

```python
# ❌ WRONG - from rph_core/steps/step3_opt/
from ...utils.file_io import read_xyz
# Resolves to: rph_core.steps.utils
# Actual location: rph_core/utils
# Result: ModuleNotFoundError
```

**Fix:** Use absolute import
```python
# ✅ CORRECT
from rph_core.utils.file_io import read_xyz
```

### Pitfall 2: Inconsistent Mixing of Styles

```python
# ❌ INCONSISTENT - Mix of absolute and multi-dot relative
from rph_core.utils.log_manager import LoggerMixin
from ...utils.file_io import read_xyz  # Different style!
```

**Fix:** Be consistent - prefer absolute for utils
```python
# ✅ CONSISTENT
from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.file_io import read_xyz
```

### Pitfall 3: Assuming Directory Structure

```python
# ❌ FRAGILE - Assumes 4 dots reach utils
from ....utils.file_io import read_xyz  # Breaks if directory depth changes
```

**Fix:** Use absolute import (survives restructure)
```python
# ✅ ROBUST
from rph_core.utils.file_io import read_xyz
```

## Testing Import Chains

Run import smoke tests before running full pipeline:

```bash
# Run all import smoke tests
pytest tests/test_imports_step4_features.py -v

# Quick manual check
python -c "import rph_core.steps.step4_features.feature_miner as m; print('OK')"
```

## Migration from Relative to Absolute

### Case 1: From `steps/step3_opt/` or `steps/step4_features/`

**Before:**
```python
from ...utils.log_manager import LoggerMixin
from ...utils.file_io import write_xyz, read_xyz
from ...utils.orca_interface import ORCAInterface
```

**After:**
```python
from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.file_io import write_xyz, read_xyz
from rph_core.utils.orca_interface import ORCAInterface
```

### Case 2: From `steps/step4_features/extractors/`

**Before:**
```python
from ...utils.file_io import read_xyz
```

**After:**
```python
from rph_core.utils.file_io import read_xyz
```

## Enforcement

### CI/CD Checks

1. **Import Smoke Tests** - Run on every commit
   ```bash
   pytest tests/test_imports_*.py
   ```

2. **CI Import Style Check** - Block multi-dot relative imports
   ```bash
   # Run the automated import checker
   python ci/check_imports.py rph_core
   ```

   The `ci/check_imports.py` script scans all Python files and blocks:
   - `from ...utils.*` patterns (resolves to incorrect `rph_core.steps.utils`)
   - `from ....utils.*` patterns (ambiguous, should use absolute)

   **Exit codes:**
   - 0: Pass - No violations
   - 1: Fail - Forbidden patterns found

   **Why not use Ruff F401?**
   `F401` (unused imports) does NOT catch multi-dot relative imports. This creates a false sense of security. Use the dedicated `ci/check_imports.py` script instead, which uses grep-based pattern matching tailored to your project structure.

3. **Static Analysis (Future Enhancement)** - Consider adding linter rules:
   - Flag multi-dot relative imports (`from ...` or `from ....`)
   - Enforce absolute imports for utils
   - Tools: `ruff`, `flake8` with custom plugins, or `pylint`

### Code Review Checklist

When reviewing PRs:
- [ ] No `from ...utils` imports (should be `from rph_core.utils`)
- [ ] No `from ....utils` imports
- [ ] Same-subpackage imports use single-dot (`.module`)
- [ ] Import order follows PEP 8 (std → third-party → local)
- [ ] Import smoke tests pass

## References

- [PEP 8 - Imports](https://peps.python.org/pep-0008/#imports)
- [Python Module Import System](https://docs.python.org/3/tutorial/modules.html)

## History

- **2026-01-19**: Created guidelines after fixing `ModuleNotFoundError` in Step 4
- Issue: Multi-dot relative imports causing `rph_core.steps.utils` ghost module error
- Solution: Standardized on absolute imports for `rph_core/utils`
