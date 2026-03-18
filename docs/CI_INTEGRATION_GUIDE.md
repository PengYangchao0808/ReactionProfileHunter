# CI Integration Guide

## Import Style Checks in CI Pipeline

This document explains how to integrate import style checks into CI/CD workflows.

---

## 1. Import Smoke Tests (pytest)

Run on every commit and pull request to catch import chain errors early.

### Execution
```bash
# Run all import smoke tests
pytest tests/test_imports_*.py -v
```

### GitHub Actions Example
```yaml
name: CI

on: [push, pull_request]

jobs:
  import-checks:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.12'
      - name: Install dependencies
        run: |
          pip install -e .
          pip install pytest
      - name: Run import smoke tests
        run: |
          pytest tests/test_imports_*.py -v
```

---

## 2. CI Import Style Checker (grep-based script)

Blocks multi-dot relative imports that would cause ModuleNotFoundError.

### Execution
```bash
# Check all Python files in rph_core/
python scripts/ci/check_imports.py rph_core
```

### Exit Codes
- `0`: Pass - No forbidden import patterns found
- `1`: Fail - Multi-dot relative imports detected

### What It Checks

The script scans for these forbidden patterns:

| Pattern | Why Forbidden | Correct Fix |
|---------|---------------|--------------|
| `from ...utils` | Resolves to `rph_core.steps.utils` (doesn't exist) | `from rph_core.utils` |
| `from ....utils` | Ambiguous depth, breaks on restructure | `from rph_core.utils` |

### GitHub Actions Example
```yaml
name: CI

on: [push, pull_request]

jobs:
  import-style-check:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.12'
      - name: Run CI import style checker
        run: |
          python scripts/ci/check_imports.py rph_core
```

### Pre-commit Hook (Optional)

Add to `.git/hooks/pre-commit` (or use pre-commit framework):

```bash
#!/bin/bash
# .git/hooks/pre-commit

echo "Running import style checks..."

# Check import style
python scripts/ci/check_imports.py rph_core
if [ $? -ne 0 ]; then
    echo "❌ Import style check failed"
    exit 1
fi

# Run smoke tests
pytest tests/test_imports_*.py -q
if [ $? -ne 0 ]; then
    echo "❌ Import smoke tests failed"
    exit 1
fi

exit 0
```

Make executable:
```bash
chmod +x .git/hooks/pre-commit
```

---

## 3. CI Pipeline Example (Complete)

Combine both checks for comprehensive coverage:

```yaml
name: CI

on: [push, pull_request]

jobs:
  import-checks:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.12'
      - name: Install dependencies
        run: |
          pip install -e .
          pip install pytest
      - name: Run import smoke tests
        run: |
          pytest tests/test_imports_*.py -v
      - name: Run CI import style checker
        run: |
          python scripts/ci/check_imports.py rph_core
```

---

## Verification

### Test the CI Checker Manually

```bash
# Should PASS (no violations)
python scripts/ci/check_imports.py rph_core
echo $?  # Expected: 0

# Create test file with violations
cat > /tmp/test_violation.py << 'EOF'
from ...utils.file_io import read_xyz
EOF

# Should FAIL (detects violations)
python scripts/ci/check_imports.py /tmp 2>&1 | head -20
echo $?  # Expected: 1
```

### Expected Output

**When PASS:**
```
✅ No forbidden import patterns found
📁 Scanned 59 Python files
✅ PASSED: Import style check
```

**When FAIL:**
```
❌ IMPORT VIOLATIONS FOUND
================================================================================

📄 rph_core/steps/step4_features/extractors/geometry.py
   Line 15: from ...utils.file_io
   ⚠️  from ...utils (resolves to rph_core.steps.utils)
   ✅ FIX: Use 'from rph_core.utils...' instead

================================================================================

📊 SUMMARY: 1 violation(s) in 1 file(s)

❌ FAILED: Multi-dot relative imports detected
```

---

## Troubleshooting

### Issue: "Import style checker fails but code is correct"

**Possible cause:** False positive due to legitimate use case.

**Solution:** Add to `ALLOWED_EXCEPTIONS` list in `scripts/ci/check_imports.py`:

```python
ALLOWED_EXCEPTIONS = [
    # Add file:line pairs with justification
    "rph_core/special_module.py:123",  # Legitimate use of from ...utils
]
```

---

## Related Documentation

- [tests/test_imports_step4_features.py](../tests/test_imports_step4_features.py) - Import smoke tests

---

## History

- **2026-01-19**: Created CI integration guide with both pytest smoke tests and grep-based checker
- **2026-03-11**: Moved from `ci/` to `docs/`; script path updated to `scripts/ci/check_imports.py`
