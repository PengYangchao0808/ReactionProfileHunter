#!/usr/bin/env python3
"""
CI Check Script: Block Multi-Dot Relative Imports to rph_core.utils
==================================================================

This script scans Python files for dangerous multi-dot relative import patterns
that would resolve to incorrect paths like 'rph_core.steps.utils'.

Purpose: Prevent ModuleNotFoundError before code enters CI pipeline.
Execution: Fast (grep-based, < 5 seconds for entire codebase).
Exit Codes:
    0: No violations found (PASS)
    1: Forbidden import patterns found (FAIL)

Usage:
    python ci/check_imports.py [directory]

    Default directory: rph_core/

Author: QC Descriptors Team
Date: 2026-01-19
"""

import sys
import os
import re
from pathlib import Path
from typing import List, Tuple


# Forbidden patterns: Multi-dot imports targeting utils
# These are dangerous because they resolve to wrong directory structure
FORBIDDEN_PATTERNS = [
    (r'from \.\.\..*utils', "from ...utils (resolves to rph_core.steps.utils)"),
    (r'from \.\.\.\..*utils', "from ....utils (ambiguous, should use absolute)"),
]

# Allowed exceptions (documented here to prevent false positives)
# Example: If there's a legitimate use case, add it here with rationale
ALLOWED_EXCEPTIONS = []


def find_python_files(root_dir: Path) -> List[Path]:
    """Find all Python files in directory recursively."""
    return list(root_dir.rglob("*.py"))


def check_file_for_violations(file_path: Path) -> List[Tuple[int, str, str]]:
    """
    Check a single Python file for forbidden import patterns.

    Returns:
        List of tuples: (line_number, matched_text, violation_description)
    """
    violations = []

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                # Skip comments
                stripped = line.strip()
                if stripped.startswith('#'):
                    continue

                for pattern, description in FORBIDDEN_PATTERNS:
                    match = re.search(pattern, stripped)
                    if match:
                        # Check if this is an allowed exception
                        exception_key = f"{file_path}:{line_num}"
                        if exception_key in ALLOWED_EXCEPTIONS:
                            continue

                        violations.append((line_num, match.group(0), description))

    except Exception as e:
        print(f"ERROR: Could not read {file_path}: {e}", file=sys.stderr)

    return violations


def check_imports(directory: Path) -> Tuple[int, int]:
    """
    Scan all Python files in directory for import violations.

    Returns:
        (total_violations, total_files_with_violations)
    """
    python_files = find_python_files(directory)
    all_violations = []

    for py_file in python_files:
        # Skip backup files and test files (optional, adjust as needed)
        if 'backup' in py_file.name.lower():
            continue
        if 'test_' in py_file.name and py_file.name != 'check_imports.py':
            continue

        violations = check_file_for_violations(py_file)
        if violations:
            all_violations.append((py_file, violations))

    # Report results
    if all_violations:
        print("❌ IMPORT VIOLATIONS FOUND", file=sys.stderr)
        print("=" * 80, file=sys.stderr)

        for file_path, violations in all_violations:
            rel_path = file_path.relative_to(directory.parent)
            print(f"\n📄 {rel_path}", file=sys.stderr)

            for line_num, matched_text, description in violations:
                print(
                    f"   Line {line_num}: {matched_text}",
                    file=sys.stderr
                )
                print(
                    f"   ⚠️  {description}",
                    file=sys.stderr
                )

                # Provide fix suggestion
                print(
                    f"   ✅ FIX: Use 'from rph_core.utils...' instead",
                    file=sys.stderr
                )

        print("\n" + "=" * 80, file=sys.stderr)
        total_violations = sum(len(v) for _, v in all_violations)
        total_files = len(all_violations)
        print(
            f"\n📊 SUMMARY: {total_violations} violation(s) in {total_files} file(s)",
            file=sys.stderr
        )
        return total_violations, total_files

    else:
        print("✅ No forbidden import patterns found", file=sys.stdout)
        print(f"📁 Scanned {len(python_files)} Python files", file=sys.stdout)
        return 0, 0


def main():
    # Determine directory to scan
    if len(sys.argv) > 1:
        target_dir = Path(sys.argv[1])
    else:
        # Default: scan rph_core/ directory
        target_dir = Path(__file__).parent.parent / "rph_core"

    if not target_dir.exists():
        print(f"❌ ERROR: Directory not found: {target_dir}", file=sys.stderr)
        sys.exit(2)

    # Run checks
    total_violations, total_files = check_imports(target_dir)

    # Exit with appropriate code
    if total_violations > 0:
        print("\n❌ FAILED: Multi-dot relative imports detected", file=sys.stderr)
        print(
            "Please follow IMPORT_GUIDELINES.md and use absolute imports for rph_core.utils",
            file=sys.stderr
        )
        sys.exit(1)
    else:
        print("✅ PASSED: Import style check", file=sys.stdout)
        sys.exit(0)


if __name__ == "__main__":
    main()
