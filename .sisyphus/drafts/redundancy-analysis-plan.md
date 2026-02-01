# Draft: RPH Redundancy Analysis Report Plan

## User Request Summary
Generate comprehensive redundancy analysis report for ReactionProfileHunter codebase after multiple iterations. NO IMPLEMENTATION - report only with file/line references, redundancy types, impact assessment, and safe removal notes.

## Key Areas to Analyze
- rph_core/orchestrator.py
- rph_core/steps/* (step2_retro, step3_opt, step4_features)
- rph_core/utils/* (qc_interface, orca_interface, xtb_runner, checkpoint_manager)
- checkpoint/resume mechanisms
- post_qc_enrichment
- SPMatrixReport

## Recent Context
- Step3 energy handling refactored
- UI updates applied
- Notify system additions
- Possible leftover code from iterations

## Open Questions (Awaiting User Response)
1. **Redundancy types priority**: Dead code, duplication, deprecated patterns, overlapping functionality, or all?
2. **Backup directory**: Include rph_core_backup_20260115/ for comparison or exclude?
3. **Safety threshold**: Conservative, moderate, or aggressive for "safe to remove" classification?
4. **Recent changes**: Are they complete or might have leftovers?
5. **Report format**: Markdown tables, categorized, prioritized, or all?

## Research Completed ✓
- Agent bg_a1abf022: Codebase structure and redundancy patterns ✓
- Agent bg_ea10dcaa: Recent changes and deprecated patterns ✓
- Agent bg_0fde8866: Duplicate code via AST analysis ✓

## Key Findings Summary

### Critical Redundancies Found:
1. **Duplicate imports** in step3_opt/__init__.py (3x TSValidator, 2x post_qc_enrichment)
2. **Duplicate IRCResult classes** (irc_driver.py + validator.py - name conflict!)
3. **Triple LogParser implementations** (geometry_tools + qc_interface + log_parser.py)
4. **Duplicate toxic path functions** (path_compat.py + qc_interface.py)
5. **Duplicate hashlib import** in result_inspector.py

### High-Impact Patterns:
6. **Executable finding duplicated 3x** (xtb_runner, orca_interface should use resource_utils)
7. **Directory creation 50x** (.mkdir(parents=True, exist_ok=True) everywhere)
8. **Subprocess execution 19x** (qst2_rescue, irc_driver bypass qc_interface)
9. **Energy parsing duplicated 7x** (SCF Done: regex copy-pasted)
10. **old_checkpoint threading** through 22 files, 5 layers deep

### Deprecated/Legacy Code:
11. **Step2 dual API** (run() vs run_with_precursor() - backward compat only)
12. **Schema migration scaffolding** (mech_index_v1, well-isolated)
13. **Empty QC artifact dicts** (V6 removed NBO, code remains)
14. **Backup directory** (rph_core_backup_20260115/ - 18,754 lines frozen copy)

### Recent Changes Status:
- **Notify system**: Clean implementation, coexists with logger (intended)
- **Step3 energy handling**: L2 energy override added, but parsing scattered
- **UI updates**: rich console.print used appropriately

## Redundancy Types to Detect
- **Dead Code**: Unused functions/classes with zero references
- **Code Duplication**: Same logic in multiple locations
- **Deprecated Patterns**: Old implementations kept alongside new
- **Overlapping Functionality**: Multiple ways to achieve same goal
- **Commented-out Code**: Large blocks of disabled code
- **Orphaned Imports**: Unused import statements
- **Test Gaps**: Production code with no test coverage (potential dead code indicator)

## Report Structure (Proposed)
```markdown
# ReactionProfileHunter Redundancy Analysis

## Executive Summary
- Total redundancies found: N
- High-impact removal candidates: M
- Safe to remove: K

## 1. Dead Code
| File:Line | Function/Class | References | Safe to Remove? | Notes |
|-----------|----------------|------------|-----------------|-------|

## 2. Code Duplication
| Original | Duplicate | Similarity | Consolidation Target | Impact |
|----------|-----------|------------|---------------------|--------|

## 3. Deprecated Patterns
| File:Line | Pattern | Replacement | Migration Needed? | Notes |
|-----------|---------|-------------|-------------------|-------|

## 4. Overlapping Functionality
| Feature | Implementations | Recommended Keep | Removal Impact |
|---------|-----------------|------------------|----------------|

## 5. Commented-Out Code
| File:Line Range | LOC | Context | Safe to Delete? |
|-----------------|-----|---------|-----------------|

## 6. Cleanup Recommendations
Prioritized list by impact...
```

## Next Steps
1. Wait for explore agent results
2. Clarify open questions with user
3. Finalize report structure
4. Generate comprehensive redundancy analysis plan
