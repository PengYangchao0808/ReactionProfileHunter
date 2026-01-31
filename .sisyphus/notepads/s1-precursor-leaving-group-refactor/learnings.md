# Learnings

- RDKit is the primary library for molecular processing.
- Optional[T] is preferred over T | None for compatibility with Python 3.9 as per AGENTS.md.
- Chem.MolFromSmiles returns None on failure, requiring explicit checks despite some LSP warnings.
- SmallMoleculeCache provides a clean interface for managing global molecule artifacts using SMILES-derived keys.

## S1_Product References Scan Results
Total matches found: 48 (excluding binary and output dirs)

### 1. Path Strings (Hardcoded literals)
- **Core Logic**:
    - `rph_core/orchestrator.py`: `"S1_Product"` used in `s1_work_dir` assignment and path joins.
    - `rph_core/steps/step4_features/mech_packager.py`: `S1_DIR_ALIASES = ["S1_Product"]`
    - `rph_core/utils/result_inspector.py`: Multiple joins with `"S1_Product"`.
- **Tests**:
    - Extensive use in `tests/` (e.g., `test_m2_precursor_fallback.py`, `test_s4_mech_packager.py`) usually as `tmp_path / "S1_Product"`.

### 2. Variable Names & Aliases
- `s1_work_dir` (orchestrator.py)
- `S1_DIR_ALIASES` (mech_packager.py)
- `s1_dir` (Commonly used in function signatures and tests)

### 3. Comments & Documentation
- **Markdown**: `README.md`, `rph_core/steps/anchor/AGENTS.md`, `rph_core/steps/conformer_search/AGENTS.md`.
- **Docstrings**:
    - `rph_core/orchestrator.py`
    - `rph_core/steps/anchor/handler.py`
    - `rph_core/steps/conformer_search/engine.py`
    - `rph_core/steps/step2_retro/retro_scanner.py`
    - `rph_core/steps/step4_features/mech_packager.py`

### Observations
- Most hardcoded strings are for directory creation or artifact lookup.
- `S1_DIR_ALIASES` in `mech_packager.py` suggests a mechanism for directory name flexibility already exists in Step 4.

## AnchorPhase Small Molecule Cache Integration
- Integrated `SmallMoleculeCache` into `AnchorPhase` to skip redundant conformer searches for small molecules.
- Cache root defaults to `base_work_dir.parent / "SmallMolecules"`.
- Cache hit detection logic:
  - If `is_small_molecule(smiles)` and `small_mol_cache.exists(smiles)`:
    - Populates local `[molecule_name]/dft/` directory from cache.
    - Extracts SP energy from the cached XYZ comment.
    - Skips `ConformerEngine.run()`.
- Cache save logic:
  - After successful `ConformerEngine.run()`, if molecule is small and was not a cache hit, saves results to global cache.
- Verified with integration tests mocking the QC engine.

## Orchestrator Multi-Molecule Integration (2026-01-31)
- Modified `rph_core/orchestrator.py` to support precursor and leaving group molecules in S1.
- Initialized `SmallMoleculeCatalog` in `ReactionProfileHunter.__init__`.
- `run_pipeline` now accepts `precursor_smiles` and `leaving_group_key`.
- `leaving_group_key` is resolved to SMILES using the catalog; if resolution fails, a warning is logged and the molecule is skipped.
- `_run_tasks` extracts these fields from `task.meta`.
- Verified with unit tests in `tests/test_orchestrator_multi_molecule.py`.

### Multi-Molecule Directory Structure (v3.0)
- AnchorPhase correctly iterates through all molecules (product, precursor, leaving_group).
- Each molecule gets its own subdirectory: `S1_Product/[Molecule_Name]/`.
- ConformerEngine manages internal subdirs (`xtb2/`, `cluster/`, `dft/`) under the molecule directory.
- Global minimum files are named `[Molecule_Name]_global_min.xyz`.
- All paths are resolved to absolute paths to prevent I/O issues in QC programs.

## S1_Product to S1_ConfGeneration Rename
- Successfully renamed `S1_Product` to `S1_ConfGeneration` across the codebase.
- Maintained backward compatibility:
    - `rph_core/orchestrator.py`: `s1_candidates` list now includes both old and new names when skipping Step 1.
    - `rph_core/steps/step4_features/mech_packager.py`: `S1_DIR_ALIASES` now includes `["S1_ConfGeneration", "S1_Product"]`.
    - `rph_core/utils/result_inspector.py`: All search candidates now include both old and new paths.
- Updated all test files in `tests/` using `sed`.
- Updated all `README.md` and `AGENTS.md` files (except in backup directory).
- Updated Docstrings and comments in source files.
- Robust path resolution in S2: It now checks for product_min.xyz in both the root of S1_ConfGeneration and in molecule subdirectories.
