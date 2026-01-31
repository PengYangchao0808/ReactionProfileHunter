# Decisions

- Implemented molecule_utils.py in rph_core/utils/ to centralize SMILES canonicalization and formula calculation.
- Used rdMolDescriptors.CalcMolFormula for stable molecular formula generation.
- Defined small molecules based on heavy atom count threshold (default 10) to support global small molecule directory caching.
- Used get_molecule_key as the primary directory naming strategy to ensure consistency with other parts of the pipeline that might search for these artifacts.

## Small Molecule Cache in S1
- **Decision**: Skip `ConformerEngine` entirely upon cache hit in `AnchorPhase`.
- **Rationale**: `ConformerEngine` instantiation and call overhead, though small, can be avoided. Populating the local directory manually from cache provides full control and ensures downstream steps (which expect files in a specific local path) still work.
- **Decision**: Extract energy from XYZ comment in `AnchorPhase` for cached molecules.
- **Rationale**: Reuses the "Gold Standard" extraction logic from `ConformerEngine` without needing to run the engine itself.

### Directory Atomicity
- Decided to maintain strict separation of molecule data in Step 1.
- Each molecule's calculations are isolated to avoid concurrency conflicts (if parallelization is added later) and to improve organization.
- Used role-based keys ("product", "precursor", "leaving_group") as default directory names in the orchestrator.
- Updated RetroScanner path resolution to support v6.1 flat structure (product_min.xyz directly in S1_ConfGeneration) while maintaining v3.0 subdirectory support and v2.1 file support.
