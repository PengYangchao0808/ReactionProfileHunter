# rph_core/steps/anchor/AGENTS.md

## OVERVIEW
S1 anchor phase: conformer search + DFT OPT/SP for each molecule in the substrate pool (product, precursor, leaving group, small molecules). Outputs the global minimum-energy structure and SP energy per molecule.

## WHERE TO LOOK
- Coordinator: `handler.py` → `AnchorPhase.run()`
- Conformer engine: `rph_core/steps/conformer_search/engine.py` → `ConformerEngine`
- QC facade: `rph_core/utils/qc_interface.py` → `try_formchk()` for .chk → .fchk conversion
- Small molecule cache: `rph_core/utils/small_molecule_cache.py` → avoids re-running for common small molecules

## OUTPUTS
Written to `S1_ConfGeneration/<molecule_name>/` (managed by `ConformerEngine`):
```
<molecule_name>/
├── xtb2/          # xTB conformer search (stage1_gfn0/, stage2_gfn2/ in two-stage mode)
├── cluster/       # ISOSTAT clustering output
└── dft/           # Gaussian/ORCA OPT+SP (conformer_thermo.csv, *.fchk, *.log)
```

Returned in `anchor_result.anchored_molecules[name]`:
- `xyz` (Path): final minimum geometry
- `e_sp` (float): single-point energy (Hartree)
- `log`, `fchk`, `qm_output`, `checkpoint` (optional Paths)
- `product_thermo` (optional Path): `dft/conformer_thermo.csv`

## CONVENTIONS
- All paths via `pathlib.Path`; never hand-concatenate string paths here.
- `.chk` → `.fchk` conversion uses `qc_interface.try_formchk()` — never call `formchk` directly.
- `SmallMoleculeCache` (utils) should be checked before running conformer search on common fragments.

## ANTI-PATTERNS
- Reimplementing QC runner / sandbox inside S1 — use `utils/qc_interface.py`.
- Renaming output files without syncing orchestrator + downstream step contracts (S2/S3/S4 path lookups will break).
