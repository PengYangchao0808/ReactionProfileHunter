# step2_retro/AGENTS.md

## OVERVIEW
Retro scan: identifies forming bonds from product SMILES via SMARTS matching, then generates TS guess (bonds stretched to TS distances + constrained opt) and intermediate (bonds elongated to dissociation + relaxed opt). **Supports both retro_scan and forward_scan strategies.**

## WHERE TO LOOK
| File | Role |
|------|------|
| `retro_scanner.py` | Main engine; `run()` (retro) + `run_forward_scan()` (forward); handles v2.1/v3.0/v6.1 layouts |
| `smarts_matcher.py` | SMARTS pattern matching + SMARTSTemplate registry ([5+2], [4+3], [4+2], [3+2]) |
| `bond_stretcher.py` | Geometry manipulation: stretch_bonds() for arbitrary bond sets |

## KEY METHODS
- `RetroScanner.run(product_xyz, output_dir)` — legacy retro scan
- `RetroScanner.run_forward_scan(product_xyz, forming_bonds, config)` — xTB forward scan (NEW)
- `SMARTSMatcher.find_reactive_bonds(product_xyz, cleaner_data=None)` — template registry + cleaner-first

## REQUIRED OUTPUTS (both mandatory)
- `ts_guess.xyz` — initial TS geometry for S3 Berny
- `intermediate.xyz` — reaction intermediate for S3 optimization and S4 analysis (legacy alias: `reactant_complex.xyz`)
- Forming bond atom indices — passed to `utils/forming_bonds_resolver.py` for S4 fragment split

## CONVENTIONS
- Internal atom indices are **0-based**; convert to **1-based** when writing constraint files for Gaussian.
- Two input layout modes are supported inside `retro_scanner.py` (branch on presence of `S1_ConfGeneration/`); do not add a third layout without updating the orchestrator.
- `run_forward_scan` parameters driven by `reaction_profiles` config: `scan_start_distance`, `scan_end_distance`, `scan_steps`, `scan_mode`, `scan_force_constant`.
- Backward compatibility: both `intermediate.xyz` (new) and `reactant_complex.xyz` (legacy) are written

## FORMING BONDS NOTE
The forming bond indices produced here flow through `utils/forming_bonds_resolver.py` (called in orchestrator between S3 and S4) and are the **single authoritative source** for fragment splitting in S4. Do not recompute or re-derive forming bonds inside S4 extractors.

## ANTI-PATTERNS
- Producing only `ts_guess.xyz` without `intermediate.xyz` — both are required for S3/S4.
- Implementing batch scheduling inside this module — orchestrator/batch layer owns that.
