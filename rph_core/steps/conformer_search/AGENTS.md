# rph_core/steps/conformer_search/AGENTS.md

## OVERVIEW
S1 Unified Conformer Engine (UCE) v3.1: two-stage xTB conformer search (GFN0 coarse → GFN2 fine) followed by DFT OPT/SP coupling. Backward-compatible single-stage mode also supported.

## TWO-STAGE WORKFLOW (default, v3.1)
1. **Stage 1**: GFN0-xTB rapid sampling → ISOSTAT clustering → `stage1_gfn0/cluster/cluster.xyz`
2. **Stage 2**: GFN2-xTB fine optimization of stage-1 cluster → ISOSTAT clustering → `stage2_gfn2/cluster/cluster.xyz`
3. **DFT**: OPT+SP on stage-2 ensemble → `finalDFT/` (Gaussian or ORCA)

**Single-stage mode** (set `two_stage_enabled: false`): direct GFN2-xTB → cluster → DFT.

## WHERE TO LOOK
- Engine: `engine.py` → `ConformerEngine`
- QC facade: `rph_core/utils/qc_interface.py` (xTB/CREST/Gaussian factory)
- OPT/SP loop: `rph_core/utils/qc_task_runner.py`
- Geometry/log parsing: `rph_core/utils/geometry_tools.py`

## CONFIG SURFACES
| Key | Type | Default | Notes |
|-----|------|---------|-------|
| `step1.conformer_search.two_stage_enabled` | bool | `true` | Master switch |
| `step1.conformer_search.stage1_gfn0.gfn_level` | int | `0` | 0=GFN0 |
| `step1.conformer_search.stage1_gfn0.energy_window_kcal` | float | `10.0` | Stage-1 energy window |
| `step1.conformer_search.stage2_gfn2.gfn_level` | int | `2` | 2=GFN2 |
| `step1.conformer_search.stage2_gfn2.energy_window_kcal` | float | `3.0` | Stage-2 energy window |
| `step1.crest` | dict | — | Single-stage fallback config |
| `executables.isostat`, `.shermo`, `.gaussian.wrapper_path` | str | — | Binary paths |

## DIRECTORY LAYOUT

**Two-stage mode:**
```
S1_ConfGeneration/<molecule_name>/
├── crest/
│   └── ensemble.xyz
├── xtb/
│   ├── stage1_gfn0/
│   │   ├── crest_conformers.xyz
│   │   └── cluster/
│   │       ├── cluster.xyz        # input to Stage 2
│   │       └── isostat.log
│   ├── stage2_gfn2/
│   │   ├── crest_ensemble.xyz
│   │   └── cluster/
│   │       └── cluster.xyz        # final ensemble → DFT
├── cluster/                       # single-stage clustering output
├── prescan/                       # full-protocol fast-SP pre-screening
├── fastsp/                        # full/lite fast-SP screening
└── finalDFT/                      # DFT OPT/SP outputs
```

**Single-stage mode:** `crest/crest_conformers.xyz` → `cluster/cluster.xyz` → `finalDFT/`

## GOTCHAS
- `engine.py` uses `subprocess.run(..., shell=True, cwd=...)` — path/escaping sensitive; prefer `utils` sandbox/toxic-path helpers when touching those calls.
- GFN0 is ~10× faster than GFN2 but lower accuracy; two-stage strategy covers conformer space cheaply then refines.
- Two-stage default; single-stage only for backward compatibility with older configs.

## ANTI-PATTERNS
- Hardcoding `g16`/`orca`/`xtb` binary paths inside engine — use `config['executables']` or utils lookup.
- Using `rph_core_backup_20260115/` behavior as a reference for new logic.
