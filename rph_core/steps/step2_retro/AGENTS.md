# step2_retro/AGENTS.md

## OVERVIEW

V4 S2 PEB + xTB PATH pipeline. A coarse xTB scan (0.20 A) locates the
topology drift boundary, then xTB2 meta-dynamics PATH runs from the
last-valid frame before drift to the product. Every PATH frame is
refined at ORCA B97-3c SP, and TS/INT are selected from the refined
profile (PATH TS estimate first, Kneedle fallback).

Only the canonical PEB+PATH flow is supported. Its sole escalation is the
single-ended ORCA B97-3c relaxed-scan rescue, and it is used only when
topology drift persists from the first PATH frame so no usable path segment
exists. The historical NEB-assisted flow and corridor refinement remain
removed.

## WHERE TO LOOK

| File | Role |
|---|---|
| `peb_scanner.py` | Canonical V4 entrypoint |
| `peb_engine.py` | Coarse scan, anchor detection, xTB PATH, B97-3c SP, node selection |
| `scan_trajectory.py` | Immutable `ScanAttempt` records (composite stitching no longer used) |
| `energy_refinement.py` | Geometry-hash-cached B97-3c SP refinement on PATH frames |
| `geometry_guard.py` | Topology and risky-contact diagnostics for the coarse scan |
| `scan_policies.py` | Concerted/asymmetric xTB scan policies |
| `bond_stretcher.py` | Standalone utility (no longer wired into the engine) |

## PIPELINE

```
coarse xTB scan (0.20 A, policy_c)
   |
   +-- topology guard + anchor detection (P, M, D, E)
   |
   v
xTB2 meta-dynamics PATH
   start = last_valid_before_drift frame
   end   = product frame
   npoint = 28 (configurable)
   |
   +-- per-frame xTB SP to populate PATH energies
   |
   v
B97-3c SP on ALL PATH frames (full coverage, hash-cached)
   |
   v
TS/INT seed selection from refined profile
   TS  = B97-3c local maximum on the PATH energy curve
   INT = pre-TS B97-3c basin | pre-TS arc-length midpoint search seed
   |
   v
plot_scan_profile (xTB coarse overlay + B97-3c PATH curve)
   |
   v
full-path-distortion gate -> S3 TS/INT seeds | B97-3c relaxed-scan rescue
   |
   v
qualified ts_guess.xyz / intermediate.xyz, scan_profile.json -> S3
```

## NODE SELECTION CONTRACT

- **TS** is the highest topology-valid internal local maximum on the
  method-consistent refined PATH curve. xTB PATH estimates are diagnostics only.
- **INT search seed** is a significant local minimum between the reactant-side
  endpoint and TS. If no basin is resolved, use the nearest existing PATH frame
  to their arc-length midpoint. It never claims to be a stationary point.
- B97-3c coverage is evaluated across the entire PATH. A partial B97-3c
  failure cascades to both TS and INT selections and triggers xTB fallback.

## OUTPUTS

- `ts_guess.xyz` - TS guess frame for S3 OptTS
- `intermediate.xyz` - INT search seed for S3 MIN optimization; it may share
  the TS frame and never claims a stationary point
- `xtb_path/` - xTB PATH run directory with `xtbpath_NNN.xyz` and per-frame SP logs
- `energy_refinement/` - B97-3c SP cache (SHA256-keyed)
- `scan_profile.json` - schema `s2_scan_profile_v10`; contains coarse scan summary,
  PATH metadata (barriers, estimated TS), all PATH frames with xTB and B97-3c
  energies, anchor payload, and selection records
- `scan_profile.png` - the sole root-level S2 energy-profile figure
- `rescue/scan_profile.json` and `rescue/scan_profile.png` - B97-3c
  relaxed-scan profile and figure, written only for a full-path-distortion
  rescue and rendered through the same S2 plotting entrypoint

## CONVENTIONS

- Indices in manifests are 0-based.
- Forming bonds are unordered edges via `canonicalize_bond_pairs`.
- Every external xTB/ORCA job goes through the existing QC interfaces.
- All scientific parameters live in `config/defaults.yaml` under `step2.*`.
- The coarse scan still uses the PEB scanner (topology guard, policy_c).
- PATH frames use mean forming-bond distance as the reaction coordinate so
  the PATH curve overlays the coarse PEB curve on the same x-axis.
- B97-3c SP is run via `ScanEnergyRefiner` with SHA256 geometry caching.

## ANTI-PATTERNS

- Reintroducing the corridor refinement or independent TS/INT windows.
- Reintroducing NEB-assisted flow (`S2NEBCoordinator`, `IntermediateBasinFinder`).
- Treating a TS/INT seed as proof of a stationary point, or using S2 seed
  labels to decide an S3 endpoint/NEB workflow.
- Running B97-3c SP only on a subset of PATH frames (must be full coverage).
- Treating the absolute highest-energy point as the TS guess.
- Direct QC subprocess calls or hard-coded scientific parameters.
- Re-enabling the `refinement_corridor` config section (removed; no effect).
