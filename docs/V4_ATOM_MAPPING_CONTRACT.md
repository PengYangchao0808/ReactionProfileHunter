# V4 atom-mapping contract

V4 uses one executable atom order from S0 through S4: the RDKit atom order of
the exact mapped branch SMILES passed to S1. Canonical/unmapped SMILES may be
used for identity and display, but their local atom indices must not be used by
QC constraints.

## Stage artifacts

- S0 writes one `S0_Mechanism/atom_mappings/<variant>.json` per product branch.
  It records map number, branch-SMILES index, forming bonds in both spaces, the
  exact reference SMILES and a mapping digest. `variant_registry.json` points
  each branch to its own table.
- S1 writes `initial_atom_mapping.json` when RDKit creates the geometry and
  `atom_mapping.json` for `selected.xyz`. The table includes every heavy atom
  and hydrogen. Atom count, element sequence and XYZ hash are verified.
- S2 resolves forming bonds only by composing the S0 and S1 tables. It writes
  `atom_mapping.json` for the PEB map/SMILES/XYZ audit and complete atom tables
  for the product, intermediate and TS geometries.
- S3 and S4 verify their input atom table before launching QC, bind a new table
  to the optimization/SP geometry, and propagate the mapping fields in their
  manifests. A mapping-required structure cannot be `usable_for_ml` unless its
  mapping status is `verified`.

## Index bases

- Atom-map numbers are the positive 1-based labels embedded in mapped SMILES.
- RDKit SMILES and internal XYZ indices are 0-based.
- xTB/ORCA input conversion remains the responsibility of the existing QC
  interfaces; stage manifests and `forming_bonds` remain 0-based.

Missing sidecars, hash mismatches, atom-count changes, element-order changes,
unresolved map numbers, duplicate indices, additive atoms, or a scanner that
changes verified forming bonds are fail-fast mapping errors. S0/S1/S2 schema
and signature versions were advanced so old V4 checkpoints cannot silently
reuse artifacts produced before this contract.

## Artifact-first checkpoint reconciliation

`pipeline.state` is a rebuildable cache index, not the authority for whether a
QC result exists. On S1 resume, V4 validates the existing manifest,
`selected.xyz`, original generation SMILES, complete element order and current
atom mapping. A scientifically unchanged legacy S1 result is migrated without
running QC: the old manifest is preserved under `migrations/`, a current atom
mapping sidecar is generated, the current manifest is written atomically, and
the checkpoint is materialized from that validated evidence.

Scientific signature changes are not migrated. Explicit `--recompute-from`
also takes precedence over reconciliation. This permits S0 and invalid S2/S3
nodes to be recalculated while retaining expensive, provably compatible S1
results.

Selective refresh is available when particular stages must be rerun without
declaring every downstream artifact invalid:

```bash
bin/rph_run ... --refresh-stages s0,s2,s3 --stop-after s3
```

Each named stage is rerun. Unnamed stages still pass through artifact
validation and checkpoint reconciliation, so a compatible legacy S1 is
migrated and reused without QC. `--refresh-stages` is intentionally mutually
exclusive with `--recompute-from` and `--start-from`.

## S2 path and node contract

Forming bonds are unordered graph edges. Every stage canonicalizes them as
sorted `(min(i, j), max(i, j))` pairs before comparison, hashing or manifest
serialization. Reversing the two atoms can therefore never invalidate an S2
result by itself.

S2 publishes explicit nodes rather than assuming that every path contains an
intermediate. `ts_guess.xyz` is selected from the topology-valid energy peak.
`intermediate.xyz` exists only when at least three contiguous topology-valid
post-TS frames define a local minimum or sufficiently flat basin. When no basin
exists, the manifest records `intermediate.status=not_found`; S3 receives the
product and TS but no invented minimum.

The mechanically stretched geometry is diagnostic only and is written as
`diagnostics/stretched_seed.xyz`. S2 no longer writes the duplicate legacy
`reactant_complex.xyz`. Before a refreshed scan, previously published nodes
and manifests are moved under `superseded/<run-id>/`, preventing a failed run
from exposing stale structures as current results.

The default S2 trajectory uses a maximum coarse spacing of 0.08 Angstrom and a
0.03 Angstrom local refinement around the xTB peak. If the peak lacks enough
dissociation-side coverage, the upper endpoint is extended in bounded 0.25
Angstrom increments until the peak and a topology-valid post-TS segment are
bracketed, persistent topology drift stops extension, or the configured hard
limit is reached. Persistent drift can also trigger a keep-away constraint and
policy retry. Fixed keep-away distances are distinct from scan coordinates in
the xTB input.

Every coarse, endpoint-extension, refinement and topology-retry calculation is
retained as an immutable scan attempt. A full-range topology-valid attempt is
the composite backbone; accepted fine attempts replace only their overlapping
region after coordinate, atom-order, aligned-RMSD and xTB-energy continuity
checks. The published reaction path is this multi-resolution composite, never
merely the last narrow scan.

When `step2.energy_refinement.enabled` is true, every topology-valid point in
the composite receives a cached ORCA B97-3c/CPCM single-point calculation
through `qc_jobs.py`. Cache identity is the geometry hash plus QC specification,
so frame renumbering does not trigger duplicate work. A bracketed B97-3c curve
selects the TS and post-TS basin. Missing/partial DFT coverage or method
disagreement can use an explicitly labelled low-confidence xTB fallback, and a
B97-3c boundary point is never silently promoted as a verified TS. xTB mRRHO is
deliberately not applied to arbitrary non-stationary scan frames.

`s2_scan_profile_v4` stores explicit `attempts`, the assembled
`composite_profile`, independent `energy_curves.xtb` and
`energy_curves.b973c`, and the actual `selection_policy`. Plotting always shows
the complete composite xTB profile, overlays available B97-3c points, and marks
refinement coverage and topology drift. The ambiguous legacy
`energies_hartree` field remains only as a compatibility view and is not the
authoritative curve contract.
