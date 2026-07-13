# Lewis Acid Module — Design Document (v2)

**Date**: 2026-06-11  
**Status**: Final — Ready for Implementation  
**Target**: RPH v3.0.0+

---

## 0. Design Decision

> **LiCl is adopted as the default neutral lithium-salt chelation surrogate.**

This is **not** a "universal Lewis acid replacement." LiCl is a minimal probe that generates Lewis-acid-sensitive conformational, coordination, and electronic-response features at negligible computational cost. Real Lewis acid identity enters only through an independent descriptor registry consumed in ML/post-analysis.

---

## 1. Chemical Rationale

### 1.1 Why not naked Li⁺

Li⁺'s first coordination shell is dramatically altered by its counteranion. In THF, LiCl does not exist as free Li⁺/Cl⁻ — it forms a dynamic ensemble of (LiX)ₙ aggregates, THF-coordinated species, and contact ion pairs (2024 *Chemical Science* LiCl/THF study). Naked Li⁺ would overestimate coordination strength and miss anion-induced steric/preorganization effects.

### 1.2 Why LiCl is preferred

| Criterion | LiCl | Naked Li⁺ | BF₃ | BH₃ |
|-----------|------|-----------|-----|-----|
| Cost (atoms beyond organic) | 2 | 1 | 4 | 4 |
| Neutral system? | Yes (charge 0) | No (+1) | Yes | Yes |
| Captures ion-pair dynamics? | Partially | No | N/A | N/A |
| CREST/xTB stable? | Yes (neutral) | Requires counterion | Yes | Yes |
| Known aggregation behavior | Yes (literature) | Unphysical | Well-characterized | Well-characterized |

LiCl is chosen as the **pragmatic minimum**: neutral, 2 atoms, cheap, and chemically grounded in lithium-salt solution chemistry.

### 1.3 Why BF₃/BH₃ are not default

BF₃ and BH₃ are molecular Lewis acids, not lithium salts. They would require different precomputed descriptors and cannot serve as general probes for salt-based Lewis acid conditions (LiClO₄/Et₂O, ZnCl₂, etc.). They remain valid as **additional surrogate modes** if a specific reaction calls for them.

### 1.4 LiClO₄ explicit counteranion mode

For reactions where the experimental condition is explicitly LiClO₄ in ether (Grieco 1990: Diels–Alder acceleration in 5 M LiClO₄/Et₂O), a dedicated `special_salts.LiClO4` mode is available. Default remains LiCl due to lower conformational complexity (ClO₄⁻ adds 5 atoms + 4 O coordination sites).

---

## 2. Architecture: Dual-Part Design

```
┌─────────────────────────────────────────────────────────────────────┐
│                    LEWIS ACID MODULE — v2                           │
│                                                                     │
│  ┌─────────────────────────────┐    ┌────────────────────────────┐ │
│  │  Part A: Descriptor Registry│    │  Part B: LiCl Additive     │ │
│  │  (independent precompute)   │    │  (pipeline integration)    │ │
│  │                             │    │                            │ │
│  │  Per LA in library:         │    │  Config-driven:            │ │
│  │  ├─ acetone_affinity        │    │  ├─ geometry-only LiCl     │ │
│  │  ├─ thf_affinity (opt.)     │    │  ├─ additive atom tracking │ │
│  │  ├─ lumo_energy             │    │  ├─ CREST NCI sampling     │ │
│  │  ├─ coordination_number     │    │  ├─ S2 quality gate        │ │
│  │  └─ counterion_type         │    │  ├─ S3 imaginary mode gate │ │
│  │                             │    │  └─ S4 LiCl feature mining │ │
│  │  NOT a ΔG‡ corrector.       │    │                             │ │
│  │  Used as ML features only.  │    │  Canonical SMILES NEVER     │ │
│  │                             │    │  mutated. Li/Cl always      │ │
│  └─────────────────────────────┘    │  tracked as additive atoms. │
│                                      └────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────┘
```

---

## 3. Part A — Lewis Acid Descriptor Registry

### 3.1 Purpose

Not a ΔG‡ corrector. Lewis acid descriptors are **independent ML/post-analysis features** that interact with LiCl-derived geometry/electronic-response features through learned interaction terms.

```python
# ML target ~ reaction outcome
# Features:
#   + organic geometry features (existing S4)
#   + LiCl chelation-response features (new, Section 7)
#   + real LA acetone_affinity (Part A)
#   + real LA specific affinities (Part A, optional)
#   + interaction terms (learned)
```

### 3.2 v1 Descriptor Set

```json
{
    "lewis_acid": "AlCl3",
    "descriptors": {
        "acetone_affinity_kcal_mol": -12.34,
        "thf_affinity_kcal_mol": null,
        "acetonitrile_affinity_kcal_mol": null,
        "lumo_energy_ev": -1.23,
        "preferred_coordination_number": 4,
        "formal_charge": 0,
        "counterion_type": "none"
    }
}
```

v1 only computes **acetone affinity** as a first-order carbonyl Lewis acidity descriptor. Additional donor affinities (THF, MeCN) can be added later.

### 3.3 Computation Protocol

```
For each LA in library:
  1. DFT opt of acetone        → E_acetone
  2. DFT opt of free LA        → E_LA  
  3. DFT opt of acetone·LA     → E_complex
  4. ΔE = E_complex − E_acetone − E_LA  (Hartree → kcal/mol)
```

Uses the same theory level as pipeline optimization (`theory.optimization`).

### 3.4 Storage

`lewis_acid_registry.json` (cached globally):

```json
{
    "schema_version": "la_registry_v1",
    "theory": {"method": "B3LYP", "basis": "def2-SVP", ...},
    "acids": {
        "LiCl":  {"acetone_affinity_kcal_mol": -0.06, ...},
        "AlCl3": {"acetone_affinity_kcal_mol": -12.34, ...},
        "BF3":   {"acetone_affinity_kcal_mol": -15.67, ...},
        "ZnCl2": {"acetone_affinity_kcal_mol": -8.90, ...}
    }
}
```

### 3.5 No Linear ΔG‡ Correction

**Removed.** The v1 design contained:

> ΔG‡(LA) ≈ ΔG‡(LiCl) + α × [ΔE_desc(LA) − ΔE_desc(LiCl)]

This assumption is too strong. Real Lewis acid identity affects both electronic structure and conformational ensembles in ways a single scalar cannot capture. Registry descriptors enter ML as independent features, not hard-coded corrections.

---

## 4. Core Data Model: Additive Atoms

### 4.1 The Central Rule

| What | Value |
|------|-------|
| `canonical_product_smiles` | Original organic product SMILES (NEVER mutated) |
| `geometry_product_smiles` | Original SMILES + `.[Li]Cl` (for S1 geometry only) |
| `checkpoint_identity_hash` | `hash(canonical SMILES + LA config + theory)` — NOT the geometry SMILES |
| `SMARTS matching input` | Built from original product SMILES, NOT from XYZ |
| `XYZ / DFT geometry` | Organic product + LiCl atoms |
| `additive_atom_indices` | Explicitly tracked Li/Cl indices |

### 4.2 LewisAcidAdditive Dataclass

```python
@dataclass
class LewisAcidAdditive:
    """Tracks LiCl additive atoms. NEVER mutates canonical SMILES."""
    enabled: bool = False

    # Surrogate identity
    mode: str = "neutral_lithium_salt_chelation_surrogate"
    surrogate_name: str = "LiCl"
    surrogate_smiles: str = "[Li]Cl"

    # Charge / multiplicity (overall system)
    charge: int = 0
    multiplicity: int = 1

    # Tracking
    additive_atom_indices: Tuple[int, ...] = ()      # 0-based in XYZ
    additive_elements: Tuple[str, ...] = ("Li", "Cl")

    # Coordination info (filled post-S1)
    li_coordination_site: Optional[str] = None
    li_o_distance_initial: Optional[float] = None
    li_cl_distance_initial: Optional[float] = None

    # Config snapshot for checkpoint hashing
    config_snapshot: Optional[Dict[str, Any]] = None
```

### 4.3 PipelineResult Extension

```python
@dataclass
class PipelineResult:
    # ... existing fields ...
    lewis_acid_additive: Optional[LewisAcidAdditive] = None
```

No `lewis_acid_enabled`, no `lewis_acid_surrogate`, no `lewis_acid_name`, no `original_product_smiles` — these are redundant information already in `LewisAcidAdditive`.

---

## 5. Checkpoint Identity (Rewritten)

### 5.1 What NOT to do

**Do NOT store `product.[Li]Cl` as the canonical product identity.**

The checkpoint identity must distinguish LA-on from LA-off without corrupting the product SMILES:

```python
checkpoint_identity = hash({
    "product_smiles": product_smiles,          # Original, clean
    "precursor_smiles": precursor_smiles,      # Original, clean
    "reaction_type": reaction_type,
    "lewis_acid": {
        "enabled": True,
        "mode": "neutral_lithium_salt_chelation_surrogate",
        "surrogate_name": "LiCl",
        "surrogate_smiles": "[Li]Cl",
    },
    "theory": theory_config,
    "solvent": solvent_config,
})
```

This hash is appended to the checkpoint metadata. A previous run without LA will produce a different hash → correctly treated as a different computation.

### 5.2 Pipeline State

```python
# In pipeline.state (CheckpointManager):
{
    "product_smiles": "CC1CC(=O)C(=O)O",   # Clean
    "precursor_smiles": "C=C(C)C(=O)OC",   # Clean
    "lewis_acid_additive": {
        "enabled": True,
        "mode": "neutral_lithium_salt_chelation_surrogate",
        "surrogate_name": "LiCl",
        "surrogate_smiles": "[Li]Cl",
    },
    "checkpoint_identity_hash": "a1b2c3d4..."
}
```

---

## 6. Pipeline Integration (Part B)

### 6.1 Config

```yaml
# ==================== Lewis Acid Module ====================
lewis_acid:
  enabled: false
  mode: "neutral_lithium_salt_chelation_surrogate"

  default_surrogate:
    name: "LiCl"
    smiles: "[Li]Cl"
    charge: 0
    multiplicity: 1

  special_salts:
    LiClO4:
      enabled: false
      smiles: "[Li+].[O-]Cl(=O)(=O)=O"
      charge: 0
      multiplicity: 1
      use_only_when_experimental_salt_is: "LiClO4"

  # S1 sampling controls
  s1_sampling:
    use_nci_mode: true
    use_wall_potential: true
    weak_li_o_restraint: optional
    max_li_cl_distance: 3.5        # Å
    max_li_substrate_distance: 4.0 # Å

  # S2 quality gate thresholds
  s2_quality_gate:
    li_cl_max: 3.5                 # Å
    li_o_max: 2.6                  # Å
    cl_reaction_center_min: 3.2    # Å
    coordination_switch_allowed: false

  # S3 quality gate thresholds
  s3_quality_gate:
    reaction_mode_projection_min: 0.35
    additive_mode_projection_max: 0.25
    require_irc_for_flagged_ts: true

  # Append to precursor too?
  append_to_precursor: true
```

### 6.2 Orchestrator — Geometry Preparation

```python
def _prepare_lewis_acid_system(
    self,
    product_smiles: str,
    precursor_smiles: Optional[str] = None,
) -> Tuple[str, Optional[str], LewisAcidAdditive]:
    """
    Prepare combined SMILES for S1 geometry construction.
    Canonical product_smiles is NEVER mutated.
    
    Returns:
        (geometry_product_smiles, geometry_precursor_smiles, additive_info)
    """
    la_config = self.config.get("lewis_acid", {})
    if not la_config.get("enabled", False):
        return product_smiles, precursor_smiles, LewisAcidAdditive()

    surrogate = la_config.get("default_surrogate", {})
    la_smiles = surrogate.get("smiles", "[Li]Cl")
    la_name = surrogate.get("name", "LiCl")
    
    geometry_product = f"{product_smiles}.{la_smiles}"
    geometry_precursor = None
    if precursor_smiles and la_config.get("append_to_precursor", True):
        geometry_precursor = f"{precursor_smiles}.{la_smiles}"

    additive = LewisAcidAdditive(
        enabled=True,
        surrogate_name=la_name,
        surrogate_smiles=la_smiles,
        charge=surrogate.get("charge", 0),
        multiplicity=surrogate.get("multiplicity", 1),
        config_snapshot=la_config,
    )

    self.logger.info(f"[LA] Enabled: {la_name} ({la_smiles})")
    self.logger.info(f"[LA] Canonical SMILES unchanged: {product_smiles}")
    self.logger.info(f"[LA] Geometry SMILES: {geometry_product}")

    return geometry_product, geometry_precursor, additive
```

### 6.3 run_pipeline() Insertion

```python
# ≈ line 1885: after PipelineResult creation, before S0
la_product, la_precursor, la_additive = self._prepare_lewis_acid_system(
    product_smiles=product_smiles,
    precursor_smiles=precursor_smiles,
)
result.lewis_acid_additive = la_additive if la_additive.enabled else None

# ≈ line 2211: molecules dict uses geometry SMILES (for S1 only)
molecules = {"product": la_product}
if la_precursor:
    molecules["precursor"] = la_precursor

# ≈ line 2364: S2 receives additive info for forming bond filtering
step2_artifacts = run_step2(
    ...,
    lewis_acid_additive=la_additive,
)

# ≈ line 2426: S3 receives additive info for mode validation
step3_artifacts = run_step3(
    ...,
    lewis_acid_additive=la_additive,
)
```

### 6.4 Checkpoint Register

```python
# When saving checkpoint, include LA config in identity hash:
checkpoint_mgr.mark_step_completed("s1", output_files={...}, metadata={
    "lewis_acid_additive": asdict(la_additive) if la_additive.enabled else None,
})
```

---

## 7. S1 — NCI-Aware Conformer Generation

### 7.1 CREST Cannot Be Naïve

**Do NOT assume "CREST naturally finds the complex."** Standard CREST metadynamics can dissociate non-covalent fragments. For LiCl-enabled systems:

- **Use CREST NCI mode** (non-covalent interactions): NCI-MTD with ellipsoidal wall potential prevents fragment dissociation
- **Use wall potential** to constrain LiCl within a reasonable coordination envelope
- **Weak Li–O restraint** (optional, config-controlled): if LiCl repeatedly dissociates

### 7.2 Implementation

```python
# In ConformerEngine, when LA mode is active:
crest_kwargs = {}
if la_additive.enabled:
    crest_kwargs["nci"] = True
    crest_kwargs["wall"] = la_config["s1_sampling"]["max_li_substrate_distance"]
    if la_config["s1_sampling"]["weak_li_o_restraint"] == "optional":
        crest_kwargs["constraints"] = {"distance": {"Li": "O_receptor", "max": 4.0}}
```

### 7.3 Post-S1 Validation

After conformer search + DFT relaxation:
1. Verify Li–O contact exists (< 2.6 Å to nearest substrate O)
2. Verify Li–Cl distance is stable (~2.0–2.4 Å)
3. If LiCl dissociated → fall back to no-LA structure + mark degraded

```python
def _validate_la_geometry(xyz_path: Path, la_additive: LewisAcidAdditive) -> bool:
    atoms, coords = read_xyz(xyz_path)
    li_idx = la_additive.additive_atom_indices[0]
    cl_idx = la_additive.additive_atom_indices[1]
    
    d_licl = distance(coords[li_idx], coords[cl_idx])
    d_li_o = min(distance(coords[li_idx], coords[o_idx]) 
                 for o_idx in get_oxygen_indices(atoms, exclude=cl_idx))
    
    if d_licl > 3.5 or d_li_o > 4.0:
        return False  # Dissociated
    return True
```

---

## 8. S2 — Retro Scan with Additive Monitoring

### 8.1 Forming Bond Resolution (Critical)

**Forming bonds must be resolved from the CANONICAL product, not the geometry.**

```python
def _resolve_forming_bonds_for_s2(self, ..., lewis_acid_additive=None):
    if lewis_acid_additive and lewis_acid_additive.enabled:
        # Force SMARTS fallback on canonical product SMILES
        # Do NOT use cleaner-data bonds (indices don't account for LiCl)
        # Do NOT use geometry-difference (will pick up Li-Cl)
        return self._resolve_via_smarts_fallback(
            product_smiles=product_smiles,  # Canonical, NOT geometry
            ...,
        )
    
    # Normal resolution chain (unchanged)
```

**Hard filter applied to ALL resolution paths:**

```python
def _filter_additive_atoms(forming_bonds, additive_indices):
    """Remove any forming bond candidate involving additive atoms."""
    additive_set = set(additive_indices)
    return [(i, j) for (i, j) in forming_bonds
            if i not in additive_set and j not in additive_set]
```

### 8.2 S2 Quality Gate (Mandatory)

A LiCl-enabled S2 scan is **accepted only if**:

| Condition | Threshold |
|-----------|-----------|
| Forming bonds contain no Li/Cl atom | Must pass |
| LiCl remains within coordination envelope | Li–Cl ≤ 3.5 Å, Li–O ≤ 2.6 Å |
| Li–Cl distance does not exceed dissociation threshold | Li–Cl ≤ 3.5 Å |
| Li does not jump to different coordination site | `coordination_switch_allowed: false` |
| Cl does not approach reaction center | Cl–C_reactive ≥ 3.2 Å |

```python
def _check_s2_quality_gate(scan_result, la_additive, config) -> Tuple[bool, List[str]]:
    """Returns (pass, degraded_reasons)."""
    if not la_additive.enabled:
        return True, []
    
    violations = []
    trace = scan_result.li_coordination_trace
    
    if max(trace.li_cl_distances) > config["s2_quality_gate"]["li_cl_max"]:
        violations.append("LiCl_dissociation")
    if min(trace.li_o_distances) > config["s2_quality_gate"]["li_o_max"]:
        violations.append("LiCl_not_coordinated")
    if trace.coordination_switched and not config["s2_quality_gate"]["coordination_switch_allowed"]:
        violations.append("Li_coordination_switch")
    
    return len(violations) == 0, violations
```

If violations are found → scan status set to `degraded`, reasons logged.

### 8.3 Li–O / Li–Cl Distance Trace

The xTB scan module logs Li–O and Li–Cl distances at each scan step:

```python
# In XTBRunner.run_scan(), when additive atoms are present:
scan_result.li_coordination_trace = LiCoordinationTrace(
    step_indices=[0, 1, ..., n],
    li_cl_distances=[...],
    li_o_distances=[...],
    li_coordination_sites=[...],
)
```

---

## 9. S3 — TS Validation with Additive Mode Exclusion

### 9.1 Imaginary Mode Projection Check

After TS optimization, validate that the imaginary mode corresponds to the intended reaction, not LiCl motion:

```python
def _validate_ts_imaginary_mode(
    ts_log_path: Path,
    forming_bonds: Tuple[Tuple[int, int], ...],
    la_additive: LewisAcidAdditive,
    config: Dict[str, Any],
) -> Tuple[bool, str]:
    """
    Validate that the dominant imaginary mode projects onto forming bonds
    and not onto LiCl motion.
    """
    if not la_additive.enabled:
        return True, ""
    
    modes = extract_imaginary_modes(ts_log_path)
    if len(modes) == 0:
        return False, "no_imaginary_mode"
    
    # Skip modes dominated by additive atoms
    relevant_modes = []
    for mode in modes:
        projection_on_additive = mode.atom_projection(la_additive.additive_atom_indices)
        if projection_on_additive < config["s3_quality_gate"]["additive_mode_projection_max"]:
            relevant_modes.append(mode)
    
    if len(relevant_modes) == 0:
        return False, "all_modes_are_additive"
    
    # Primary mode must project onto forming bonds
    primary = relevant_modes[0]
    fb_projection = primary.bond_projection(forming_bonds)
    
    if fb_projection < config["s3_quality_gate"]["reaction_mode_projection_min"]:
        return False, f"low_reaction_mode_projection:{fb_projection:.2f}"
    
    return True, f"ok_projection:{fb_projection:.2f}"
```

### 9.2 S3 Quality Gate (Mandatory)

A LiCl-enabled TS is **accepted only if**:

| Condition | Threshold |
|-----------|-----------|
| ≥1 chemically relevant imaginary mode after excluding LiCl floppy modes | ≥ 1 |
| Dominant mode projects onto intended forming bonds | ≥ 0.35 |
| LiCl motion contributes below threshold | ≤ 0.25 |
| IRC for flagged cases | required if projection < 0.50 |

If gate fails → TS marked `degraded` with reason.

---

## 10. S4 — LiCl Feature Mining

### 10.1 Coordination Features

| Feature | Description |
|---------|-------------|
| `d_LiCl` | Li–Cl distance (Å) |
| `d_LiO_min` | Distance to nearest substrate O (Å) |
| `d_LiO_second` | Distance to second-nearest substrate O (Å) |
| `CN_Li_hetero` | Li coordination number (heteroatoms only) |
| `CN_Li_total` | Total Li coordination number |
| `angle_O_Li_O` | O–Li–O angle if 2+ O contacts (degrees) |
| `angle_Cl_Li_O` | Cl–Li–O angle (degrees) |
| `chelation_ring_size` | Size of chelation ring if applicable |

### 10.2 Activation Features

| Feature | Description |
|---------|-------------|
| `d_CO_elongation` | C=O bond elongation vs no-LA structure |
| `delta_q_carbonyl_O` | Carbonyl O NBO/Mulliken charge shift |
| `delta_q_carbonyl_C` | Carbonyl C NBO/Mulliken charge shift |
| `delta_Wiberg_CO` | C=O Wiberg bond order shift |
| `delta_dipole` | Total dipole moment shift |

### 10.3 Paired ΔFeatures (CRITICAL for ML)

The most important features are **differences from the no-LA baseline**:

```python
# For a reaction run with LiCl:
Δfeature = feature(substrate + LiCl) − feature(substrate_without_LiCl)

# ML uses Δfeatures to learn "how the system responds to chelation"
```

This requires running a matched no-LA reference for each reaction. For batch mode, the no-LA reference can be computed once per unique product SMILES and reused.

### 10.4 Quality Flags

| Flag | Source |
|------|--------|
| `LiCl_dissociation_flag` | S2 quality gate |
| `coordination_switch_flag` | S2 coordination trace |
| `Cl_near_reaction_center_flag` | S2 Cl–C distance |
| `additive_imaginary_mode_flag` | S3 mode projection |
| `ts_validated_flag` | S3 quality gate passed |

---

## 11. Validation Gates

### Gate A — Geometry Sanity

- Product/precursor + LiCl must form stable Li–O or Li–heteroatom contact
- LiCl must not drift away (> 3.5 Å Li–Cl or > 4.0 Å Li–substrate)
- Cl must not coordinate directly to reaction center (< 3.2 Å)

### Gate B — Organic Reaction Identity

- SMARTS-derived forming bonds (from canonical product SMILES) must match no-LA run
- Li/Cl must never enter forming-bond list
- TS imaginary mode must remain organic reaction mode

### Gate C — ML Usefulness

- LiCl features must be non-constant across dataset
- LiCl features must improve LORO/nested CV over no-LA baseline
- Improvement must disappear or shrink under shuffled target / shuffled LA control

### Gate D — Limited Physical Benchmark

- On 2–3 representative reactions: compare no-LA vs LiCl vs one real Lewis acid (at intermediate/product level only)
- Do **not** require exact energy match
- Require **same direction** of coordination/preorganization effect

---

## 12. Modified File Summary

| File | Change | Lines |
|------|--------|-------|
| `config/defaults.yaml` | Add `lewis_acid:` section (all config keys) | +40 |
| `rph_core/lewis_acid/__init__.py` | Create | +1 |
| `rph_core/lewis_acid/acidity.py` | Part A — descriptor registry | ~200 |
| `rph_core/lewis_acid/lewis_acid_registry.json` | Part A — output data | auto |
| `rph_core/orchestrator.py` | Add `LewisAcidAdditive`, `_prepare_lewis_acid_system()`, modify molecules dict, PipelineResult, pass additive to S2/S3, checkpoint hash | ~100 |
| `rph_core/steps/conformer_search/engine.py` | CREST NCI mode + wall potential for LA | ~20 |
| `rph_core/steps/step2_retro/retro_scanner.py` | Add `lewis_acid_additive` param, Li–O/Li–Cl trace, quality gate | ~50 |
| `rph_core/steps/step2_retro/smarts_matcher.py` | Build Mol from canonical SMILES, not XYZ; accept additive_indices | ~15 |
| `rph_core/steps/step3_opt/ts_optimizer.py` | Pass additive to validator | ~5 |
| `rph_core/steps/step3_opt/validator.py` | Imaginary mode projection check; skip additive modes | ~40 |
| `rph_core/utils/forming_bonds_resolver.py` | LA-hardened fallback: force SMARTS + element filter | +25 |
| **Total** | | **~500** |

---

## 13. Implementation Order (Phase 1)

```
Week 1:
  Day 1:  Config + LewisAcidAdditive dataclass + PipelineResult
  Day 2:  Orchestrator geometry preparation + checkpoint hashing
  Day 3:  S1 CREST NCI mode + post-S1 validation
  Day 4:  S2 forming bonds resolution + quality gate
  Day 5:  S3 mode projection check + quality gate

Week 2:
  Day 1:  S4 LiCl features (coordination + activation + Δfeatures)
  Day 2:  Part A descriptor registry (acetone affinity)
  Day 3:  Validation gates A–D
  Day 4:  Integration test + edge cases
  Day 5:  Documentation + cleanup
```

---

## 14. Risk Register (Final)

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| SMARTS poisoned by Li–O bond perception | Low | High | Build Mol from canonical SMILES, not XYZ |
| Li–Cl as false forming bond | Medium | High | Force SMARTS fallback; hard element filter |
| CREST dissociates LiCl | Medium | Medium | NCI mode + wall potential + weak restraint |
| xTB scan noisy with LiCl drift | Medium | Medium | Monitor trace; mark degraded on dissociation |
| TS optimizer chases LiCl modes | Medium | Medium | Validate mode projection; exclude additive modes |
| LiCl coordinates to wrong site | Medium | Medium | CREST ensemble finds global min; validate |
| Checkpoint identity confusion | Low | High | Hash includes canonical SMILES + full LA config |
| ML features not useful | Low | Medium | Gate C validates improvement over no-LA baseline |
| Real LA correlation not linear | Certain | Low | Descriptors are ML features, not ΔG‡ correctors |
