# step4_features/AGENTS.md

## OVERVIEW
S4: extracts geometric/electronic/thermochemical features from S1/S2/S3 QC outputs via a plugin-based extractor system. V6.2 adds mechanism-aware features (Step1 activation, Step2 cyclization) and Multiwfn integration.

## WHERE TO LOOK
| File | Role |
|------|------|
| `feature_miner.py:31` | `FeatureMiner` — discovers plugins, builds context, runs each, writes 3 output files |
| `extractors/` | 14 plugin files + `base.py` registry |
| `context.py` | `FeatureContext`, `PluginTrace`, `FeatureResult` → serialized to `feature_meta.json` |
| `schema.py` | `FIXED_COLUMNS`, `DEFAULT_MLR_COLUMNS`, CSV writers, `JOB_SPEC_REQUIRED_KEYS` |
| `status.py` | `FeatureResultStatus`, `FeatureStatus`, `aggregate_plugin_status()` |
| `path_accessor.py` | `PathAccessor` — path validation + fingerprinting used by extractors |
| `mech_packager.py` | Resolves mechanism context (1574 lines) when explicit inputs not provided |
| `fragment_extractor.py` | Fragment split helpers (1071 lines); centralized — do not reimplement in extractors |
| `log_parser.py` | Gaussian/ORCA log parsing helpers |
| `config/` | `multiwfn_recipes.yaml`, `nbo_templates.yaml` — loaded at runtime, never hardcode |
| `multiwfn/recipes.py` | Multiwfn recipe loader/renderer |
| `nbo/e2_parser.py` | NBO E(2) parsing helpers |

## PLUGIN INTERFACE CONTRACT
Subclass `BaseExtractor` (`extractors/base.py:40`):

```python
class MyExtractor(BaseExtractor):
    def get_plugin_name(self) -> str: return "my_plugin"
    def get_required_inputs(self) -> List[str]: return ["ts_xyz", "ts_fchk"]
    def extract(self, context: FeatureContext) -> Dict[str, Any]:
        # return flat dict of feature_name → value; use NaN for missing
        # append warnings to context.get_plugin_trace("my_plugin").warnings

register_extractor(MyExtractor())    # module-level; causes registration on import
```

**Discovery**: `extractors/__init__.py` imports all modules (side-effect registration). `FeatureMiner` calls `list_extractors()`. Filter via `config.step4.enabled_plugins`.

**Degradation**: catch domain errors inside `extract()`, return NaN, add structured warnings to trace — do NOT raise. `BaseExtractor.run()` catches exceptions and sets `trace.status = FAILED`.

**Job specs** (Phase B/C): append dicts to `trace.job_specs`. Required keys: `engine`, `command`, `workdir`, `input_files`, `expected_outputs`, `cache_key`. In extract-only runs `context.job_run_policy = "disallow"`.

## EXTRACTOR INVENTORY
| File | Feature prefix | Domain |
|------|---------------|--------|
| `thermo.py` | `thermo.*` | ΔG‡, ΔE, method, solvent, energy source |
| `geometry.py` | `geom.*` | r1, r2, asynch, rg_ts, close contacts, r_avg, dr |
| `qc_checks.py` | `qc.*` | Validation flags (has_gibbs, forming_bonds_valid…) |
| `ts_quality.py` | `ts.*` | n_imag, imag freq (cm⁻¹), dipole |
| `step1_activation.py` | `s1_*` | Boltzmann/Gibbs weighted energy, Nconf_eff, Sconf, leaving-group geometry, α-H gating |
| `step2_cyclization.py` | `s2_*` | CDFT (μ, η, ω, HOMO/LUMO), GEDT, forming-bond distances, TS validity flags |
| `multiwfn_features.py` | `mw_*` | Fukui f+/f−, dual descriptor, QTAIM BCPs — via `utils/multiwfn_runner.py` |
| `nbo_e2.py` | `nbo.e2.*` | NBO E(2) interactions (or job_specs in Phase B) |
| `nics.py` | `arom.*` | NICS aromaticity via Multiwfn (or job_specs) |
| `interaction_analysis.py` | `eda.*` | Fragment EDA (or job_specs for fragment SPs) |
| `precursor_geometry.py` | `prec_geom.*` | Precursor molecular geometry |
| `asm_enrichment.py` | `asm.*` | Distortion/interaction energies from S3 `enrichment.json` |
| `fmo_cdft_dipolar.py` | `fmo_cdft_dipolar.*` | HOMO/LUMO/gap/ω from S3 dipolar outputs |

## V6.2 FEATURE CONTRACT
- **Step1 activation** (P0): Boltzmann/Gibbs weighted energies, conformer entropy/flexibility, leaving-group geometry, α-H gating
- **Step2 cyclization** (P0): kinetics/thermochemistry, TS geometry, CDFT (Fukui/Dual Descriptor), GEDT (graph-partition), imaginary-freq validation
- **Multiwfn** (P2): Tier-1 Fukui/Dual Descriptor/QTAIM, non-interactive batch, failure-tolerant, caching
- **Meta**: `feature_meta.json` records `config_snapshot` + `provenance` + per-plugin `trace.plugins`

## OUTPUT CONTRACT (3 files always written)
```
S4_Data/
├── features_raw.csv     # FIXED_COLUMNS (in order) + all dynamic columns (alpha-sorted)
├── features_mlr.csv     # DEFAULT_MLR_COLUMNS or config.step4.mlr.columns
└── feature_meta.json    # meta.{schema_version, feature_status, config_snapshot, provenance}
                         # trace.{inputs_fingerprint, plugins.<name>.{status, runtime_ms,
                         #   missing_fields, missing_paths, errors, warnings, job_specs}}
                         # warnings (top-level deduped list)
```

Tests generating/asserting `feature_meta.json`: `test_s4_v62_final_verification.py`, `test_s4_meta_warnings_and_weights.py`, `test_s4_extractor_degrade_behavior.py`, `test_s4_no_qc_execution.py`.

## CONVENTIONS
- Extractors MUST use `context`/`path_accessor` for all paths — never hardcode directory names.
- Missing artifacts: produce NaN + record in `trace.warnings`; never silently succeed.
- Forming bonds from S2 (via `forming_bonds_resolver`) are the **sole authority** for fragment split — do not recompute in extractors.
- Multiwfn calls exclusively through `rph_core/utils/multiwfn_runner.py`; disable interactive input; handle timeout + failure.

## ANTI-PATTERNS
- Reimplementing fragment split logic in individual extractors — use `fragment_extractor.py`.
- Silent "success" when fchk/log is missing — must record NaN + warning.
- Hardcoding config parameters inside extractor — read from `context.config`.
