# extractors/AGENTS.md

## OVERVIEW
Plugin registry + 14 extractor implementations. Each module registers an instance of its `BaseExtractor` subclass at import time. `FeatureMiner` discovers all registered extractors via `list_extractors()`.

## WHERE TO LOOK
| File | Class | Feature prefix |
|------|-------|---------------|
| `base.py:40` | `BaseExtractor` | Registry: `register_extractor()`, `list_extractors()` |
| `thermo.py` | `ThermoExtractor` | `thermo.*` — ΔG‡/ΔE/ΔG_rxn, method, solvent; Gibbs priority with Shermo fallback |
| `geometry.py` | `GeometryExtractor` | `geom.*` — forming bond r1/r2, asynch index, rg_ts, close contacts density |
| `qc_checks.py` | `QCChecksExtractor` | `qc.*` — validation flags: has_gibbs, sp_report_validated, forming_bonds_valid |
| `ts_quality.py` | `TSQualityExtractor` | `ts.*` — n_imag, imag1_cm1_abs, dipole_debye; warns if n_imag ≠ 1 |
| `step1_activation.py` | `Step1ActivationExtractor` | `s1_*` — Boltzmann/Gibbs weighted energy, Nconf_eff, Sconf, leaving-group geom, α-H |
| `step2_cyclization.py` | `Step2CyclizationExtractor` | `s2_*` — CDFT (μ/η/ω), GEDT, s2_d_forming_1/2, TS validity, imag freq |
| `multiwfn_features.py` | `MultiwfnFeaturesExtractor` | `mw_*` — Fukui f+/f−, dual descriptor, QTAIM BCP ρ/∇²ρ; mw_cache_hit |
| `nbo_e2.py` | `NBOE2Extractor` | `nbo.e2.*` — NBO E(2) interactions; generates job_specs in Phase B |
| `nics.py` | `NICSExtractor` | `arom.*` — NICS aromaticity via Multiwfn; generates job_specs in Phase B |
| `interaction_analysis.py` | `InteractionAnalysisExtractor` | `eda.*` — fragment EDA; generates job_specs for fragment SPs |
| `precursor_geometry.py` | `PrecursorGeometryExtractor` | `prec_geom.*` — precursor molecular geometry |
| `asm_enrichment.py` | `ASMEnrichmentExtractor` | `asm.*` — distortion_total_kcal, interaction_kcal from S3 `enrichment.json` |
| `fmo_cdft_dipolar.py` | `FmoCdftDipolarParser` | `fmo_cdft_dipolar.*` — HOMO/LUMO/gap/ω from S3 dipolar outputs |

## HOW TO ADD A PLUGIN
1. Create `extractors/my_plugin.py`
2. Implement `class MyExtractor(BaseExtractor)` with `get_plugin_name()`, `get_required_inputs()`, `extract(context) -> Dict`
3. Add `register_extractor(MyExtractor())` at module level
4. Import your module in `extractors/__init__.py` so registration fires on package load
5. Add test coverage to `tests/test_s4_no_qc_execution.py` or similar

## BASE EXTRACTOR CONTRACT
```python
# Required
get_plugin_name(self) -> str          # unique name key
get_required_inputs(self) -> List[str] # checked by validate_inputs()
extract(self, context: FeatureContext) -> Dict[str, Any]  # returns flat feature dict

# Optional
get_required_inputs_for_context(self, context) -> List[str]  # dynamic requirements
can_submit_jobs(self) -> bool  # default False; True = plugin may populate trace.job_specs
```

`run(context)` in `BaseExtractor` handles: `validate_inputs()` → `extract()` → populate `PluginTrace` (status, runtime_ms, errors, warnings, job_specs, _extracted_features). Exceptions are caught and recorded as `FAILED`.

## JOB SPEC KEYS (for Phase B/C plugins)
Required: `engine`, `command`, `workdir`, `input_files`, `expected_outputs`, `cache_key`
Optional: `stdin_lines`, `recipe_name`, `recipe_version`, `template_yaml_fingerprint`

## ANTI-PATTERNS
- Raising uncaught exceptions in `extract()` — catch domain errors, return NaN + warnings.
- Calling `subprocess` or writing files in `extract()` — use `context` and `trace.job_specs`.
- Using paths from `context` by string manipulation — use `PathAccessor` methods.
- Reimplementing fragment split — use `fragment_extractor.py`.
