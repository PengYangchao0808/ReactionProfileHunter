from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Tuple

from rph_core.steps.contracts import Step2Artifacts, Step3Artifacts, Step4Artifacts
from rph_core.utils.data_types import MolIdx


def _unpack_s2_engine_result(
    result: Sequence[Any],
) -> Tuple[Path, Path, Path, Tuple[Tuple[int, int], ...], Path, str, str, Tuple[str, ...], Optional[Path]]:
    n = len(result)
    if n == 8:
        (
            ts_guess_xyz,
            _,  # position 1: duplicate of intermediate_xyz (legacy)
            intermediate_xyz,
            returned_forming_bonds,
            scan_profile_json,
            status,
            ts_guess_confidence,
            degraded_reasons,
        ) = result
        ts_guess_gau_xtb = None
    elif n == 9:
        (
            ts_guess_xyz,
            _,  # position 1: duplicate of intermediate_xyz (legacy)
            intermediate_xyz,
            returned_forming_bonds,
            scan_profile_json,
            status,
            ts_guess_confidence,
            degraded_reasons,
            ts_guess_gau_xtb,
        ) = result
    else:
        raise ValueError(f"Unexpected S2 engine return length: expected 8 or 9, got {n}")

    normalized_bonds: list[tuple[int, int]] = []
    for pair in returned_forming_bonds:
        if not isinstance(pair, (tuple, list)) or len(pair) != 2:
            raise ValueError(f"Invalid forming bond pair from S2 engine: {pair!r}")
        normalized_bonds.append((int(pair[0]), int(pair[1])))

    return (
        Path(ts_guess_xyz),
        Path(intermediate_xyz),  # retained for backward compat with callers expecting 2 values
        Path(intermediate_xyz),
        tuple(normalized_bonds),
        Path(scan_profile_json),
        str(status),
        str(ts_guess_confidence),
        tuple(str(x) for x in degraded_reasons),
        Path(ts_guess_gau_xtb) if ts_guess_gau_xtb else None,
    )


def _read_s1_provenance(work_dir: Path) -> Dict[str, Any]:
    import json

    provenance_file = Path(work_dir) / "S1_ConfGeneration" / "provenance.json"
    if not provenance_file.exists():
        return {}

    try:
        with open(provenance_file, encoding="utf-8") as f:
            data = json.load(f)
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def _write_step2_provenance(
    *,
    work_dir: Path,
    version: str,
    product_xyz: Path,
    product_xyz_hash: Optional[str],
    step2_signature: Dict[str, Any],
    reaction_profile: Optional[str],
    forming_bonds: list,
    scan_config: Dict[str, Any],
    path_search_mode: str,
    s2_strategy: str,
    metadata: Dict[str, Any],
) -> None:
    import json as _json
    s2_dir = work_dir / "S2_Retro"
    s2_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": "step2_provenance_v1",
        "rph_version": version,
        "product_xyz": str(product_xyz),
        "product_xyz_hash": product_xyz_hash,
        "step2_signature": step2_signature,
        "reaction_profile": reaction_profile,
        "forming_bonds": forming_bonds,
        "scan_config": {k: v for k, v in scan_config.items() if k != "output_dir"},
        "path_search": {"mode": path_search_mode},
        "s2_strategy": s2_strategy,
        "metadata": metadata,
    }
    prov_path = s2_dir / "step2_provenance.json"
    with open(prov_path, "w", encoding="utf-8") as f:
        _json.dump(payload, f, indent=2, ensure_ascii=False, default=str)


def _adapt_product_xyz_for_s2_if_needed(
    *,
    hunter: Any,
    work_dir: Path,
    product_xyz_file: Path,
) -> Path:
    import json
    import shutil

    step2_cfg = hunter.config.get("step2", {}) if isinstance(hunter.config, dict) else {}
    adapter_cfg = step2_cfg.get("pes_adapter", {}) if isinstance(step2_cfg, dict) else {}
    if not bool(adapter_cfg.get("enabled", False)):
        return product_xyz_file

    mode = str(adapter_cfg.get("mode", "fallback")).lower()
    if mode != "fallback":
        return product_xyz_file

    provenance = _read_s1_provenance(work_dir)
    has_opt = bool(provenance.get("has_geometry_optimization", False))
    if has_opt:
        return product_xyz_file

    protocol = str(provenance.get("protocol", "unknown"))
    provenance_schema_version = str(provenance.get("schema_version", "unknown"))
    trigger_class = "exception_path" if protocol in {"ext", "default", "lite", "zero"} else "legacy_or_unknown"

    fallback_dir = Path(work_dir) / "S2_Retro" / "pes_adapter_fallback"
    fallback_dir.mkdir(parents=True, exist_ok=True)
    adapted_xyz = fallback_dir / "product_relaxed.xyz"
    if not Path(product_xyz_file).exists():
        raise FileNotFoundError(f"PES adapter source not found: {product_xyz_file}")

    shutil.copy2(product_xyz_file, adapted_xyz)

    meta = {
        "mode": "fallback",
        "triggered": True,
        "reason": "s1_geometry_not_dft_optimized",
        "protocol": protocol,
        "provenance_schema_version": provenance_schema_version,
        "trigger_source": "s1_provenance.has_geometry_optimization=false",
        "trigger_class": trigger_class,
        "source": str(product_xyz_file),
        "adapted": str(adapted_xyz),
    }
    with open(fallback_dir / "adapter_meta.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)

    hunter.logger.warning(
        "[S2] S1 geometry is not DFT-optimized; applying fallback pes_adapter path"
    )
    return adapted_xyz


def run_step2(
    hunter: Any,
    product_xyz: Path,
    work_dir: Path,
    reaction_profile: Optional[str],
    cleaner_data: Optional[Dict[str, Any]],
) -> Step2Artifacts:
    import json
    import numpy as np
    from rph_core.utils.file_io import read_xyz
    from rph_core.utils.geometry_tools import GeometryUtils
    from rph_core.utils.scan_profile_plotter import plot_scan_profile
    product_xyz_file_raw = hunter._resolve_product_xyz_for_s2(product_xyz)
    product_xyz_file = _adapt_product_xyz_for_s2_if_needed(
        hunter=hunter,
        work_dir=work_dir,
        product_xyz_file=product_xyz_file_raw,
    )
    adapter_triggered = product_xyz_file != product_xyz_file_raw
    adapter_trigger_class = "unknown"
    if adapter_triggered:
        adapter_meta_file = work_dir / "S2_Retro" / "pes_adapter_fallback" / "adapter_meta.json"
        if adapter_meta_file.exists():
            try:
                adapter_meta = json.loads(adapter_meta_file.read_text(encoding="utf-8"))
                adapter_trigger_class = str(adapter_meta.get("trigger_class", "unknown"))
            except Exception:
                adapter_trigger_class = "unknown"
    profile_key = hunter._resolve_profile_key(reaction_profile=reaction_profile, cleaner_data=cleaner_data)
    forming_bonds = hunter._resolve_forming_bonds_for_s2(
        cleaner_data=cleaner_data,
        product_xyz_file=product_xyz_file,
        work_dir=work_dir,
    )
    hunter.logger.info(f"[S2] Using executable XYZ forming_bonds: {forming_bonds}")

    reaction_profiles = hunter.config.get("reaction_profiles", {}) if isinstance(hunter.config, dict) else {}
    profile_cfg = (
        reaction_profiles.get(str(profile_key), {}) if isinstance(reaction_profiles, dict) and profile_key else {}
    )
    step2_cfg = hunter.config.get("step2", {}) if isinstance(hunter.config, dict) else {}
    s2_scan_cfg = dict(step2_cfg.get("scan", {}) or {})
    if isinstance(profile_cfg, dict):
        s2_scan_cfg.update(dict(profile_cfg.get("scan", {}) or {}))
    s2_scan_cfg["output_dir"] = work_dir / "S2_Retro"

    step2_signature = hunter._build_step2_signature(
        work_dir=work_dir,
        product_xyz_file=product_xyz_file,
        forming_bonds=forming_bonds,
        reaction_profile=profile_key,
        scan_config=s2_scan_cfg,
    )

    path_search_cfg = step2_cfg.get("path_search", {}) or {}
    if path_search_cfg.get("enabled") is True:
        path_search_mode = "always"
    elif path_search_cfg.get("enabled") is False:
        path_search_mode = "never"
    else:
        path_search_mode = str(path_search_cfg.get("mode", "rescue")).strip().lower()

    s2_strategy = str(profile_cfg.get("s2_strategy", "retro_scan")).strip().lower()
    hunter.logger.debug(f"[S2] Ignoring configured s2_strategy={s2_strategy!r}; retro_scan is the only supported workflow")
    hunter.logger.info("[S2] Running retro_scan workflow")
    engine_result = hunter.s2_engine.run_retro_scan(
        product_xyz=product_xyz_file,
        output_dir=s2_scan_cfg["output_dir"],
        forming_bonds=forming_bonds,
        scan_config=s2_scan_cfg,
    )
    (
        ts_guess_xyz,
        _,
        intermediate_xyz,
        returned_forming_bonds,
        scan_profile_json,
        status,
        ts_guess_confidence,
        degraded_reasons,
        ts_guess_gau_xtb,
    ) = _unpack_s2_engine_result(engine_result)
    generation_method = "retro_scan"
    
    # Capture S2.1 results BEFORE path_search potentially overwrites scan_profile_json
    s1_scan_profile_json = scan_profile_json
    s1_gau_xtb_xyz = ts_guess_gau_xtb
    s1_gau_xtb_energy = None
    s1_gau_xtb_distance = None
    
    # Calculate S2.1 Gau_XTB energy and distance if available
    if s1_gau_xtb_xyz and Path(s1_gau_xtb_xyz).exists():
        try:
            from rph_core.utils.gau_xtb_interface import GauXTBInterface
            gxtb = GauXTBInterface(config=hunter.config, nproc=1)
            log_file = Path(s1_gau_xtb_xyz).parent / "attempt_1" / "input.log"
            if log_file.exists():
                s1_gau_xtb_energy = gxtb._parse_energy(log_file)
            if returned_forming_bonds:
                coords, _ = read_xyz(s1_gau_xtb_xyz)
                dist_candidates = []
                for bond in returned_forming_bonds:
                    idx0 = int(bond[0])
                    idx1 = int(bond[1])
                    dist = GeometryUtils.calculate_distance(coords, idx0, idx1)
                    dist_candidates.append(dist)
                if dist_candidates:
                    s1_gau_xtb_distance = min(dist_candidates)
        except Exception as e:
            hunter.logger.warning(f"[S2] Failed to get S2.1 Gau_XTB data: {e}")
    
    retro_scan_profile_path = scan_profile_json
    retro_scan_energies_hartree = None
    ts_distance_from_scan = None
    dipole_distance_from_scan = None
    if scan_profile_json.exists():
        try:
            with open(scan_profile_json) as f:
                retro_data = json.load(f)
            retro_scan_energies_hartree = retro_data.get("energies_hartree")
            knee_data = retro_data.get("knee_point_algorithm", {})
            ts_distance_from_scan = knee_data.get("ts_distance")
            dipole_distance_from_scan = knee_data.get("dipole_distance")
        except Exception:
            pass

    ts_guess_gau_xtb_from_path = None
    returned_forming_bonds_from_path = None

    should_run_path_search = False
    if path_search_mode == "always":
        should_run_path_search = True
    elif path_search_mode == "rescue":
        retro_status_upper = str(status).upper() if status else "UNKNOWN"
        should_run_path_search = retro_status_upper in ("DEGRADED", "FAILED", "UNKNOWN")
        if should_run_path_search:
            hunter.logger.info(
                f"[S2] S2.1 retro_scan status: {retro_status_upper} — "
                f"triggering path_search as rescue"
            )
    elif path_search_mode == "never":
        should_run_path_search = False

    if should_run_path_search:
        reason = "rescue (retro_scan degraded)" if path_search_mode == "rescue" else "always-enabled"
        hunter.logger.info(f"[S2] S2.2: Running path search ({reason})")
        
        (
            ts_guess_from_path,
            reactant_xyz_from_path,
            reactant_complex_xyz_from_path,
            returned_forming_bonds_from_path,
            scan_profile_json_from_path,
            status_from_path,
            ts_guess_confidence_from_path,
            degraded_reasons_from_path,
            ts_guess_gau_xtb_from_path,
        ) = hunter.s2_engine.run_path_search(
            start_xyz=intermediate_xyz,
            end_xyz=product_xyz_file,
            output_dir=s2_scan_cfg["output_dir"],
            forming_bonds=forming_bonds,
        )
        
        ts_guess_xyz = ts_guess_from_path
        scan_profile_json = scan_profile_json_from_path
        status = status_from_path
        ts_guess_confidence = ts_guess_confidence_from_path
        degraded_reasons = degraded_reasons_from_path
        generation_method = "xtb_path_search"
        
        path_ts_distance = None
        if ts_guess_from_path and Path(ts_guess_from_path).exists() and returned_forming_bonds_from_path:
            try:
                coords, symbols = read_xyz(ts_guess_from_path)
                if coords is not None and len(returned_forming_bonds_from_path) > 0:
                    dist_candidates = []
                    for bond in returned_forming_bonds_from_path:
                        idx0 = int(bond[0])
                        idx1 = int(bond[1])
                        dist = GeometryUtils.calculate_distance(coords, idx0, idx1)
                        dist_candidates.append((dist, idx0, idx1))

                    if dist_candidates:
                        path_ts_distance, idx0, idx1 = min(dist_candidates, key=lambda x: x[0])
                        hunter.logger.info(f"[S2] Path TS bond distance: {path_ts_distance:.3f} Å (using xyz indices {idx0}-{idx1})")
                    
                    if scan_profile_json.exists():
                        with open(scan_profile_json) as f:
                            profile_data = json.load(f)
                        
                        profile_data["path_ts_distance"] = path_ts_distance
                        profile_data["energies_hartree"] = retro_scan_energies_hartree
                        profile_data["knee_point_algorithm"] = {
                            "ts_distance": ts_distance_from_scan,
                            "dipole_distance": dipole_distance_from_scan,
                        }
                        
                        with open(scan_profile_json, "w") as f:
                            json.dump(profile_data, f, indent=2)
            except Exception as e:
                hunter.logger.warning(f"[S2] Failed to calculate path TS distance: {e}")
    
    # Final plotting: use S2.1 scan_profile which has energies_hartree, add both Gau_XTB results
    if s1_scan_profile_json and s1_scan_profile_json.exists():
        try:
            with open(s1_scan_profile_json) as f:
                profile_data = json.load(f)
            
            # Check if we have energies_hartree (from S2.1 retro_scan)
            if profile_data.get("energies_hartree"):
                # Get S2.2 Gau_XTB data
                s2_gau_xtb_xyz = ts_guess_gau_xtb_from_path if should_run_path_search else None
                s2_gau_xtb_energy = None
                s2_gau_xtb_distance = None
                
                if s2_gau_xtb_xyz and Path(s2_gau_xtb_xyz).exists():
                    try:
                        from rph_core.utils.gau_xtb_interface import GauXTBInterface
                        gxtb = GauXTBInterface(config=hunter.config, nproc=1)
                        log_file = Path(s2_gau_xtb_xyz).parent / "attempt_1" / "input.log"
                        if log_file.exists():
                            s2_gau_xtb_energy = gxtb._parse_energy(log_file)
                        if returned_forming_bonds_from_path:
                            coords, _ = read_xyz(s2_gau_xtb_xyz)
                            dist_candidates = []
                            for bond in returned_forming_bonds_from_path:
                                idx0 = int(bond[0])
                                idx1 = int(bond[1])
                                dist = GeometryUtils.calculate_distance(coords, idx0, idx1)
                                dist_candidates.append(dist)
                            if dist_candidates:
                                s2_gau_xtb_distance = min(dist_candidates)
                    except Exception:
                        pass
                
                # Update scan_profile with both Gau_XTB results
                profile_data["gau_xtb_s2.1_energy"] = s1_gau_xtb_energy
                profile_data["gau_xtb_s2.1_distance"] = s1_gau_xtb_distance
                profile_data["gau_xtb_s2.2_energy"] = s2_gau_xtb_energy
                profile_data["gau_xtb_s2.2_distance"] = s2_gau_xtb_distance
                
                with open(s1_scan_profile_json, "w") as f:
                    json.dump(profile_data, f, indent=2)
                
                path_ts_distance_for_plot = profile_data.get("path_ts_distance")
                plot_scan_profile(
                    s1_scan_profile_json,
                    energy_unit="kcal",
                    show_peak=True,
                    ts_distance=ts_distance_from_scan,
                    dipole_distance=dipole_distance_from_scan,
                    path_ts_distance=path_ts_distance_for_plot,
                    path_ts_energy=None,
                    gau_xtb_distance=s1_gau_xtb_distance,
                    gau_xtb_energy=s1_gau_xtb_energy,
                    gau_xtb2_distance=s2_gau_xtb_distance,
                    gau_xtb2_energy=s2_gau_xtb_energy,
                )
                hunter.logger.info(f"[S2] Final scan profile plot saved with Gau_XTB results")
        except Exception as e:
            hunter.logger.warning(f"[S2] Failed to create final plot: {e}")

    fallback_reasons: Tuple[str, ...] = tuple()
    if adapter_triggered:
        fallback_reasons = ("pes_adapter_fallback_applied",)
        if adapter_trigger_class == "exception_path":
            fallback_reasons += ("pes_adapter_fallback_exception_path",)
        elif adapter_trigger_class == "legacy_or_unknown":
            fallback_reasons += ("pes_adapter_fallback_legacy_or_unknown",)

    _write_step2_provenance(
        work_dir=work_dir,
        version=hunter.config.get("_rph_version") or __import__("rph_core.version", fromlist=["__version__"]).__version__,
        product_xyz=product_xyz_file,
        product_xyz_hash=hunter._build_step2_signature(
            work_dir=work_dir,
            product_xyz_file=product_xyz_file,
            forming_bonds=forming_bonds,
            reaction_profile=profile_key,
            scan_config=s2_scan_cfg,
        ).get("product_xyz_hash") if step2_signature else None,
        step2_signature=step2_signature,
        reaction_profile=profile_key,
        forming_bonds=list(forming_bonds),
        scan_config=s2_scan_cfg,
        path_search_mode=path_search_mode,
        s2_strategy=s2_strategy,
        metadata={
            "generation_method": generation_method,
            "status": status,
            "ts_guess_confidence": ts_guess_confidence,
            "degraded_reasons": list(degraded_reasons) + list(fallback_reasons),
        },
    )

    return Step2Artifacts(
        ts_guess_xyz=ts_guess_xyz,
        intermediate_xyz=intermediate_xyz,
        forming_bonds=tuple((MolIdx(int(i)), MolIdx(int(j))) for (i, j) in returned_forming_bonds),
        generation_method=generation_method,
        status=status,
        ts_guess_confidence=ts_guess_confidence,
        degraded_reasons=tuple(degraded_reasons) + fallback_reasons,
        step2_signature=step2_signature,
        scan_profile_json=scan_profile_json,
    )


def run_step3(
    hunter: Any,
    ts_guess_xyz: Path,
    intermediate_xyz: Path,
    product_xyz: Path,
    work_dir: Path,
    e_product_l2: Optional[float],
    product_thermo: Optional[Path],
    forming_bonds,
    old_checkpoint: Optional[Path],
) -> Step3Artifacts:
    s3_result = hunter.s3_engine.run(
        ts_guess=ts_guess_xyz,
        intermediate=intermediate_xyz,
        product=product_xyz,
        output_dir=work_dir / "S3_TS",
        e_product_l2=e_product_l2,
        product_thermo=product_thermo,
        forming_bonds=forming_bonds,
        old_checkpoint=old_checkpoint,
    )

    return Step3Artifacts(
        ts_final_xyz=s3_result.ts_final_xyz,
        sp_report=s3_result.sp_report,
        ts_fchk=s3_result.ts_fchk,
        ts_log=s3_result.ts_log,
        ts_qm_output=s3_result.ts_qm_output,
        intermediate_fchk=s3_result.intermediate_fchk,
        intermediate_log=s3_result.intermediate_log,
        intermediate_qm_output=s3_result.intermediate_qm_output,
        # V7.1: forward S3-optimized intermediate data
        intermediate_xyz=getattr(s3_result, 'intermediate_xyz', None),
        intermediate_l2_energy=getattr(s3_result, 'intermediate_l2_energy', None),
        intermediate_opt_output=getattr(s3_result, 'intermediate_opt_output', None),
        intermediate_sp_output=getattr(s3_result, 'intermediate_sp_output', None),
    )


def run_step4(
    hunter: Any,
    ts_final_xyz: Path,
    intermediate_xyz: Path,
    product_xyz: Path,
    work_dir: Path,
    forming_bonds,
    sp_matrix_report,
    ts_fchk: Optional[Path],
    intermediate_fchk: Optional[Path],
    product_fchk: Optional[Path],
    ts_log: Optional[Path],
    intermediate_log: Optional[Path],
    product_log: Optional[Path],
    ts_qm_output: Optional[Path],
    intermediate_qm_output: Optional[Path],
    product_qm_output: Optional[Path],
    product_smiles: Optional[str] = None,
    s3_intermediate_xyz: Optional[Path] = None,
    s3_intermediate_l2_energy: Optional[float] = None,
) -> Step4Artifacts:
    raise ImportError(
        "Step 4 feature extraction has been moved to RPH_Postprocess/rph_features/. "
        "This function is no longer functional.\n"
        "Use the standalone CLI:\n"
        "  rph-features extract --rph-run <work_dir> --output <features_dir>"
    )
