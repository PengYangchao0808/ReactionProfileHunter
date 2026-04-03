from pathlib import Path
from typing import Any, Dict, Optional

from rph_core.steps.contracts import Step2Artifacts, Step3Artifacts, Step4Artifacts


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
    product_xyz_file = hunter._resolve_product_xyz_for_s2(product_xyz)
    profile_key = hunter._resolve_profile_key(reaction_profile=reaction_profile, cleaner_data=cleaner_data)
    forming_bonds = hunter._resolve_forming_bonds_for_s2(
        cleaner_data=cleaner_data,
        product_xyz_file=product_xyz_file,
        work_dir=work_dir,
    )
    hunter.logger.info(f"[S2] Using MolIdx forming_bonds: {forming_bonds}")

    forming_bonds_map = None
    atom_map = None
    if cleaner_data:
        raw_map_pairs = cleaner_data.get("formed_bond_map_pairs") or cleaner_data.get("raw", {}).get("formed_bond_map_pairs")
        if raw_map_pairs:
            from rph_core.utils.cleaner_adapter import parse_pairs_text
            forming_bonds_map = tuple(tuple(x) for x in parse_pairs_text(str(raw_map_pairs)))
            
        mapped_smiles = cleaner_data.get("mapped_product_smiles") or cleaner_data.get("raw", {}).get("mapped_product_smiles")
        if mapped_smiles and product_xyz_file.exists():
            from rph_core.utils.cleaner_adapter import get_map_to_xyz_dict
            atom_map = get_map_to_xyz_dict(mapped_smiles, product_xyz_file)
        hunter.logger.info(f"[S2] Using MapId forming_bonds: {forming_bonds_map}")

    s2_scan_cfg = hunter._resolve_forward_scan_config(reaction_profile=reaction_profile, cleaner_data=cleaner_data)
    s2_scan_cfg["output_dir"] = work_dir / "S2_Retro"

    reaction_profiles = hunter.config.get("reaction_profiles", {}) if isinstance(hunter.config, dict) else {}
    profile_cfg = (
        reaction_profiles.get(str(profile_key), {}) if isinstance(reaction_profiles, dict) and profile_key else {}
    )
    step2_cfg = hunter.config.get("step2", {}) if isinstance(hunter.config, dict) else {}

    step2_signature = hunter._build_step2_signature(
        work_dir=work_dir,
        product_xyz_file=product_xyz_file,
        forming_bonds=forming_bonds,
        reaction_profile=profile_key,
        scan_config=s2_scan_cfg,
    )

    path_search_enabled = step2_cfg.get("path_search", {}).get("enabled", False)

    hunter.logger.info("[S2] S2.1: Running stretch search (retro_scan) to generate dipole intermediate")
    
    (
        ts_guess_from_scan,
        substrate_xyz,
        intermediate_xyz,
        returned_forming_bonds,
        scan_profile_json,
        status,
        ts_guess_confidence,
        degraded_reasons,
        ts_guess_gau_xtb,
    ) = hunter.s2_engine.run_retro_scan(
        product_xyz=product_xyz_file,
        output_dir=s2_scan_cfg["output_dir"],
        forming_bonds=forming_bonds_map if forming_bonds_map else forming_bonds,
        scan_config=s2_scan_cfg,
        atom_map=atom_map,
    )
    
    ts_guess_xyz = ts_guess_from_scan
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

    if path_search_enabled:
        hunter.logger.info("[S2] S2.2: Running path search (xtb_path_search) using dipole intermediate")
        
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
            forming_bonds=forming_bonds_map if forming_bonds_map else forming_bonds,
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
                        if atom_map and len(atom_map) > 0:
                            mapped0 = atom_map.get(idx0, idx0)
                            mapped1 = atom_map.get(idx1, idx1)
                            if mapped0 is None or mapped1 is None:
                                continue
                            idx0 = int(mapped0)
                            idx1 = int(mapped1)
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
                s2_gau_xtb_xyz = ts_guess_gau_xtb_from_path if path_search_enabled else None
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
                                if atom_map:
                                    mapped0 = atom_map.get(idx0, idx0)
                                    mapped1 = atom_map.get(idx1, idx1)
                                    if mapped0 is None or mapped1 is None:
                                        continue
                                    idx0 = int(mapped0)
                                    idx1 = int(mapped1)
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

    return Step2Artifacts(
        ts_guess_xyz=ts_guess_xyz,
        substrate_xyz=substrate_xyz,
        intermediate_xyz=intermediate_xyz,
        forming_bonds=tuple(returned_forming_bonds),
        forming_bonds_map=None,
        generation_method=generation_method,
        status=status,
        ts_guess_confidence=ts_guess_confidence,
        degraded_reasons=tuple(degraded_reasons),
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
    )


def run_step4(
    hunter: Any,
    ts_final_xyz: Path,
    substrate_xyz: Path,
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
) -> Step4Artifacts:
    s1_artifacts = hunter._resolve_s1_artifacts(work_dir)
    features_csv = hunter.s4_engine.run(
        ts_final=ts_final_xyz,
        reactant=substrate_xyz,
        product=product_xyz,
        output_dir=work_dir / "S4_Data",
        s1_dir=s1_artifacts.get("s1_dir"),
        s1_shermo_summary_file=s1_artifacts.get("s1_shermo_summary_file"),
        s1_hoac_thermo_file=s1_artifacts.get("s1_hoac_thermo_file"),
        s1_conformer_energies_file=s1_artifacts.get("s1_conformer_energies_file"),
        s1_precursor_xyz=s1_artifacts.get("s1_precursor_xyz"),
        forming_bonds=forming_bonds,
        sp_matrix_report=sp_matrix_report,
        ts_fchk=ts_fchk,
        reactant_fchk=intermediate_fchk,
        product_fchk=product_fchk,
        ts_log=ts_log,
        reactant_log=intermediate_log,
        product_log=product_log,
        ts_orca_out=ts_qm_output,
        reactant_orca_out=intermediate_qm_output,
        product_orca_out=product_qm_output,
    )
    return Step4Artifacts(features_csv=features_csv)
