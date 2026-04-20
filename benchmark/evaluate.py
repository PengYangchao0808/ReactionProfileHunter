#!/usr/bin/env python3
# pyright: reportDeprecated=false, reportUnusedCallResult=false, reportAny=false, reportExplicitAny=false, reportUnknownMemberType=false, reportUnknownVariableType=false, reportUnknownArgumentType=false, reportUnknownLambdaType=false, reportReturnType=false, reportOptionalOperand=false
"""Post-run evaluator for S1 conformer-search benchmark outputs."""

import argparse
import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional


HARTREE_TO_KCAL = 627.509474
DEFAULT_PROTOCOLS = ["ext", "full", "lite", "zero"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate S1 benchmark outputs across protocols.",
    )
    parser.add_argument(
        "--rx-ids",
        nargs="+",
        required=True,
        help="One or more reaction ids, e.g. 1 3 8",
    )
    parser.add_argument(
        "--protocols",
        nargs="+",
        default=list(DEFAULT_PROTOCOLS),
        help="Protocols to evaluate (default: ext full lite zero)",
    )
    parser.add_argument(
        "--output-base",
        default="./Output/benchmark",
        help="Base output directory (default: ./Output/benchmark)",
    )
    return parser.parse_args()


def read_text_file(path: Path) -> str:
    data = path.read_bytes()
    for encoding in ("utf-8-sig", "utf-8", "utf-16", "utf-16-le", "utf-16-be", "latin-1"):
        try:
            return data.decode(encoding)
        except UnicodeDecodeError:
            continue
    return data.decode("utf-8", errors="replace")


def load_json_file(path: Path) -> Any:
    return json.loads(read_text_file(path))


def safe_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def safe_int(value: Any) -> Optional[int]:
    parsed = safe_float(value)
    if parsed is None:
        return None
    return int(parsed)


def nested_get(data: Dict[str, Any], keys: Iterable[str], default: Any = None) -> Any:
    current = data
    for key in keys:
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current


def find_named_files(root: Path, filename: str) -> List[Path]:
    if not root.exists():
        return []
    try:
        return sorted((path for path in root.rglob(filename) if path.is_file()), key=lambda p: str(p))
    except OSError:
        return []


def locate_s1_dir_from_file(path: Path) -> Optional[Path]:
    for parent in (path,) + tuple(path.parents):
        if parent.name == "S1_ConfGeneration":
            return parent
    return None


def find_s1_provenance(run_dir: Path) -> Optional[Path]:
    for path in find_named_files(run_dir, "provenance.json"):
        if locate_s1_dir_from_file(path) is not None:
            return path
    return None


def pick_preferred_file(run_dir: Path, s1_dir: Optional[Path], filename: str) -> Optional[Path]:
    if s1_dir is not None:
        matches = find_named_files(s1_dir, filename)
        if matches:
            return matches[0]
    matches = find_named_files(run_dir, filename)
    return matches[0] if matches else None


def estimate_wall_time_seconds(paths: List[Path]) -> Optional[float]:
    timestamps = []
    for path in paths:
        try:
            stat = path.stat()
        except OSError:
            continue
        timestamps.extend([stat.st_ctime, stat.st_mtime])
    if not timestamps:
        return None
    return max(timestamps) - min(timestamps)


def parse_conformer_energies(path: Optional[Path]) -> List[float]:
    if path is None or not path.exists():
        return []
    try:
        raw = load_json_file(path)
    except (OSError, ValueError, json.JSONDecodeError):
        return []
    if not isinstance(raw, list):
        return []
    values = []
    for item in raw:
        value = safe_float(item)
        if value is not None:
            values.append(value)
    return values


def choose_gibbs_hartree(row: Dict[str, Any]) -> Optional[float]:
    """Try multiple column names for Gibbs free energy (handles schema evolution)."""
    for key in ("g_used", "g_hartree", "g_sum", "g_conc"):
        value = safe_float(row.get(key))
        if value is not None:
            return value
    return None


def parse_conformer_thermo(path: Optional[Path]) -> Dict[str, Optional[float]]:
    result = {
        "best_sp_energy_hartree": None,
        "best_gibbs_kcal": None,
    }
    if path is None or not path.exists():
        return result
    try:
        text = read_text_file(path)
    except OSError:
        return result
    rows = list(csv.DictReader(text.splitlines()))
    if not rows:
        return result

    # Detect which column holds SP energy (handles schema evolution)
    sample_row = rows[0]
    sp_col = None
    for col in ("sp_energy", "sp_energy_hartree"):
        if col in sample_row:
            sp_col = col
            break

    sp_values = []
    gibbs_values = []
    for row in rows:
        if sp_col is not None:
            sp_value = safe_float(row.get(sp_col))
            if sp_value is not None:
                sp_values.append(sp_value)
        gibbs_value = choose_gibbs_hartree(row)
        if gibbs_value is not None:
            gibbs_values.append(gibbs_value)

    if sp_values:
        result["best_sp_energy_hartree"] = min(sp_values)
    if gibbs_values:
        result["best_gibbs_kcal"] = min(gibbs_values) * HARTREE_TO_KCAL
    return result


def build_base_metrics(protocol: str) -> Dict[str, Any]:
    return {
        "status": "failed",
        "protocol": protocol,
        "provenance_found": False,
        "product_xyz_found": False,
        "two_stage_enabled": None,
        "freq_requested": None,
        "final_sp_requested": None,
        "selection_mode": None,
        "handoff_mode_requested": None,
        "handoff_mode_effective": None,
        "fallback_triggered": None,
        "funnel_candidates": None,
        "funnel_survivors": None,
        "selected_candidates": None,
        "final_dft_count": None,
        "best_sp_energy_hartree": None,
        "best_gibbs_kcal": None,
        "energy_range_kcal": None,
        "mean_energy_kcal": None,
        "product_energy_hartree": None,
        "molecule_names": [],
        "per_molecule": {},
        "wall_time_seconds": None,
        "paths": {
            "run_dir": None,
            "s1_dir": None,
            "provenance_json": None,
            "product_min_xyz": None,
            "conformer_energies_json": None,
            "conformer_thermo_csv": None,
        },
        "notes": [],
    }


def evaluate_molecule(mol_dir: Path, _mol_name: str) -> Dict[str, Any]:
    result: Dict[str, Any] = {
        "status": "failed",
        "final_dft_count": None,
        "best_sp_energy_hartree": None,
        "best_gibbs_kcal": None,
        "energy_range_kcal": None,
        "is_small_molecule": None,
        "mean_energy_kcal": None,
    }
    if not mol_dir.exists():
        result["status"] = "missing_dir"
        return result

    dft_dir = mol_dir / "finalDFT"
    if not dft_dir.exists():
        dft_dir = mol_dir / "dft"
    result["is_small_molecule"] = not ((mol_dir / "crest").exists() or (mol_dir / "xtb2").exists())

    energies = parse_conformer_energies(dft_dir / "conformer_energies.json")
    if energies:
        result["final_dft_count"] = len(energies)
        result["energy_range_kcal"] = max(energies) - min(energies) if len(energies) > 1 else 0.0
        result["best_sp_energy_hartree"] = min(energies) / HARTREE_TO_KCAL
        result["mean_energy_kcal"] = sum(energies) / len(energies)

    thermo = parse_conformer_thermo(dft_dir / "conformer_thermo.csv")
    if thermo["best_gibbs_kcal"] is not None:
        result["best_gibbs_kcal"] = thermo["best_gibbs_kcal"]

    result["status"] = "ok"
    return result


def evaluate_protocol(rx_dir: Path, protocol: str) -> Dict[str, Any]:
    metrics = build_base_metrics(protocol)
    run_dir = rx_dir / protocol
    metrics["paths"]["run_dir"] = str(run_dir)

    provenance_path = find_s1_provenance(run_dir)
    s1_dir = locate_s1_dir_from_file(provenance_path) if provenance_path is not None else None
    product_xyz_path = pick_preferred_file(run_dir, s1_dir, "product_min.xyz")

    metrics["provenance_found"] = provenance_path is not None
    metrics["product_xyz_found"] = product_xyz_path is not None
    metrics["paths"]["s1_dir"] = str(s1_dir) if s1_dir is not None else None
    metrics["paths"]["provenance_json"] = str(provenance_path) if provenance_path is not None else None
    metrics["paths"]["product_min_xyz"] = str(product_xyz_path) if product_xyz_path is not None else None

    time_paths = [path for path in (run_dir, s1_dir, provenance_path, product_xyz_path) if path is not None]
    metrics["wall_time_seconds"] = estimate_wall_time_seconds(time_paths)

    if provenance_path is None:
        metrics["notes"].append("Missing provenance.json under S1_ConfGeneration.")
        return metrics

    try:
        provenance = load_json_file(provenance_path)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        metrics["notes"].append("Failed to parse provenance.json: {0}".format(exc))
        return metrics

    if not isinstance(provenance, dict):
        metrics["notes"].append("provenance.json does not contain a JSON object.")
        return metrics

    metrics["two_stage_enabled"] = nested_get(provenance, ["conformer_search", "two_stage_enabled"])
    metrics["freq_requested"] = nested_get(provenance, ["final_opt_sp", "freq_requested"])
    metrics["final_sp_requested"] = nested_get(provenance, ["final_opt_sp", "final_sp_requested"])
    metrics["selection_mode"] = nested_get(provenance, ["final_opt_sp", "selection_mode"])
    metrics["handoff_mode_requested"] = nested_get(provenance, ["handoff", "mode_requested"], nested_get(provenance, ["handoff", "mode"]))
    metrics["handoff_mode_effective"] = nested_get(provenance, ["handoff", "mode_effective"], nested_get(provenance, ["handoff", "mode"]))
    metrics["fallback_triggered"] = nested_get(provenance, ["handoff", "fallback_triggered"])
    metrics["funnel_candidates"] = safe_int(nested_get(provenance, ["funnel", "candidate_total"]))
    metrics["funnel_survivors"] = safe_int(nested_get(provenance, ["funnel", "survivor_count"]))
    metrics["selected_candidates"] = safe_int(
        nested_get(provenance, ["handoff", "selected_candidate_count"], nested_get(provenance, ["final_opt_sp", "selected_candidate_count"]))
    )

    molecules = provenance.get("molecules", {})
    if not isinstance(molecules, dict):
        molecules = {}

    mol_names = list(molecules.keys())
    metrics["molecule_names"] = mol_names

    per_molecule = {}
    for mol_name in mol_names:
        mol_provenance = molecules.get(mol_name, {})
        if not isinstance(mol_provenance, dict):
            mol_provenance = {}
        mol_dir = s1_dir / mol_name if s1_dir is not None else None
        mol_metrics: Dict[str, Any] = {
            "status": mol_provenance.get("status", "unknown"),
            "smiles": mol_provenance.get("smiles", ""),
            "energy_hartree": safe_float(mol_provenance.get("energy_hartree")),
        }
        if mol_dir is not None and mol_dir.exists():
            mol_eval = evaluate_molecule(mol_dir, mol_name)
            mol_metrics.update(mol_eval)
        else:
            mol_metrics["status"] = "missing_dir"
        per_molecule[mol_name] = mol_metrics

    metrics["per_molecule"] = per_molecule

    product_mol = per_molecule.get("product", {})
    metrics["product_energy_hartree"] = product_mol.get("energy_hartree")
    metrics["best_sp_energy_hartree"] = product_mol.get("best_sp_energy_hartree")
    metrics["best_gibbs_kcal"] = product_mol.get("best_gibbs_kcal")
    metrics["energy_range_kcal"] = product_mol.get("energy_range_kcal")
    metrics["mean_energy_kcal"] = product_mol.get("mean_energy_kcal")
    metrics["final_dft_count"] = product_mol.get("final_dft_count")

    product_dir = s1_dir / "product" if s1_dir is not None else None
    if product_dir is not None:
        conformer_energies_path = product_dir / "finalDFT" / "conformer_energies.json"
        conformer_thermo_path = product_dir / "finalDFT" / "conformer_thermo.csv"
        if not conformer_energies_path.exists():
            conformer_energies_path = product_dir / "dft" / "conformer_energies.json"
        if not conformer_thermo_path.exists():
            conformer_thermo_path = product_dir / "dft" / "conformer_thermo.csv"
        metrics["paths"]["conformer_energies_json"] = str(conformer_energies_path)
        metrics["paths"]["conformer_thermo_csv"] = str(conformer_thermo_path)

    if "product" not in per_molecule:
        metrics["notes"].append("Product molecule not present in provenance.json molecules.")
    elif product_dir is None or not product_dir.exists():
        metrics["notes"].append("Missing product molecule directory.")
    else:
        product_energies_path = product_dir / "finalDFT" / "conformer_energies.json"
        product_thermo_path = product_dir / "finalDFT" / "conformer_thermo.csv"
        if not product_energies_path.exists():
            product_energies_path = product_dir / "dft" / "conformer_energies.json"
        if not product_thermo_path.exists():
            product_thermo_path = product_dir / "dft" / "conformer_thermo.csv"
        if not product_energies_path.exists():
            metrics["notes"].append("Missing product conformer_energies.json.")
        elif metrics["final_dft_count"] is None:
            metrics["notes"].append("Product conformer_energies.json missing data or unreadable.")
        if not product_thermo_path.exists():
            metrics["notes"].append("Missing product conformer_thermo.csv.")

    metrics["status"] = "ok"
    return metrics


def delta_hartree_to_kcal(other: Optional[float], anchor: Optional[float]) -> Optional[float]:
    if other is None or anchor is None:
        return None
    return (other - anchor) * HARTREE_TO_KCAL


def safe_ratio(numerator: Optional[float], denominator: Optional[float]) -> Optional[float]:
    if numerator is None or denominator in (None, 0):
        return None
    return float(numerator) / float(denominator)


def build_cross_protocol(protocol_metrics: Dict[str, Dict[str, Any]]) -> Dict[str, Dict[str, Optional[float]]]:
    anchor = protocol_metrics.get("ext")
    if not anchor or anchor.get("status") != "ok":
        return {}

    comparisons = {}
    for protocol, metrics in protocol_metrics.items():
        if protocol == "ext":
            continue
        comparisons[protocol] = {
            "delta_best_sp_energy_kcal": delta_hartree_to_kcal(
                metrics.get("best_sp_energy_hartree"),
                anchor.get("best_sp_energy_hartree"),
            ),
            "delta_best_gibbs_kcal": None
            if metrics.get("best_gibbs_kcal") is None or anchor.get("best_gibbs_kcal") is None
            else metrics.get("best_gibbs_kcal") - anchor.get("best_gibbs_kcal"),
            "nconf_ratio": safe_ratio(metrics.get("final_dft_count"), anchor.get("final_dft_count")),
            "funnel_candidate_ratio": safe_ratio(metrics.get("funnel_candidates"), anchor.get("funnel_candidates")),
        }
    return comparisons


def build_cross_protocol_per_molecule(protocol_metrics: Dict[str, Dict[str, Any]]) -> Dict[str, Dict[str, Dict[str, Optional[float]]]]:
    anchor = protocol_metrics.get("ext")
    if not anchor or anchor.get("status") != "ok":
        return {}

    anchor_mols = anchor.get("per_molecule", {})
    mol_names = anchor.get("molecule_names", [])

    comparisons: Dict[str, Dict[str, Dict[str, Optional[float]]]] = {}
    for protocol, metrics in protocol_metrics.items():
        if protocol == "ext":
            continue
        if metrics.get("status") != "ok":
            continue
        other_mols = metrics.get("per_molecule", {})
        proto_comp: Dict[str, Dict[str, Optional[float]]] = {}
        for mol_name in mol_names:
            ext_mol = anchor_mols.get(mol_name, {})
            other_mol = other_mols.get(mol_name, {})
            proto_comp[mol_name] = {
                "delta_sp_energy_kcal": delta_hartree_to_kcal(
                    other_mol.get("best_sp_energy_hartree"),
                    ext_mol.get("best_sp_energy_hartree"),
                ),
                "delta_gibbs_kcal": None
                if other_mol.get("best_gibbs_kcal") is None or ext_mol.get("best_gibbs_kcal") is None
                else other_mol.get("best_gibbs_kcal") - ext_mol.get("best_gibbs_kcal"),
                "nconf_ratio": safe_ratio(other_mol.get("final_dft_count"), ext_mol.get("final_dft_count")),
            }
        comparisons[protocol] = proto_comp
    return comparisons


def format_bool(value: Any) -> str:
    if value is True:
        return "yes"
    if value is False:
        return "no"
    return "—"


def format_status(value: str) -> str:
    return "✅" if value == "ok" else "❌"


def format_int(value: Any) -> str:
    return "—" if value is None else str(value)


def format_float(value: Any, decimals: int) -> str:
    if value is None:
        return "—"
    return "{0:.{1}f}".format(value, decimals)


def format_signed_float(value: Any, decimals: int) -> str:
    if value is None:
        return "—"
    return "{0:+.{1}f}".format(value, decimals)


def protocol_label(value: Optional[str]) -> str:
    return value if value not in (None, "") else "—"


def render_protocol_table(protocols: List[str], protocol_metrics: Dict[str, Dict[str, Any]]) -> List[str]:
    rows = [
        ("Status", lambda m: format_status(m.get("status"))),
        ("Two-stage", lambda m: format_bool(m.get("two_stage_enabled"))),
        ("Freq", lambda m: format_bool(m.get("freq_requested"))),
        ("Final SP", lambda m: format_bool(m.get("final_sp_requested"))),
        ("Handoff mode", lambda m: protocol_label(m.get("handoff_mode_effective") or m.get("handoff_mode_requested"))),
        ("Fallback", lambda m: format_bool(m.get("fallback_triggered"))),
        ("Funnel candidates", lambda m: format_int(m.get("funnel_candidates"))),
        ("Funnel survivors", lambda m: format_int(m.get("funnel_survivors"))),
        ("DFT optimized", lambda m: format_int(m.get("final_dft_count"))),
        ("Best E (Ha)", lambda m: format_float(m.get("best_sp_energy_hartree"), 6)),
        ("Best G (kcal/mol)", lambda m: format_float(m.get("best_gibbs_kcal"), 3)),
        ("E range (kcal)", lambda m: format_float(m.get("energy_range_kcal"), 3)),
    ]

    lines = [
        "| Metric | {0} |".format(" | ".join(protocols)),
        "|--------|{0}|".format("|".join("-----" for _ in protocols)),
    ]
    for label, formatter in rows:
        values = [formatter(protocol_metrics[protocol]) for protocol in protocols]
        lines.append("| {0} | {1} |".format(label, " | ".join(values)))
    return lines


def render_cross_protocol_table(protocols: List[str], cross_protocol: Dict[str, Dict[str, Optional[float]]]) -> List[str]:
    others = [protocol for protocol in protocols if protocol != "ext" and protocol in cross_protocol]
    if not others:
        return []

    rows = [
        ("ΔE SP (kcal/mol)", lambda m: format_signed_float(m.get("delta_best_sp_energy_kcal"), 3)),
        ("ΔG (kcal/mol)", lambda m: format_signed_float(m.get("delta_best_gibbs_kcal"), 3)),
        ("DFT count ratio", lambda m: format_float(m.get("nconf_ratio"), 3)),
        ("Funnel ratio", lambda m: format_float(m.get("funnel_candidate_ratio"), 3)),
    ]

    lines = [
        "| Metric | {0} |".format(" | ".join(others)),
        "|--------|{0}|".format("|".join("------" for _ in others)),
    ]
    for label, formatter in rows:
        values = [formatter(cross_protocol[protocol]) for protocol in others]
        lines.append("| {0} | {1} |".format(label, " | ".join(values)))
    return lines


def render_per_molecule_tables(protocols: List[str], protocol_metrics: Dict[str, Dict[str, Any]]) -> List[str]:
    mol_names = protocol_metrics.get("ext", {}).get("molecule_names", [])
    if not mol_names:
        return []

    lines = []
    for mol_name in mol_names:
        lines.append("### {0}".format(mol_name))
        mol_rows = [
            ("Status", lambda m: format_status(m.get("per_molecule", {}).get(mol_name, {}).get("status", "failed"))),
            ("DFT optimized", lambda m: format_int(m.get("per_molecule", {}).get(mol_name, {}).get("final_dft_count"))),
            ("Best E (Ha)", lambda m: format_float(m.get("per_molecule", {}).get(mol_name, {}).get("best_sp_energy_hartree"), 6)),
            ("Best G (kcal/mol)", lambda m: format_float(m.get("per_molecule", {}).get(mol_name, {}).get("best_gibbs_kcal"), 3)),
            ("E range (kcal)", lambda m: format_float(m.get("per_molecule", {}).get(mol_name, {}).get("energy_range_kcal"), 3)),
            ("Small molecule", lambda m: format_bool(m.get("per_molecule", {}).get(mol_name, {}).get("is_small_molecule"))),
        ]
        lines.append("| Metric | {0} |".format(" | ".join(protocols)))
        lines.append("|--------|{0}|".format("|".join("------" for _ in protocols)))
        for label, formatter in mol_rows:
            values = [formatter(protocol_metrics[protocol]) for protocol in protocols]
            lines.append("| {0} | {1} |".format(label, " | ".join(values)))
        lines.append("")
    return lines


def render_cross_protocol_per_molecule(protocols: List[str], cross_per_mol: Dict[str, Dict[str, Dict[str, Optional[float]]]]) -> List[str]:
    others = [p for p in protocols if p != "ext" and p in cross_per_mol]
    if not others:
        return []

    mol_names = list(cross_per_mol.get(others[0], {}).keys()) if others else []
    if not mol_names:
        return []

    lines = []
    for mol_name in mol_names:
        lines.append("### {0}".format(mol_name))
        rows = [
            ("ΔE SP (kcal/mol)", lambda m: format_signed_float(m.get("delta_sp_energy_kcal"), 3)),
            ("ΔG (kcal/mol)", lambda m: format_signed_float(m.get("delta_gibbs_kcal"), 3)),
            ("DFT count ratio", lambda m: format_float(m.get("nconf_ratio"), 3)),
        ]
        lines.append("| Metric | {0} |".format(" | ".join(others)))
        lines.append("|--------|{0}|".format("|".join("------" for _ in others)))
        for label, formatter in rows:
            values = [formatter(cross_per_mol[protocol].get(mol_name, {})) for protocol in others]
            lines.append("| {0} | {1} |".format(label, " | ".join(values)))
        lines.append("")
    return lines


def render_behavior_summary(protocols: List[str], protocol_metrics: Dict[str, Dict[str, Any]]) -> List[str]:
    lines = []
    for protocol in protocols:
        metrics = protocol_metrics[protocol]
        lines.append("### {0}".format(protocol))
        lines.append("- Funnel candidates: {0}".format(format_int(metrics.get("funnel_candidates"))))
        lines.append("- Handoff: {0} (effective: {1})".format(
            protocol_label(metrics.get("handoff_mode_requested")),
            protocol_label(metrics.get("handoff_mode_effective")),
        ))
        lines.append("- Fallback: {0}".format("triggered" if metrics.get("fallback_triggered") is True else "not triggered" if metrics.get("fallback_triggered") is False else "unknown"))
        lines.append("- Selection: {0}".format(protocol_label(metrics.get("selection_mode"))))
        if metrics.get("notes"):
            lines.append("- Notes: {0}".format("; ".join(metrics["notes"])))
        lines.append("")
    return lines


def build_markdown_report(
    rx_id: str,
    protocols: List[str],
    protocol_metrics: Dict[str, Dict[str, Any]],
    cross_protocol: Dict[str, Dict[str, Optional[float]]],
    cross_per_mol: Dict[str, Dict[str, Dict[str, Optional[float]]]],
    generated_at: str,
) -> str:
    succeeded = sum(1 for protocol in protocols if protocol_metrics[protocol].get("status") == "ok")
    lines = [
        "# Benchmark Report: rx{0}".format(rx_id),
        "",
        "Generated: {0}".format(generated_at),
        "",
        "## Overall: {0}/{1} protocols succeeded".format(succeeded, len(protocols)),
        "",
        "## Per-Protocol Summary",
        "",
    ]
    lines.extend(render_protocol_table(protocols, protocol_metrics))
    lines.extend(["", "## Per-Molecule Detail", ""])
    per_molecule_tables = render_per_molecule_tables(protocols, protocol_metrics)
    if per_molecule_tables:
        lines.extend(per_molecule_tables)
    else:
        lines.append("Per-molecule detail unavailable because ext did not provide molecule metadata.")
    lines.extend(["", "## Cross-Protocol Comparison (vs ext)", ""])
    if cross_protocol:
        lines.extend(render_cross_protocol_table(protocols, cross_protocol))
    else:
        lines.append("Cross-protocol comparison skipped because ext did not succeed or was not requested.")
    lines.extend(["", "## Cross-Protocol Per-Molecule (vs ext)", ""])
    if cross_per_mol:
        lines.extend(render_cross_protocol_per_molecule(protocols, cross_per_mol))
    else:
        lines.append("Per-molecule cross-protocol comparison skipped because ext did not succeed or lacked comparable molecule data.")
    lines.extend(["", "## Protocol Behavior Summary", ""])
    lines.extend(render_behavior_summary(protocols, protocol_metrics))
    return "\n".join(lines).rstrip() + "\n"


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def evaluate_rx(output_base: Path, rx_id: str, protocols: List[str]) -> Dict[str, Any]:
    rx_dir = output_base / "rx{0}".format(rx_id)
    protocol_metrics = {}
    for protocol in protocols:
        protocol_metrics[protocol] = evaluate_protocol(rx_dir, protocol)
    cross_protocol = build_cross_protocol(protocol_metrics)
    cross_per_mol = build_cross_protocol_per_molecule(protocol_metrics)
    generated_at = datetime.now().isoformat(timespec="seconds")

    summary = {
        "rx_id": str(rx_id),
        "rx_dir": str(rx_dir),
        "generated_at": generated_at,
        "protocols_requested": list(protocols),
        "protocols": protocol_metrics,
        "cross_protocol_vs_ext": cross_protocol,
        "cross_protocol_per_molecule_vs_ext": cross_per_mol,
    }

    evaluation_dir = rx_dir / "evaluation"
    summary_json = json.dumps(summary, indent=2, sort_keys=True)
    summary_md = build_markdown_report(str(rx_id), protocols, protocol_metrics, cross_protocol, cross_per_mol, generated_at)

    write_text(evaluation_dir / "summary.json", summary_json + "\n")
    write_text(evaluation_dir / "summary.md", summary_md)
    return summary


def main() -> int:
    args = parse_args()
    output_base = Path(args.output_base).expanduser()
    protocols = [str(protocol) for protocol in args.protocols]
    rx_ids = [str(rx_id) for rx_id in args.rx_ids]

    for rx_id in rx_ids:
        evaluate_rx(output_base, rx_id, protocols)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
