"""
Forming bonds resolver
======================

Infer and load forming bond indices from optimized geometries (S1/S3).
Provides a small metadata contract (mechanism_meta.json) for S4 reuse.

Author: QCcalc Team
Date: 2026-02-04
"""

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from rph_core.utils.file_io import read_xyz
from rph_core.utils.geometry_tools import GeometryUtils
from rph_core.utils.molecular_graph import build_bond_graph

logger = logging.getLogger(__name__)

MECHANISM_META_FILENAME = "mechanism_meta.json"


@dataclass
class FormingBondsResult:
    """Result of forming bond resolution.

    Attributes:
        forming_bonds: Tuple of forming bond pairs or None
        meta: Mechanism metadata dict if available/generated
        warnings: List of warning strings
    """

    forming_bonds: Optional[Tuple[Tuple[int, int], ...]]
    meta: Optional[Dict[str, Any]]
    warnings: List[str]


def load_mechanism_meta(meta_path: Path) -> Optional[Dict[str, Any]]:
    """Load mechanism_meta.json from disk.

    Args:
        meta_path: Path to mechanism_meta.json

    Returns:
        Parsed dict or None if load fails
    """

    if not meta_path.exists():
        return None

    try:
        return json.loads(meta_path.read_text(encoding="utf-8"))
    except Exception as e:
        logger.warning(f"Failed to read {meta_path}: {e}")
        return None


def _normalize_pairs(pairs: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    normalized: List[Tuple[int, int]] = []
    for i, j in pairs:
        if i == j:
            continue
        a, b = (i, j) if i < j else (j, i)
        normalized.append((a, b))
    return sorted(set(normalized))


def parse_forming_bonds(meta: Dict[str, Any]) -> Optional[Tuple[Tuple[int, int], ...]]:
    """Parse forming bonds from mechanism meta.

    Args:
        meta: Parsed mechanism_meta.json content

    Returns:
        Tuple of forming bond pairs or None if invalid/missing
    """

    raw = meta.get("forming_bonds")
    if not raw:
        return None

    index_base = int(meta.get("index_base", 0))
    pairs: List[Tuple[int, int]] = []

    try:
        for pair in raw:
            if not isinstance(pair, (list, tuple)) or len(pair) != 2:
                return None
            i, j = int(pair[0]), int(pair[1])
            if index_base == 1:
                i -= 1
                j -= 1
            if i < 0 or j < 0:
                return None
            pairs.append((i, j))
    except Exception:
        return None

    normalized = _normalize_pairs(pairs)
    if not normalized:
        return None
    return tuple(normalized)


def _graph_to_bond_set(graph: Dict[int, List[int]]) -> set[Tuple[int, int]]:
    bonds: set[Tuple[int, int]] = set()
    for i, neighbors in graph.items():
        for j in neighbors:
            if i < j:
                bonds.add((i, j))
    return bonds


def infer_forming_bonds_from_geometries(
    product_xyz: Path,
    ts_xyz: Path,
    config: Optional[Dict[str, Any]] = None
) -> FormingBondsResult:
    """Infer forming bonds from optimized product and TS geometries.

    Args:
        product_xyz: Path to optimized product XYZ (S1)
        ts_xyz: Path to optimized TS XYZ (S3)
        config: Optional configuration dict for thresholds

    Returns:
        FormingBondsResult with inferred forming bonds and metadata
    """

    cfg = config or {}
    product_scale = float(cfg.get("product_graph_scale", 1.25))
    ts_scale = float(cfg.get("ts_graph_scale", 1.15))
    dist_min = float(cfg.get("forming_distance_min", 1.5))
    dist_max = float(cfg.get("forming_distance_max", 3.5))
    max_pairs = int(cfg.get("max_pairs", 2))

    warnings: List[str] = []

    product_coords, product_symbols = read_xyz(product_xyz)
    ts_coords, ts_symbols = read_xyz(ts_xyz)

    if product_coords is None or ts_coords is None:
        warnings.append("W_FORMING_BONDS_READ_FAILED")
        return FormingBondsResult(None, None, warnings)

    if len(product_symbols) != len(ts_symbols):
        warnings.append("W_FORMING_BONDS_ATOM_COUNT_MISMATCH")
        return FormingBondsResult(None, None, warnings)

    if product_symbols != ts_symbols:
        warnings.append("W_FORMING_BONDS_SYMBOL_MISMATCH")

    product_graph = build_bond_graph(product_coords, product_symbols, scale=product_scale)
    ts_graph = build_bond_graph(ts_coords, ts_symbols, scale=ts_scale)

    product_bonds = _graph_to_bond_set(product_graph)
    ts_bonds = _graph_to_bond_set(ts_graph)

    candidates = product_bonds - ts_bonds
    if not candidates:
        warnings.append("W_FORMING_BONDS_NO_CANDIDATES")
        return FormingBondsResult(None, None, warnings)

    distances: List[Tuple[float, Tuple[int, int]]] = []
    for i, j in candidates:
        d_ij = GeometryUtils.calculate_distance(ts_coords, i, j)
        if dist_min <= d_ij <= dist_max:
            distances.append((float(d_ij), (i, j)))

    if not distances:
        warnings.append("W_FORMING_BONDS_DISTANCE_FILTER_EMPTY")
        return FormingBondsResult(None, None, warnings)

    distances.sort(key=lambda x: x[0])
    selected = [pair for _, pair in distances[:max_pairs]]
    forming_bonds = tuple(_normalize_pairs(selected)) if selected else None

    validation_status = "pass" if forming_bonds else "fail"
    if forming_bonds and len(forming_bonds) < max_pairs:
        validation_status = "warn"
        warnings.append("W_FORMING_BONDS_FEWER_THAN_EXPECTED")

    meta = {
        "version": "1",
        "source": {
            "product_xyz": str(product_xyz),
            "ts_xyz": str(ts_xyz),
            "derived_from_steps": ["S1", "S3"],
        },
        "index_base": 0,
        "forming_bonds": [list(b) for b in forming_bonds] if forming_bonds else None,
        "method": {
            "name": "graph_diff_product_vs_ts_strict",
            "parameters": {
                "product_graph_scale": product_scale,
                "ts_graph_scale": ts_scale,
                "forming_distance_min": dist_min,
                "forming_distance_max": dist_max,
                "max_pairs": max_pairs,
            },
        },
        "validation": {
            "status": validation_status,
            "warnings": warnings,
        },
    }

    return FormingBondsResult(forming_bonds, meta, warnings)


def write_mechanism_meta(meta: Dict[str, Any], meta_path: Path) -> None:
    """Write mechanism meta JSON to disk.

    Args:
        meta: Metadata dict
        meta_path: Target path
    """

    meta_path.parent.mkdir(parents=True, exist_ok=True)
    meta_path.write_text(json.dumps(meta, indent=2, ensure_ascii=False))


def resolve_forming_bonds(
    product_xyz: Optional[Path],
    ts_xyz: Optional[Path],
    s3_dir: Optional[Path],
    s4_dir: Optional[Path],
    config: Optional[Dict[str, Any]] = None,
    write_meta: bool = True
) -> FormingBondsResult:
    """Resolve forming bonds from mechanism_meta.json or infer from geometries.

    Args:
        product_xyz: Product XYZ path
        ts_xyz: TS XYZ path
        s3_dir: Step3 directory (preferred meta location)
        s4_dir: Step4 directory (fallback meta location)
        config: Optional config dict for inference thresholds
        write_meta: Whether to write inferred mechanism_meta.json

    Returns:
        FormingBondsResult with forming bonds and metadata
    """

    meta_paths: List[Path] = []
    if s3_dir:
        meta_paths.append(Path(s3_dir) / MECHANISM_META_FILENAME)
    if s4_dir:
        meta_paths.append(Path(s4_dir) / MECHANISM_META_FILENAME)

    for meta_path in meta_paths:
        meta = load_mechanism_meta(meta_path)
        if meta:
            forming_bonds = parse_forming_bonds(meta)
            if forming_bonds:
                return FormingBondsResult(forming_bonds, meta, [])

    if product_xyz is None or ts_xyz is None:
        return FormingBondsResult(None, None, ["W_FORMING_BONDS_MISSING_INPUTS"])

    result = infer_forming_bonds_from_geometries(product_xyz, ts_xyz, config=config)
    if write_meta and result.meta:
        target_dir = Path(s3_dir) if s3_dir else Path(s4_dir) if s4_dir else None
        if target_dir is not None:
            write_mechanism_meta(result.meta, target_dir / MECHANISM_META_FILENAME)

    return result
