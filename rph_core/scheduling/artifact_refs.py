from __future__ import annotations

import json
import logging
import platform
import re
import shutil
from pathlib import Path
from typing import Any, Iterable, Optional

logger = logging.getLogger(__name__)
_SAFE_SEGMENT_RE = re.compile(r"^[A-Za-z0-9_-]+$")


def _can_symlink() -> bool:
    return platform.system() != "Windows"


def _safe_segment(value: str, label: str) -> str:
    text = str(value or "").strip()
    if not text or not _SAFE_SEGMENT_RE.match(text):
        raise ValueError(f"Unsafe {label}: {value!r}")
    return text


def _is_relative_to(path: Path, root: Path) -> bool:
    try:
        path.resolve(strict=False).relative_to(root.resolve(strict=False))
        return True
    except ValueError:
        return False


def link_or_copy_precursor_into_branch(reaction_root: Path, branch_root: Path) -> None:
    precursor_s1_dir = reaction_root / "precursor" / "S1_ConfGeneration" / "precursor"
    if not precursor_s1_dir.exists():
        return
    if not _is_relative_to(precursor_s1_dir, reaction_root):
        raise ValueError(f"Precursor source escapes reaction root: {precursor_s1_dir}")
    if not _is_relative_to(branch_root, reaction_root / "branches"):
        raise ValueError(f"Branch root escapes branches directory: {branch_root}")

    target_dir = branch_root / "S1_ConfGeneration"
    target_dir.mkdir(parents=True, exist_ok=True)
    link_path = target_dir / "precursor"

    if link_path.exists():
        return
    if link_path.is_symlink():
        link_path.unlink()

    try:
        if _can_symlink():
            link_path.symlink_to(precursor_s1_dir.resolve(strict=False), target_is_directory=True)
        else:
            shutil.copytree(precursor_s1_dir, link_path, dirs_exist_ok=True)
        logger.debug(f"Linked precursor into branch: {link_path} -> {precursor_s1_dir}")
    except Exception as exc:
        logger.warning(f"Failed to link/copy precursor into branch: {exc}. Falling back to copy.")
        try:
            if link_path.exists() or link_path.is_symlink():
                if link_path.is_symlink():
                    link_path.unlink()
                else:
                    shutil.rmtree(link_path, ignore_errors=True)
            shutil.copytree(precursor_s1_dir, link_path, dirs_exist_ok=True)
        except Exception as exc2:
            logger.error(f"Copy fallback also failed: {exc2}")


def materialize_small_molecule_refs(
    branch_root: Path,
    cache_root: Path,
    keys: Iterable[str],
    config: Optional[dict[str, Any]] = None,
) -> None:
    from rph_core.utils.small_molecule_cache import SmallMoleculeCache
    from rph_core.utils.small_molecule_catalog import SmallMoleculeCatalog

    target_dir = branch_root / "S1_ConfGeneration" / "small_molecules"
    target_dir.mkdir(parents=True, exist_ok=True)
    cache = SmallMoleculeCache(cache_root)
    catalog = SmallMoleculeCatalog(config or {}) if config is not None else None
    manifest: dict[str, Any] = {"refs": {}, "warnings": []}

    for key in keys:
        key = str(key).strip()
        if not key:
            continue
        try:
            key = _safe_segment(key, "small_molecule_key")
        except ValueError as exc:
            manifest["warnings"].append(str(exc))
            continue

        cache_entry = None
        if catalog is not None:
            molecule = catalog.get(key)
            if molecule is not None:
                cache_entry = cache.get_path(molecule.smiles)
        if cache_entry is None:
            candidate = cache_root / key
            if candidate.exists():
                cache_entry = candidate
        if cache_entry is None or not cache_entry.exists():
            manifest["warnings"].append(f"W_MISSING_SMALL_MOLECULE_CACHE:{key}")
            continue
        if not _is_relative_to(cache_entry, cache_root):
            manifest["warnings"].append(f"W_SMALL_MOLECULE_CACHE_ESCAPE:{key}")
            continue

        ref_path = target_dir / key
        if not _is_relative_to(ref_path, target_dir):
            manifest["warnings"].append(f"W_SMALL_MOLECULE_REF_ESCAPE:{key}")
            continue
        try:
            if ref_path.exists() or ref_path.is_symlink():
                continue
            if _can_symlink():
                ref_path.symlink_to(cache_entry, target_is_directory=True)
            else:
                shutil.copytree(cache_entry, ref_path, dirs_exist_ok=True)
            manifest["refs"][key] = str(cache_entry)
        except Exception as exc:
            manifest["warnings"].append(f"W_SMALL_MOLECULE_REF_FAILED:{key}:{exc}")

    with open(target_dir / "small_molecule_refs.json", "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)


def link_or_copy_s0_into_branch(reaction_root: Path, branch_root: Path) -> None:
    """Materialize reaction-level S0 planning artifacts into a branch workspace.

    The v3 scheduler runs S0 once at the reaction root during planning (Stage 0),
    writing mechanism_summary.json, mechanism_graph.json, atom_map_smiles.json,
    and dr_branch_plan.json into reaction_root/S0_Mechanism/.  Each branch
    skips S0 but still needs these artifacts for S2 forming_bonds resolution.

    This function copies the relevant JSON artifacts from the reaction-level
    S0 directory into branch_root/S0_Mechanism/, following the same pattern
    as link_or_copy_precursor_into_branch().
    """
    src = reaction_root / "S0_Mechanism"
    dst = branch_root / "S0_Mechanism"

    if not src.exists():
        logger.debug("No reaction-level S0_Mechanism directory to propagate")
        return

    dst.mkdir(parents=True, exist_ok=True)

    _S0_ARTIFACTS = (
        "mechanism_summary.json",
        "mechanism_graph.json",
        "atom_map_smiles.json",
        "dr_branch_plan.json",
    )

    for artifact in _S0_ARTIFACTS:
        src_file = src / artifact
        if not src_file.exists():
            continue

        dst_file = dst / artifact
        if dst_file.exists():
            continue

        try:
            shutil.copy2(src_file, dst_file)
            logger.debug(f"Copied S0 artifact into branch: {artifact}")
        except Exception as exc:
            logger.warning(f"Failed to copy S0 artifact {artifact} into branch: {exc}")


def read_branch_manifest(branch_root: Path) -> dict[str, Any]:
    manifest_path = branch_root / "branch_manifest.json"
    if not manifest_path.exists():
        return {}
    try:
        with open(manifest_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def write_branch_manifest(branch_root: Path, payload: dict[str, Any]) -> Path:
    branch_root.mkdir(parents=True, exist_ok=True)
    manifest_path = branch_root / "branch_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)
    return manifest_path
