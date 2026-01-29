"""
NBO E(2) Template Loader and Matcher (V4.2 Phase B/C)
======================================================

Load NBO E(2) template whitelist from YAML and match interactions.
Phase A/B: Interface only (dry-run, no matching).
Phase C: Template matching applied to parsed NBO outputs.

Author: QC Descriptors Team
Date: 2026-01-18
"""

import yaml
from pathlib import Path
from typing import Dict, List, Any, Optional
from dataclasses import dataclass


@dataclass
class NBOE2Template:
    """NBO E(2) analysis template specification."""
    name: str
    donor: str
    acceptor: str
    required_atoms: Optional[List[int]] = None
    spin: Optional[int] = None
    notes: str = ""


def load_nbo_templates(path: Path) -> List[NBOE2Template]:
    """Load NBO E(2) templates from YAML file.

    Args:
        path: Path to nbo_templates.yaml

    Returns:
        List of NBOE2Template objects
    """
    if not path.exists():
        return []

    with open(path, "r") as f:
        data = yaml.safe_load(f)

    if not data or "nbo_e2_templates" not in data:
        return []

    templates = []
    for spec in data["nbo_e2_templates"]:
        if not isinstance(spec, dict):
            continue

        name = spec.get("name", "")
        donor = spec.get("donor", "")
        acceptor = spec.get("acceptor", "")
        required_atoms = spec.get("required_atoms", None)
        spin = spec.get("spin", None)
        notes = spec.get("notes", "")

        if not name or not donor or not acceptor:
            continue

        templates.append(NBOE2Template(
            name=name,
            donor=donor,
            acceptor=acceptor,
            required_atoms=required_atoms,
            spin=spin,
            notes=notes,
        ))

    return templates


def match_e2_interactions(
    parsed_e2_list: List[Dict[str, Any]],
    templates: List[NBOE2Template],
) -> Dict[str, Optional[float]]:
    """Match parsed E(2) interactions against template whitelist (Contract 3).

    Args:
        parsed_e2_list: List of parsed E(2) interactions from NBO output
        templates: List of NBOE2Template whitelist

    Returns:
        Dict mapping template.name -> E(2) energy or None
        None if 0 or >1 matches found, or template mismatch
    """
    results = {}

    for template in templates:
        matches = []
        for interaction in parsed_e2_list:
            donor_atom = interaction.get("donor_atom", "")
            acceptor_atom = interaction.get("acceptor_atom", "")
            e2_energy = interaction.get("e2_energy")

            if not isinstance(donor_atom, str) or not isinstance(acceptor_atom, str):
                continue

            donor_match = donor_atom.strip().lower() == template.donor.lower()
            acceptor_match = acceptor_atom.strip().lower() == template.acceptor.lower()

            if donor_match and acceptor_match:
                matches.append(e2_energy)

        if len(matches) == 1:
            results[template.name] = matches[0]
        elif len(matches) > 1:
            results[template.name] = None
        else:
            results[template.name] = None

    return results
