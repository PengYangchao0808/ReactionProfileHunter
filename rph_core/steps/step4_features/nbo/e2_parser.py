"""
NBO E(2) parser for ORCA output.
"""

import re
from typing import List, Dict, Any, Optional, Tuple


def normalize_nbo_label(label: str) -> str:
    if not label:
        return ""
    return re.sub(r"\s+", " ", label).strip()


def find_nbo_e2_section(text: str) -> Optional[str]:
    anchors = [
        r"SECOND ORDER PERTURBATION THEORY ANALYSIS",
        r"E\(2\) PERTURBATION THEORY ANALYSIS",
    ]
    for pattern in anchors:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            start_idx = match.end()
            remainder = text[start_idx:]

            stop_pattern = r"\n\s*(?:NATURAL\s+(?:POPULATIONS|BOND|HYBRID)|NHO\s+DIRECTIONAL|SUMMARY\s+OF|\*{10,})"
            end_match = re.search(stop_pattern, remainder)

            if end_match:
                return remainder[:end_match.start()]
            else:
                return remainder[:50000]
    return None


def parse_nbo_e2_table(section: str) -> List[Dict[str, Any]]:
    interactions = []
    pattern = re.compile(
        r"^\s*\d+\.\s+(?P<donor>.+?)\s*(?:/|->)\s*(?P<acceptor>.+?)\s+(?P<e2>\d+\.\d+)",
        re.MULTILINE
    )

    for match in pattern.finditer(section):
        try:
            d_label = match.group("donor").strip()
            a_label = match.group("acceptor").strip()
            e2_val = float(match.group("e2"))

            interactions.append({
                "donor": d_label,
                "acceptor": a_label,
                "e2": e2_val
            })
        except ValueError:
            continue

    return interactions


def match_templates(
    interactions: List[Dict[str, Any]],
    templates: Dict[str, Dict[str, Any]]
) -> Tuple[Dict[str, float], List[str]]:
    features = {}
    warnings = []

    norm_interactions = []
    for row in interactions:
        norm_interactions.append({
            "donor": normalize_nbo_label(row["donor"]),
            "acceptor": normalize_nbo_label(row["acceptor"]),
            "e2": row["e2"],
            "raw": row
        })

    matched_any = False

    for tmpl_name, tmpl_def in templates.items():
        t_donor = normalize_nbo_label(tmpl_def.get("donor", ""))
        t_acceptor = normalize_nbo_label(tmpl_def.get("acceptor", ""))

        matches = [
            i for i in norm_interactions
            if i["donor"] == t_donor and i["acceptor"] == t_acceptor
        ]

        if len(matches) == 1:
            features[f"nbo.e2.{tmpl_name}.e2_kcal"] = matches[0]["e2"]
            matched_any = True
        elif len(matches) > 1:
            warnings.append(f"Template '{tmpl_name}' matched {len(matches)} interactions (Ambiguous). Ignored.")

    if not matched_any and templates:
        warnings.append("No templates matched any extracted interactions.")

    return features, warnings
