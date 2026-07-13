from __future__ import annotations

from typing import Dict, List

from rph_core.utils.method_registry import NormalizedMethodSpec


class OrcaInputRenderer:
    def __init__(self, spec: NormalizedMethodSpec):
        self.spec = spec

    def render_simple_keywords(self, task_type: str = "sp") -> str:
        parts: List[str] = ["!", self.spec.method]

        if self.spec.family != "composite_3c" and self.spec.basis:
            parts.append(self.spec.basis)

        accel = self.spec.scf.accel.strip()
        aux_j = self.spec.scf.aux_j
        if aux_j and self.spec.family != "composite_3c":
            parts.append(aux_j)

        if accel and accel.lower() != "none":
            parts.append(accel)

        aux_c = self.spec.correlation.aux_c
        if aux_c and self.spec.family != "composite_3c":
            parts.append(aux_c)

        task = str(task_type or "sp").strip().lower()
        if task == "opt":
            parts.append("Opt")
        elif task == "opt_freq":
            parts.extend(["Opt", "Freq"])
        elif task == "ts":
            parts.append("OptTS")
        elif task == "ts_freq":
            parts.extend(["OptTS", "Freq"])

        if self.spec.family != "composite_3c":
            parts.append("tightSCF")

        if self.spec.pno.preset:
            parts.append(self.spec.pno.preset)

        if self.spec.dispersion.mode in {"external", "required"} and self.spec.dispersion.keyword:
            parts.append(self.spec.dispersion.keyword)

        if self.spec.route_extras:
            # ORCA simple-input keywords are case-insensitive and a duplicate
            # (notably ``TightSCF``) is a fatal input error in ORCA 5.  The
            # renderer owns canonical keywords, so only append genuinely new
            # extras supplied by configuration.
            known_keywords = {keyword.lower() for keyword in parts}
            extra_keywords = [
                keyword
                for keyword in self.spec.route_extras.split()
                if keyword.lower() not in known_keywords
            ]
            parts.extend(extra_keywords)

        parts.extend(["noautostart", "miniprint", "nopop"])
        return " ".join(parts)

    def render_blocks(self) -> Dict[str, str]:
        # ORCA accepts auxiliary basis names such as ``def2/J`` and
        # ``def2-TZVPP/C`` as simple keywords.  They are deliberately not
        # rendered into %basis or %method blocks: ORCA 5 rejects ``AuxJ`` in
        # %method and %basis only accepts NewGTO/NewECP-style definitions.
        return {}
