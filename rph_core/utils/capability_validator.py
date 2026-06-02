from __future__ import annotations

from typing import List

from rph_core.utils.method_registry import NormalizedMethodSpec


class CapabilityValidator:
    @staticmethod
    def validate(spec: NormalizedMethodSpec, task_type: str) -> List[str]:
        errors: List[str] = []
        profile = spec.family
        task = str(task_type or "").strip().lower()

        if profile == "composite_3c":
            if spec.aux_basis:
                errors.append("Composite 3c methods must not define aux_basis")
            if spec.dispersion.mode not in {"forbidden", "none", "auto_builtin"}:
                errors.append("Composite 3c methods forbid external dispersion")

        if profile == "vv10_family" and spec.dispersion.keyword:
            errors.append("VV10-family methods must not use external D3/D4 dispersion")

        if profile == "double_hybrid_dft":
            if spec.correlation.model.lower() == "none":
                errors.append("Double-hybrid methods require non-none correlation model")
            if not spec.correlation.aux_c:
                errors.append("Double-hybrid methods require correlation aux basis (aux_c)")

        if task in {"opt", "opt_freq", "ts", "ts_freq", "geo_benchmark"} and spec.scf.accel.upper() == "RIJK":
            errors.append("RIJK accelerator is unsupported for ORCA gradient tasks")

        if task in {"opt", "opt_freq", "ts", "ts_freq", "geo_benchmark"} and profile == "post_hf_dlpno":
            errors.append("DLPNO/post-HF methods are not supported for geometry gradient tasks")

        if spec.dispersion.mode in {"required", "external"} and not spec.dispersion.keyword:
            errors.append("Dispersion mode requires an explicit keyword")

        return errors
