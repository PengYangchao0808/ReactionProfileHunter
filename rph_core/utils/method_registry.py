from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Mapping, Optional

from rph_core.utils.keyword_translator import KeywordTranslator


@dataclass(frozen=True)
class SCFSpec:
    accel: str = "RIJCOSX"
    aux_j: Optional[str] = "def2/J"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "accel": self.accel,
            "aux_j": self.aux_j,
        }


@dataclass(frozen=True)
class CorrelationSpec:
    model: str = "none"
    aux_c: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "model": self.model,
            "aux_c": self.aux_c,
        }


@dataclass(frozen=True)
class DispersionSpec:
    mode: str = "none"
    keyword: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "mode": self.mode,
            "keyword": self.keyword,
        }


@dataclass(frozen=True)
class PNOSpec:
    preset: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "preset": self.preset,
        }


@dataclass(frozen=True)
class MethodProfile:
    family: str
    base_method: str
    scf_default: SCFSpec
    correlation_default: CorrelationSpec
    dispersion_default: DispersionSpec
    self_contained: bool = False
    gradients_supported: bool = True
    forbidden_ri_for_gradients: tuple[str, ...] = ("RIJK",)
    pno_default: PNOSpec = PNOSpec()


@dataclass(frozen=True)
class NormalizedMethodSpec:
    engine: str
    family: str
    method: str
    basis: str
    scf: SCFSpec
    correlation: CorrelationSpec
    dispersion: DispersionSpec
    pno: PNOSpec
    route_extras: str = ""
    solvent: Optional[str] = None
    maxcore: Optional[int] = None
    aux_basis: Optional[str] = None
    is_reference: bool = False
    reuse_baseline_gcorr: bool = True
    raw: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        payload: Dict[str, Any] = {
            "engine": self.engine,
            "family": self.family,
            "method": self.method,
            "basis": self.basis,
            "scf": self.scf.to_dict(),
            "correlation": self.correlation.to_dict(),
            "dispersion": self.dispersion.to_dict(),
            "pno": self.pno.to_dict(),
            "route_extras": self.route_extras,
            "solvent": self.solvent,
            "maxcore": self.maxcore,
            "aux_basis": self.aux_basis,
            "is_reference": self.is_reference,
            "reuse_baseline_gcorr": self.reuse_baseline_gcorr,
        }
        if self.raw is not None:
            payload["raw"] = self.raw
        return payload


def _to_scf_spec(raw: Mapping[str, Any], default: SCFSpec) -> SCFSpec:
    accel = str(raw.get("accel", default.accel) or default.accel).strip() or default.accel
    aux_j_raw = raw.get("aux_j", default.aux_j)
    aux_j = None if aux_j_raw is None else str(aux_j_raw).strip() or None
    if aux_j == "auto":
        aux_j = default.aux_j
    return SCFSpec(accel=accel, aux_j=aux_j)


def _to_correlation_spec(raw: Mapping[str, Any], default: CorrelationSpec, basis: str) -> CorrelationSpec:
    model = str(raw.get("model", default.model) or default.model).strip() or "none"
    aux_c_raw = raw.get("aux_c", default.aux_c)
    aux_c = None if aux_c_raw is None else str(aux_c_raw).strip() or None
    if aux_c == "auto":
        aux_c = default.aux_c or f"{basis}/C"
    if model.lower() != "none" and not aux_c:
        aux_c = f"{basis}/C"
    return CorrelationSpec(model=model, aux_c=aux_c)


def _to_dispersion_spec(raw: Mapping[str, Any], default: DispersionSpec) -> DispersionSpec:
    mode = str(raw.get("mode", default.mode) or default.mode).strip() or "none"
    keyword_raw = raw.get("keyword", default.keyword)
    keyword = None if keyword_raw is None else str(keyword_raw).strip() or None
    return DispersionSpec(mode=mode, keyword=keyword)


def _to_pno_spec(raw: Mapping[str, Any], default: PNOSpec) -> PNOSpec:
    preset_raw = raw.get("preset", default.preset)
    preset = None if preset_raw is None else str(preset_raw).strip() or None
    return PNOSpec(preset=preset)


class MethodRegistry:
    _PROFILES: Dict[str, MethodProfile] = {
        "PWPB95": MethodProfile(
            family="double_hybrid_dft",
            base_method="PWPB95",
            scf_default=SCFSpec(accel="RIJCOSX", aux_j="def2/J"),
            correlation_default=CorrelationSpec(model="RI-MP2", aux_c=None),
            dispersion_default=DispersionSpec(mode="required", keyword="D4"),
            gradients_supported=True,
        ),
        "WB97X-2": MethodProfile(
            family="double_hybrid_dft",
            base_method="WB97X-2",
            scf_default=SCFSpec(accel="RIJCOSX", aux_j="def2/J"),
            correlation_default=CorrelationSpec(model="RI-MP2", aux_c=None),
            dispersion_default=DispersionSpec(mode="required", keyword="D4"),
            gradients_supported=True,
        ),
        "DSD-PBEP86": MethodProfile(
            family="double_hybrid_dft",
            base_method="DSD-PBEP86",
            scf_default=SCFSpec(accel="RIJCOSX", aux_j="def2/J"),
            correlation_default=CorrelationSpec(model="RI-MP2", aux_c=None),
            dispersion_default=DispersionSpec(mode="required", keyword="D4"),
            gradients_supported=True,
        ),
        "B2PLYP": MethodProfile(
            family="double_hybrid_dft",
            base_method="B2PLYP",
            scf_default=SCFSpec(accel="RIJCOSX", aux_j="def2/J"),
            correlation_default=CorrelationSpec(model="RI-MP2", aux_c=None),
            dispersion_default=DispersionSpec(mode="required", keyword="D4"),
            gradients_supported=True,
        ),
        "DSD-BLYP": MethodProfile(
            family="double_hybrid_dft",
            base_method="DSD-BLYP",
            scf_default=SCFSpec(accel="RIJCOSX", aux_j="def2/J"),
            correlation_default=CorrelationSpec(model="RI-MP2", aux_c=None),
            dispersion_default=DispersionSpec(mode="required", keyword="D4"),
            gradients_supported=True,
        ),
        "WB97M-V": MethodProfile(
            family="vv10_family",
            base_method="WB97M-V",
            scf_default=SCFSpec(accel="RIJCOSX", aux_j="def2/J"),
            correlation_default=CorrelationSpec(model="none", aux_c=None),
            dispersion_default=DispersionSpec(mode="forbidden", keyword=None),
            gradients_supported=True,
        ),
        "R2SCAN-3C": MethodProfile(
            family="composite_3c",
            base_method="r2SCAN-3c",
            scf_default=SCFSpec(accel="none", aux_j=None),
            correlation_default=CorrelationSpec(model="none", aux_c=None),
            dispersion_default=DispersionSpec(mode="forbidden", keyword=None),
            self_contained=True,
            gradients_supported=True,
        ),
        "B97-3C": MethodProfile(
            family="composite_3c",
            base_method="B97-3c",
            scf_default=SCFSpec(accel="none", aux_j=None),
            correlation_default=CorrelationSpec(model="none", aux_c=None),
            dispersion_default=DispersionSpec(mode="forbidden", keyword=None),
            self_contained=True,
            gradients_supported=True,
        ),
        "GFN2-XTB": MethodProfile(
            family="semiempirical_xtb",
            base_method="GFN2-xTB",
            scf_default=SCFSpec(accel="none", aux_j=None),
            correlation_default=CorrelationSpec(model="none", aux_c=None),
            dispersion_default=DispersionSpec(mode="forbidden", keyword=None),
            self_contained=True,
            gradients_supported=True,
        ),
        "DLPNO-CCSD(T)": MethodProfile(
            family="post_hf_dlpno",
            base_method="DLPNO-CCSD(T)",
            scf_default=SCFSpec(accel="RIJK", aux_j="def2/JK"),
            correlation_default=CorrelationSpec(model="DLPNO-CCSD(T)", aux_c=None),
            dispersion_default=DispersionSpec(mode="forbidden", keyword=None),
            gradients_supported=False,
            pno_default=PNOSpec(preset="NormalPNO"),
        ),
        "M062X": MethodProfile(
            family="hybrid_gga",
            base_method="M062X",
            scf_default=SCFSpec(accel="RIJCOSX", aux_j="def2/J"),
            correlation_default=CorrelationSpec(model="none", aux_c=None),
            dispersion_default=DispersionSpec(mode="none", keyword=None),
            gradients_supported=True,
        ),
        "WB97X-D4": MethodProfile(
            family="hybrid_gga",
            base_method="wB97X-D4",
            scf_default=SCFSpec(accel="RIJCOSX", aux_j="def2/J"),
            correlation_default=CorrelationSpec(model="none", aux_c=None),
            dispersion_default=DispersionSpec(mode="auto_builtin", keyword=None),
            gradients_supported=True,
        ),
        "WB97X-V": MethodProfile(
            family="hybrid_gga",
            base_method="wB97X-V",
            scf_default=SCFSpec(accel="RIJCOSX", aux_j="def2/J"),
            correlation_default=CorrelationSpec(model="none", aux_c=None),
            dispersion_default=DispersionSpec(mode="forbidden", keyword=None),
            gradients_supported=True,
        ),
    }

    @classmethod
    def get_profile(cls, method_alias: str, engine: str = "orca") -> MethodProfile:
        if str(engine).strip().lower() != "orca":
            raise ValueError(f"MethodRegistry currently supports ORCA only, got: {engine}")
        normalized, dispersion = KeywordTranslator.to_orca_method(method_alias)
        key = str(normalized or method_alias).upper()
        profile = cls._PROFILES.get(key)
        if profile is not None:
            if dispersion and profile.dispersion_default.mode in {"none", "external", "required"}:
                return MethodProfile(
                    family=profile.family,
                    base_method=profile.base_method,
                    scf_default=profile.scf_default,
                    correlation_default=profile.correlation_default,
                    dispersion_default=DispersionSpec(mode="external", keyword=dispersion),
                    self_contained=profile.self_contained,
                    gradients_supported=profile.gradients_supported,
                    forbidden_ri_for_gradients=profile.forbidden_ri_for_gradients,
                    pno_default=profile.pno_default,
                )
            return profile

        fallback_dispersion = DispersionSpec(mode="none", keyword=None)
        if dispersion:
            fallback_dispersion = DispersionSpec(mode="external", keyword=dispersion)
        return MethodProfile(
            family="hybrid_gga",
            base_method=str(normalized or method_alias),
            scf_default=SCFSpec(accel="RIJCOSX", aux_j="def2/J"),
            correlation_default=CorrelationSpec(model="none", aux_c=None),
            dispersion_default=fallback_dispersion,
            gradients_supported=True,
        )

    @classmethod
    def normalize_spec(
        cls,
        method_spec: Mapping[str, Any],
        *,
        default_engine: str = "orca",
        default_basis: str = "def2-TZVPP",
    ) -> NormalizedMethodSpec:
        raw = dict(method_spec)
        engine = str(raw.get("engine", default_engine) or default_engine).strip().lower()
        method_alias = str(raw.get("method", "") or "").strip()
        if not method_alias:
            raise ValueError("method is required in method spec")

        profile = cls.get_profile(method_alias, engine=engine)

        family = str(raw.get("family", profile.family) or profile.family)
        basis = str(raw.get("basis", default_basis) or default_basis).strip()
        aux_basis_raw = raw.get("aux_basis")
        aux_basis = None if aux_basis_raw is None else str(aux_basis_raw).strip() or None

        scf_any = raw.get("scf")
        corr_any = raw.get("correlation")
        disp_any = raw.get("dispersion")
        pno_any = raw.get("pno")

        scf_raw: Mapping[str, Any] = scf_any if isinstance(scf_any, Mapping) else {}
        corr_raw: Mapping[str, Any] = corr_any if isinstance(corr_any, Mapping) else {}
        disp_raw: Mapping[str, Any] = disp_any if isinstance(disp_any, Mapping) else {}
        pno_raw: Mapping[str, Any] = pno_any if isinstance(pno_any, Mapping) else {}

        scf = _to_scf_spec(scf_raw, profile.scf_default)
        correlation = _to_correlation_spec(corr_raw, profile.correlation_default, basis)
        dispersion = _to_dispersion_spec(disp_raw, profile.dispersion_default)
        pno = _to_pno_spec(pno_raw, profile.pno_default)

        # Self-contained methods (e.g., composite_3c like R2SCAN-3C) must not inherit
        # aux_basis from upstream merged config - force to None for proper ORCA method spec.
        if profile.self_contained:
            aux_basis = None

        route_extras = str(raw.get("route_extras", "") or "").strip()
        solvent_raw = raw.get("solvent")
        solvent = None if solvent_raw is None else str(solvent_raw).strip() or None

        maxcore = raw.get("maxcore")
        maxcore_int = int(maxcore) if isinstance(maxcore, int) or (isinstance(maxcore, str) and maxcore.isdigit()) else None

        return NormalizedMethodSpec(
            engine=engine,
            family=family,
            method=profile.base_method,
            basis=basis,
            scf=scf,
            correlation=correlation,
            dispersion=dispersion,
            pno=pno,
            route_extras=route_extras,
            solvent=solvent,
            maxcore=maxcore_int,
            aux_basis=aux_basis,
            is_reference=bool(raw.get("is_reference", False)),
            reuse_baseline_gcorr=bool(raw.get("reuse_baseline_gcorr", True)),
            raw=raw,
        )
