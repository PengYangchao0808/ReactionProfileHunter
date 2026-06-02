from dataclasses import dataclass
from typing import Optional, Any, Dict


@dataclass(frozen=True)
class ShermoOptions:
    temperature_k: float = 298.15
    pressure_atm: Optional[float] = None
    scl_zpe: Optional[float] = None
    ilowfreq: Optional[int] = None
    imagreal: Optional[float] = None
    conc: Optional[str] = None

    @classmethod
    def from_config(cls, config: Dict[str, Any], temperature_k: Optional[float] = None) -> "ShermoOptions":
        thermo_cfg = config.get("thermo", {}) or {}

        def _to_float(value) -> Optional[float]:
            if value is None:
                return None
            try:
                return float(str(value).strip())
            except (TypeError, ValueError):
                return None

        def _to_int(value) -> Optional[int]:
            if value is None:
                return None
            try:
                return int(str(value).strip())
            except (TypeError, ValueError):
                return None

        return cls(
            temperature_k=temperature_k if temperature_k is not None else float(thermo_cfg.get("temperature_k", 298.15) or 298.15),
            pressure_atm=_to_float(thermo_cfg.get("pressure_atm")),
            scl_zpe=_to_float(thermo_cfg.get("scl_zpe")),
            ilowfreq=_to_int(thermo_cfg.get("ilowfreq")),
            imagreal=_to_float(thermo_cfg.get("imagreal")),
            conc=thermo_cfg.get("conc"),
        )
