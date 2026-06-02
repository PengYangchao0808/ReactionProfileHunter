import logging
from typing import Tuple, Optional

logger = logging.getLogger(__name__)

class KeywordTranslator:
    """
    Translates generic computational chemistry terms into software-specific keywords.
    Current support: Gaussian 16, ORCA.
    """

    # ORCA built-in dispersion functionals (already include dispersion)
    ORCA_BUILTIN_DISPERSION = frozenset({
        "wB97X-D3BJ", "wB97X-D4",
        "wB97M-D3BJ", "wB97M-D4",
        "B97M-D3BJ", "B97M-D4",
    })

    # ORCA method name mapping: generic -> ORCA
    ORCA_METHOD_MAP = {
        "M06-2X": "M062X",
        "M06-L": "M06L",
        "wB97X-2": "WB97X-2",
        "WB97X-2": "WB97X-2",
        "wB97M-V": "WB97M-V",
        "WB97M-V": "WB97M-V",
    }

    # Dispersion suffixes to strip from method names (ordered: longest first)
    DISPERSION_SUFFIXES = ["-D3BJ", "-D3ZERO", "-D4", "-D3"]

    @staticmethod
    def to_gaussian_basis(basis_raw: str) -> str:
        """
        Gaussian requires 'Def2SVP' (no hyphen) for def2 family.
        Generic 'def2-SVP' -> 'Def2SVP'.
        """
        if not basis_raw: return "Def2SVP" # Safety default
        # Remove hyphens/underscores ONLY for def2 series
        if "def2" in basis_raw.lower():
            return basis_raw.replace("-", "").replace("_", "")
        return basis_raw

    @staticmethod
    def to_gaussian_dispersion(dispersion_raw: str) -> str:
        """
        Maps generic 'GD3BJ' -> Gaussian 'em=GD3BJ' (Concise format).
        """
        if not dispersion_raw: return ""
        d = dispersion_raw.upper()
        if d in ["GD3BJ", "D3BJ"]: return "em=GD3BJ"
        if d in ["GD3", "D3"]: return "em=GD3"
        return ""

    @staticmethod
    def to_gaussian_solvent(solvent_raw: str) -> str:
        from rph_core.utils.solvent_map import gaussian_pcm_keyword

        return gaussian_pcm_keyword(solvent_raw)

    @classmethod
    def to_orca_method(cls, method_name: str) -> Tuple[str, Optional[str]]:
        """
        Convert a generic method name to ORCA-compatible format.

        Handles:
        - Method name normalization (e.g., M06-2X -> M062X)
        - Dispersion suffix extraction (e.g., PWPB95-D4 -> PWPB95 + D4)
        - Built-in dispersion functionals (e.g., wB97X-D3BJ -> wB97X-D3BJ, no extra dispersion)

        Returns:
            Tuple of (orca_method_name, dispersion_keyword_or_None)
        """
        if not method_name:
            return method_name, None

        # 1. Check if it's a built-in dispersion functional
        if method_name in cls.ORCA_BUILTIN_DISPERSION:
            return method_name, None

        # 2. Check direct mapping first
        if method_name in cls.ORCA_METHOD_MAP:
            return cls.ORCA_METHOD_MAP[method_name], None

        # 3. Try to strip dispersion suffixes
        base_method = method_name
        dispersion = None

        for suffix in cls.DISPERSION_SUFFIXES:
            if base_method.endswith(suffix):
                base_method = base_method[:-len(suffix)]
                dispersion = suffix[1:]  # Remove leading '-'
                break

        # 4. Map base method if needed
        if base_method in cls.ORCA_METHOD_MAP:
            base_method = cls.ORCA_METHOD_MAP[base_method]

        # 5. Validate result
        if not base_method:
            logger.warning(f"Empty method name after normalization: {method_name}")
            return method_name, None

        return base_method, dispersion
