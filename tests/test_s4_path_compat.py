from pathlib import Path

from rph_core.steps.step4_features.mech_packager import MechanismContext, S1_DIR_ALIASES


def test_s4_uses_canonical_s1_alias_only(tmp_path: Path) -> None:
    s1_dir = tmp_path / "S1_ConfGeneration"
    s1_dir.mkdir(parents=True)

    context = MechanismContext(s1_dir=s1_dir)

    assert S1_DIR_ALIASES == ["S1_ConfGeneration"]
    assert context.s1_dir == s1_dir


def test_s4_rejects_legacy_s1_alias(tmp_path: Path) -> None:
    legacy_dir = tmp_path / "S1_Product"
    legacy_dir.mkdir(parents=True)

    assert "S1_Product" not in S1_DIR_ALIASES
