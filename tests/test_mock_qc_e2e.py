import pytest
import importlib

from rph_core.steps.step4_features.mech_packager import pack_mechanism_assets


def test_pack_mechanism_assets_importable() -> None:
    assert callable(pack_mechanism_assets)


def test_fake_backend_removed_from_public_api() -> None:
    module = importlib.import_module("rph_core.utils.qc_interface")
    assert not hasattr(module, "FakeBackend")
