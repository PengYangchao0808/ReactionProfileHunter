import json
import logging
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from rph_core.steps.runners import run_step2, run_step3, run_step4


def _write_xyz(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("3\nxyz\nC 0 0 0\nH 0 0 1\nH 0 1 0\n", encoding="utf-8")


def _make_hunter(tmp_path: Path) -> SimpleNamespace:
    hunter = SimpleNamespace()
    hunter.logger = logging.getLogger("test_phase3c")
    hunter.config = {
        "step2": {
            "pes_adapter": {"enabled": True, "mode": "fallback"},
            "path_search": {"enabled": False},
        }
    }
    hunter._resolve_product_xyz_for_s2 = MagicMock()
    hunter._resolve_profile_key = MagicMock(return_value="[4+3]_default")
    hunter._resolve_forming_bonds_for_s2 = MagicMock(return_value=((0, 1), (2, 3)))
    hunter._resolve_forward_scan_config = MagicMock(return_value={"scan_start_distance": 3.5, "scan_end_distance": 1.8})
    hunter._build_step2_signature = MagicMock(return_value={"sig": "ok"})
    hunter._resolve_s1_artifacts = MagicMock(
        return_value={
            "s1_dir": tmp_path / "S1_ConfGeneration",
            "s1_shermo_summary_file": None,
            "s1_hoac_thermo_file": None,
            "s1_conformer_energies_file": None,
            "s1_precursor_xyz": None,
        }
    )
    hunter.s2_engine = SimpleNamespace()
    hunter.s3_engine = SimpleNamespace()
    hunter.s4_engine = SimpleNamespace()
    return hunter


@pytest.mark.parametrize("protocol", ["lite", "zero"])
def test_phase3c_full_chain_no_fallback_when_dft_geometry(protocol: str, tmp_path: Path) -> None:
    s1_dir = tmp_path / "S1_ConfGeneration"
    s1_dir.mkdir(parents=True, exist_ok=True)
    (s1_dir / "provenance.json").write_text(
        json.dumps({"protocol": protocol, "has_geometry_optimization": True}, indent=2),
        encoding="utf-8",
    )

    product = s1_dir / "product_min.xyz"
    ts_guess = tmp_path / "S2_Retro" / "ts_guess.xyz"
    substrate = tmp_path / "S2_Retro" / "reactant_complex.xyz"
    intermediate = tmp_path / "S2_Retro" / "intermediate.xyz"
    ts_final = tmp_path / "S3_TS" / "ts_final.xyz"
    features = tmp_path / "S4_Data" / "features_raw.csv"

    for f in [product, ts_guess, substrate, intermediate, ts_final]:
        _write_xyz(f)
    features.parent.mkdir(parents=True, exist_ok=True)
    features.write_text("id,value\nrx,1\n", encoding="utf-8")

    hunter = _make_hunter(tmp_path)
    hunter._resolve_product_xyz_for_s2.return_value = product
    hunter.s2_engine.run_retro_scan = MagicMock(
        return_value=(
            ts_guess,
            substrate,
            intermediate,
            ((0, 1), (2, 3)),
            tmp_path / "S2_Retro" / "scan_profile.json",
            "COMPLETE",
            "high",
            tuple(),
            None,
        )
    )
    hunter.s3_engine.run = MagicMock(
        return_value=SimpleNamespace(
            ts_final_xyz=ts_final,
            sp_report=SimpleNamespace(method="mock"),
            ts_fchk=None,
            ts_log=None,
            ts_qm_output=None,
            intermediate_fchk=None,
            intermediate_log=None,
            intermediate_qm_output=None,
        )
    )
    hunter.s4_engine.run = MagicMock(return_value=features)

    step2 = run_step2(
        hunter=hunter,
        product_xyz=product,
        work_dir=tmp_path,
        reaction_profile="[4+3]_default",
        cleaner_data=None,
    )
    assert "pes_adapter_fallback_applied" not in step2.degraded_reasons
    assert not (tmp_path / "S2_Retro" / "pes_adapter_fallback" / "adapter_meta.json").exists()

    called_product = hunter.s2_engine.run_retro_scan.call_args.kwargs["product_xyz"]
    assert Path(called_product).name == "product_min.xyz"

    step3 = run_step3(
        hunter=hunter,
        ts_guess_xyz=step2.ts_guess_xyz,
        intermediate_xyz=step2.intermediate_xyz,
        product_xyz=product,
        work_dir=tmp_path,
        e_product_l2=-1.0,
        product_thermo=None,
        forming_bonds=step2.forming_bonds,
        old_checkpoint=None,
    )
    assert step3.ts_final_xyz == ts_final

    step4 = run_step4(
        hunter=hunter,
        ts_final_xyz=step3.ts_final_xyz,
        substrate_xyz=step2.substrate_xyz,
        product_xyz=product,
        work_dir=tmp_path,
        forming_bonds=step2.forming_bonds,
        sp_matrix_report=step3.sp_report,
        ts_fchk=step3.ts_fchk,
        intermediate_fchk=step3.intermediate_fchk,
        product_fchk=None,
        ts_log=step3.ts_log,
        intermediate_log=step3.intermediate_log,
        product_log=None,
        ts_qm_output=step3.ts_qm_output,
        intermediate_qm_output=step3.intermediate_qm_output,
        product_qm_output=None,
    )
    assert step4.features_csv == features


def test_phase3c_chain_survives_fallback_path(tmp_path: Path) -> None:
    s1_dir = tmp_path / "S1_ConfGeneration"
    s1_dir.mkdir(parents=True, exist_ok=True)
    (s1_dir / "provenance.json").write_text(
        json.dumps(
            {
                "schema_version": "s1_provenance_v1",
                "protocol": "zero",
                "has_geometry_optimization": False,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    product = s1_dir / "product_min.xyz"
    ts_guess = tmp_path / "S2_Retro" / "ts_guess.xyz"
    substrate = tmp_path / "S2_Retro" / "reactant_complex.xyz"
    intermediate = tmp_path / "S2_Retro" / "intermediate.xyz"
    ts_final = tmp_path / "S3_TS" / "ts_final.xyz"
    features = tmp_path / "S4_Data" / "features_raw.csv"

    for f in [product, ts_guess, substrate, intermediate, ts_final]:
        _write_xyz(f)
    features.parent.mkdir(parents=True, exist_ok=True)
    features.write_text("id,value\nrx,1\n", encoding="utf-8")

    hunter = _make_hunter(tmp_path)
    hunter._resolve_product_xyz_for_s2.return_value = product
    hunter.s2_engine.run_retro_scan = MagicMock(
        return_value=(
            ts_guess,
            substrate,
            intermediate,
            ((0, 1), (2, 3)),
            tmp_path / "S2_Retro" / "scan_profile.json",
            "COMPLETE",
            "high",
            tuple(),
            None,
        )
    )
    hunter.s3_engine.run = MagicMock(
        return_value=SimpleNamespace(
            ts_final_xyz=ts_final,
            sp_report=SimpleNamespace(method="mock"),
            ts_fchk=None,
            ts_log=None,
            ts_qm_output=None,
            intermediate_fchk=None,
            intermediate_log=None,
            intermediate_qm_output=None,
        )
    )
    hunter.s4_engine.run = MagicMock(return_value=features)

    step2 = run_step2(
        hunter=hunter,
        product_xyz=product,
        work_dir=tmp_path,
        reaction_profile="[4+3]_default",
        cleaner_data=None,
    )
    assert "pes_adapter_fallback_applied" in step2.degraded_reasons
    assert "pes_adapter_fallback_exception_path" in step2.degraded_reasons

    adapter_meta_file = tmp_path / "S2_Retro" / "pes_adapter_fallback" / "adapter_meta.json"
    assert adapter_meta_file.exists()
    adapter_meta = json.loads(adapter_meta_file.read_text(encoding="utf-8"))
    assert adapter_meta["protocol"] == "zero"
    assert adapter_meta["trigger_source"] == "s1_provenance.has_geometry_optimization=false"
    assert adapter_meta["provenance_schema_version"] == "s1_provenance_v1"

    called_product = hunter.s2_engine.run_retro_scan.call_args.kwargs["product_xyz"]
    assert Path(called_product).name == "product_relaxed.xyz"

    step3 = run_step3(
        hunter=hunter,
        ts_guess_xyz=step2.ts_guess_xyz,
        intermediate_xyz=step2.intermediate_xyz,
        product_xyz=product,
        work_dir=tmp_path,
        e_product_l2=-1.0,
        product_thermo=None,
        forming_bonds=step2.forming_bonds,
        old_checkpoint=None,
    )
    step4 = run_step4(
        hunter=hunter,
        ts_final_xyz=step3.ts_final_xyz,
        substrate_xyz=step2.substrate_xyz,
        product_xyz=product,
        work_dir=tmp_path,
        forming_bonds=step2.forming_bonds,
        sp_matrix_report=step3.sp_report,
        ts_fchk=step3.ts_fchk,
        intermediate_fchk=step3.intermediate_fchk,
        product_fchk=None,
        ts_log=step3.ts_log,
        intermediate_log=step3.intermediate_log,
        product_log=None,
        ts_qm_output=step3.ts_qm_output,
        intermediate_qm_output=step3.intermediate_qm_output,
        product_qm_output=None,
    )
    assert step4.features_csv == features
