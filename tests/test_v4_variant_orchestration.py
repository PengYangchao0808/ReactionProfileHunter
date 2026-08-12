from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from pathlib import Path

import pytest
from _pytest.monkeypatch import MonkeyPatch

import rph_core.v4_orchestrator as v4_orchestrator
from rph_core.steps.mechanism_classifier.s0_record import S0ReactionRecord
from rph_core.steps.refinement import RefinementEngine
from rph_core.steps.refinement.manifest_io import write_refinement_manifest
from rph_core.utils.variant_checkpoint import VariantCheckpoint


class StubRecord(S0ReactionRecord):
    def __init__(self, rx_id: str = "RXN_P2", product_smiles: str = "C1CCC1"):
        super().__init__(
            rx_id=rx_id,
            product_smiles=product_smiles,
            reaction_type="4+3",
            mapped_product_smiles="[CH2:1]1[CH2:2][CH2:3][CH2:4]1",
            mapped_forming_bonds=((1, 2), (3, 4)),
            forming_bonds=((0, 1), (2, 3)),
            mapping_confidence=None,
            mapping_trusted=True,
            source_csv=Path("dataset.csv"),
            canonical_precursor_smiles="",
            mapped_precursor_smiles="",
        )


def _orchestrator() -> v4_orchestrator.V4Orchestrator:
    orchestrator = object.__new__(v4_orchestrator.V4Orchestrator)
    orchestrator.config = {
        "step1": {"protocol": "censo_lite"},
        "step2": {},
        "theory": {
            "s3_low_level": {},
            "s4_high_precision": {},
        },
    }
    orchestrator.color = True
    orchestrator.ui_mode = None
    orchestrator.resume_policy = "strict"
    orchestrator.start_from = None
    orchestrator.recompute_from = None
    orchestrator.refresh_stages = set()
    orchestrator.rerun_failed_s3 = False
    orchestrator.s3_structure_ids = ()
    orchestrator.run_id = None
    return orchestrator


def _fake_s0_writer(variants: Sequence[Mapping[str, object]]):
    def _write_s0(work_dir: Path, record: StubRecord) -> Path:
        stage_dir = Path(work_dir) / "S0_Mechanism"
        _ = stage_dir.mkdir(parents=True, exist_ok=True)
        mechanism_path = stage_dir / "mechanism.json"
        mapping_dir = stage_dir / "atom_mappings"
        mapping_dir.mkdir(parents=True, exist_ok=True)
        registry_variants = []
        for variant in variants:
            variant_id = str(variant["variant_id"])
            mapping_ref = f"atom_mappings/{variant_id}.json"
            (stage_dir / mapping_ref).write_text(
                json.dumps({
                    "schema_version": "s0_atom_map_smiles_v2",
                    "product_smiles_idx_space": "geometry_product_smiles_idx",
                    "forming_bonds_map_space": [[1, 2], [3, 4]],
                    "map_to_product_smiles": {"1": 0, "2": 1, "3": 2, "4": 3},
                }),
                encoding="utf-8",
            )
            registry_variants.append({**variant, "atom_mapping_ref": mapping_ref})
        _ = mechanism_path.write_text(
            json.dumps(
                {
                    "schema_version": "s0_mechanism_v2",
                    "stage": "S0",
                    "source": "trusted_reaction_record",
                    "rx_id": record.rx_id,
                    "forming_bonds": [list(pair) for pair in record.forming_bonds],
                    "variants": [variant["variant_id"] for variant in variants],
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        _ = (stage_dir / "variant_registry.json").write_text(
            json.dumps(
                {
                    "schema_version": "s0_variant_registry_v1",
                    "reaction_id": record.rx_id,
                    "variants": registry_variants,
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        return mechanism_path

    return _write_s0


def _install_variant_fakes(
    monkeypatch: MonkeyPatch,
    *,
    s2_fail_variant: str | None = None,
    shared_ts_seed_variant: str | None = None,
) -> tuple[list[str], list[str], list[dict[str, object]]]:
    s1_calls: list[str] = []
    s2_calls: list[str] = []
    s3_structures: list[dict[str, object]] = []

    class FakeCensoLiteEngine:
        def __init__(self, _config: object, stage_dir: Path, molecule_name: str):
            self.stage_dir: Path = Path(stage_dir)
            self.molecule_name: str = molecule_name

        def run(self, smiles: str) -> dict[str, object]:
            s1_calls.append(self.molecule_name)
            molecule_dir = self.stage_dir / self.molecule_name
            candidate_dir = molecule_dir / "candidates"
            _ = candidate_dir.mkdir(parents=True, exist_ok=True)
            xyz = candidate_dir / "conf_0001.xyz"
            _ = xyz.write_text(
                f"4\n{smiles}\nC 0 0 0\nC 1 0 0\nC 1 1 0\nC 0 1 0\n",
                encoding="utf-8",
            )
            mapping = molecule_dir / "atom_mapping.json"
            mapping.write_text(json.dumps({
                "schema_version": "s1_atom_mapping_v1",
                "mapping_source": "rdkit_generation_sidecar",
                "confidence": "high",
                "product_smiles_idx_space": "geometry_product_smiles_idx",
                "atoms": [
                    {"smiles_idx": i, "xyz_idx": i, "element": "C", "type": "organic_heavy"}
                    for i in range(4)
                ],
            }), encoding="utf-8")
            manifest = molecule_dir / "manifest.json"
            payload = {
                "schema_version": "s1_censo_light_ranking_v4",
                "selected": f"{self.molecule_name}_conf_0001",
                "atom_mapping_ref": "atom_mapping.json",
                "atom_mapping_ref": "atom_mapping.json",
                "candidates": [
                    {
                        "id": f"{self.molecule_name}_conf_0001",
                        "xyz": f"{self.molecule_name}/candidates/conf_0001.xyz",
                    }
                ],
            }
            _ = manifest.write_text(json.dumps(payload, indent=2), encoding="utf-8")
            return {"manifest": manifest, "selected_xyz": xyz, "data": payload}

    class FakePEBScanner:
        def __init__(self, _config: object, molecule_name: str | None = None):
            self.molecule_name: str = molecule_name or "unknown"

        def run(
            self,
            product_xyz: Path,
            output_dir: Path,
            forming_bonds: Sequence[tuple[int, int]],
            scan_config: dict[str, object] | None = None,
        ) -> tuple[Path, Path, Path, tuple[tuple[int, int], ...], Path, str, str, tuple[()]]:
            del scan_config
            s2_calls.append(self.molecule_name)
            if self.molecule_name == s2_fail_variant:
                raise RuntimeError(f"forced S2 failure for {self.molecule_name}")
            output_dir = Path(output_dir)
            _ = output_dir.mkdir(parents=True, exist_ok=True)
            ts_seed = output_dir / "ts_guess.xyz"
            intermediate_seed = output_dir / "intermediate_seed.xyz"
            profile = output_dir / "scan_profile.json"
            contents = Path(product_xyz).read_text(encoding="utf-8")
            _ = ts_seed.write_text(contents, encoding="utf-8")
            _ = intermediate_seed.write_text(contents, encoding="utf-8")
            shared_ts_seed = self.molecule_name == shared_ts_seed_variant
            _ = profile.write_text(
                json.dumps(
                    {
                        "selections": {
                            "intermediate": {
                                "selection_mode": (
                                    "shared_ts_fallback" if shared_ts_seed else "stable_basin"
                                ),
                                "selection_status": (
                                    "shared_with_ts" if shared_ts_seed else "selected"
                                ),
                                "rule": "no_resolved_pre_ts_basin_shared_ts_seed",
                            }
                        }
                    }
                ),
                encoding="utf-8",
            )
            return (
                ts_seed,
                intermediate_seed,
                None if shared_ts_seed else intermediate_seed,
                tuple((int(pair[0]), int(pair[1])) for pair in forming_bonds),
                profile,
                "COMPLETE",
                "high",
                (),
            )

    class FakeRefinementEngine(RefinementEngine):
        def run(
            self,
            structures: Sequence[dict[str, object]],
            output_dir: Path,
            **_kwargs: object,
        ) -> Path:
            rows = list(structures)
            s3_structures.extend(rows)
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            return write_refinement_manifest(
                output_dir / "manifest.json",
                stage=self.profile.stage,
                fidelity=self.profile.fidelity,
                profile_id=self.profile.profile_id,
                structures=[
                    {
                        "id": row["id"],
                        "kind": row["kind"],
                        "role": row.get("role"),
                        "variant_id": row.get("variant_id"),
                        "branch_id": row.get("branch_id"),
                        "pathway_id": row.get("pathway_id"),
                        "parent_structure_id": row.get("parent_structure_id"),
                        "source_stage": row.get("source_stage"),
                        "forming_bonds": row.get("forming_bonds"),
                        "input_xyz": row["input_xyz"],
                        "opt_xyz": row["input_xyz"],
                        "opt_status": "complete",
                        "frequency_status": "complete" if row.get("kind") == "ts" else "complete",
                        "sp_status": "complete",
                        "status": "complete",
                        "usable_for_ml": True,
                    }
                    for row in rows
                ],
                run_id=self.run_id,
                extra={"summary": {"complete": len(rows)}},
            )

    monkeypatch.setattr(v4_orchestrator, "CensoLiteEngine", FakeCensoLiteEngine)
    monkeypatch.setattr(v4_orchestrator, "PEBScanner", FakePEBScanner)
    monkeypatch.setattr(v4_orchestrator, "RefinementEngine", FakeRefinementEngine)
    return s1_calls, s2_calls, s3_structures


def test_variant_checkpoint_tracks_per_variant_state(tmp_path: Path):
    checkpoint = VariantCheckpoint(tmp_path)
    checkpoint.set_reaction_id("RXN_P2")

    manifest = tmp_path / "S1_ConfSearch" / "product_major" / "manifest.json"
    _ = manifest.parent.mkdir(parents=True)
    _ = manifest.write_text("{}", encoding="utf-8")
    signature = checkpoint.signature({"stage": "s1", "variant": "product_major"})
    checkpoint.mark("product_major", "s1", signature, manifest)
    checkpoint.mark_failed("product_major", "s2", "sig_s2", "boom")

    saved = json.loads((tmp_path / "pipeline.state").read_text(encoding="utf-8"))
    assert checkpoint.reusable("product_major", "s1", signature, manifest)
    assert saved["reaction_id"] == "RXN_P2"
    assert saved["stages"]["s1"]["scopes"]["product_major"]["status"] == "complete"
    assert saved["stages"]["s2"]["scopes"]["product_major"]["status"] == "failed"
    assert not (tmp_path / ".rph" / "checkpoint.json").exists()


def test_single_variant_run_writes_run_manifest(monkeypatch: MonkeyPatch, tmp_path: Path):
    variants = [
        {
            "variant_id": "product_major",
            "directory_name": "product_major",
            "branch_id": "BR_MAJOR",
            "branch_product_smiles": "C1CCC1",
        }
    ]
    monkeypatch.setattr(
        v4_orchestrator.V4Orchestrator,
        "_write_s0_from_record",
        staticmethod(_fake_s0_writer(variants)),
    )
    _install_variant_fakes(monkeypatch)

    result = _orchestrator().run(StubRecord(), tmp_path / "run", stop_after="s2")

    run_manifest = json.loads((tmp_path / "run" / "run.manifest.json").read_text(encoding="utf-8"))
    assert result["completed_through"] == "s2"
    assert "s1_manifest" in result
    assert "s2_manifest" in result
    assert run_manifest["schema_version"] == "rph_run_manifest_v1"
    assert run_manifest["variants"] == ["product_major"]
    assert run_manifest["variant_status"]["product_major"] == "s2_complete"
    assert run_manifest["stages"]["s1"]["product_major"]["status"] == "complete"
    assert run_manifest["stages"]["s2"]["product_major"]["status"] == "complete"


def test_variant_registry_drives_all_variants_into_s3_with_prefixed_ids(monkeypatch: MonkeyPatch, tmp_path: Path):
    variants = [
        {
            "variant_id": "product_major",
            "directory_name": "product_major",
            "branch_id": "BR_MAJOR",
            "branch_product_smiles": "C1CCC1",
        },
        {
            "variant_id": "product_minor_001",
            "directory_name": "product_minor_001",
            "branch_id": "BR_DR_001",
            "branch_product_smiles": "C/C=C\\C",
        },
    ]
    monkeypatch.setattr(
        v4_orchestrator.V4Orchestrator,
        "_write_s0_from_record",
        staticmethod(_fake_s0_writer(variants)),
    )
    s1_calls, s2_calls, s3_structures = _install_variant_fakes(monkeypatch)

    result = _orchestrator().run(StubRecord(), tmp_path / "run", stop_after="s3")

    assert result["completed_through"] == "s3"
    assert s1_calls == ["product_major", "product_minor_001"]
    assert s2_calls == ["product_major", "product_minor_001"]
    assert [row["id"] for row in s3_structures] == [
        "product_major",
        "product_major_int",
        "product_major_ts",
        "product_minor_001",
        "product_minor_001_int",
        "product_minor_001_ts",
    ]
    assert all(
        "forming_bonds" not in row
        for row in s3_structures
        if row.get("role") == "product"
    )
    assert all(
        row.get("forming_bonds") == [[0, 1], [2, 3]]
        for row in s3_structures
        if row.get("role") == "intermediate"
    )
    assert [row["forming_bonds"] for row in s3_structures if row["kind"] == "ts"] == [
        [[0, 1], [2, 3]],
        [[0, 1], [2, 3]],
    ]

    run_manifest = json.loads((tmp_path / "run" / "run.manifest.json").read_text(encoding="utf-8"))
    assert run_manifest["variants"] == ["product_major", "product_minor_001"]
    assert run_manifest["variant_status"]["product_major"] == "s2_complete"
    assert run_manifest["variant_status"]["product_minor_001"] == "s2_complete"


def test_shared_ts_seed_is_promoted_to_s3_intermediate(monkeypatch: MonkeyPatch, tmp_path: Path):
    variants = [
        {
            "variant_id": "product_major",
            "directory_name": "product_major",
            "branch_id": "BR_MAJOR",
            "branch_product_smiles": "C1CCC1",
        }
    ]
    monkeypatch.setattr(
        v4_orchestrator.V4Orchestrator,
        "_write_s0_from_record",
        staticmethod(_fake_s0_writer(variants)),
    )
    _s1_calls, _s2_calls, s3_structures = _install_variant_fakes(
        monkeypatch,
        shared_ts_seed_variant="product_major",
    )

    result = _orchestrator().run(StubRecord(), tmp_path / "run", stop_after="s3")

    assert result["completed_through"] == "s3"
    structures = {str(row["id"]): row for row in s3_structures}
    shared_int = structures["product_major_int"]
    ts = structures["product_major_ts"]
    assert shared_int["role"] == "intermediate"
    assert shared_int["kind"] == "minimum"
    assert shared_int["input_xyz"] == ts["input_xyz"]
    assert shared_int["seed_state"] == "shared_ts_seed_for_intermediate"


def test_variant_failure_does_not_abort_other_variants(monkeypatch: MonkeyPatch, tmp_path: Path):
    variants = [
        {
            "variant_id": "product_major",
            "directory_name": "product_major",
            "branch_id": "BR_MAJOR",
            "branch_product_smiles": "C1CCC1",
        },
        {
            "variant_id": "product_minor_001",
            "directory_name": "product_minor_001",
            "branch_id": "BR_DR_001",
            "branch_product_smiles": "C/C=C\\C",
        },
    ]
    monkeypatch.setattr(
        v4_orchestrator.V4Orchestrator,
        "_write_s0_from_record",
        staticmethod(_fake_s0_writer(variants)),
    )
    _s1_calls, s2_calls, s3_structures = _install_variant_fakes(
        monkeypatch,
        s2_fail_variant="product_minor_001",
    )

    result = _orchestrator().run(StubRecord(), tmp_path / "run", stop_after="s3")

    assert result["completed_through"] == "s3"
    assert s2_calls == ["product_major", "product_minor_001"]
    assert [row["id"] for row in s3_structures] == [
        "product_major",
        "product_major_int",
        "product_major_ts",
    ]

    run_manifest = json.loads((tmp_path / "run" / "run.manifest.json").read_text(encoding="utf-8"))
    assert run_manifest["overall_status"] == "partial_failure"
    assert run_manifest["variant_status"]["product_major"] == "s2_complete"
    assert run_manifest["variant_status"]["product_minor_001"] == "failed"
    assert run_manifest["stages"]["s2"]["product_minor_001"]["status"] == "failed"


def test_unchanged_s1_is_reused_when_continuing_to_s2(
    monkeypatch: MonkeyPatch, tmp_path: Path
):
    variants = [{
        "variant_id": "product_major",
        "directory_name": "product_major",
        "branch_product_smiles": "C1CCC1",
    }]
    monkeypatch.setattr(
        v4_orchestrator.V4Orchestrator,
        "_write_s0_from_record",
        staticmethod(_fake_s0_writer(variants)),
    )
    s1_calls, s2_calls, _structures = _install_variant_fakes(monkeypatch)
    output = tmp_path / "run"

    _orchestrator().run(StubRecord(), output, stop_after="s1")
    _orchestrator().run(StubRecord(), output, stop_after="s2")
    _orchestrator().run(StubRecord(), output, stop_after="s2")

    assert s1_calls == ["product_major"]
    assert s2_calls == ["product_major"]
    assert (output / "run.config.json").exists()


def test_scientific_s1_change_fails_before_qc_under_strict_policy(
    monkeypatch: MonkeyPatch, tmp_path: Path
):
    variants = [{
        "variant_id": "product_major",
        "directory_name": "product_major",
        "branch_product_smiles": "C1CCC1",
    }]
    monkeypatch.setattr(
        v4_orchestrator.V4Orchestrator,
        "_write_s0_from_record",
        staticmethod(_fake_s0_writer(variants)),
    )
    s1_calls, _s2_calls, _structures = _install_variant_fakes(monkeypatch)
    output = tmp_path / "run"
    _orchestrator().run(StubRecord(), output, stop_after="s1")

    changed = _orchestrator()
    changed.config["step1"]["censo_lite"] = {
        "xtb_thermo": {"gfn_level": 2}
    }
    with pytest.raises(v4_orchestrator.ResumeCheckpointError, match="signature mismatch"):
        changed.run(StubRecord(), output, stop_after="s2")

    assert s1_calls == ["product_major"]


def test_explicit_start_from_s2_accepts_existing_s1_signature(
    monkeypatch: MonkeyPatch, tmp_path: Path
):
    variants = [{
        "variant_id": "product_major",
        "directory_name": "product_major",
        "branch_product_smiles": "C1CCC1",
    }]
    monkeypatch.setattr(
        v4_orchestrator.V4Orchestrator,
        "_write_s0_from_record",
        staticmethod(_fake_s0_writer(variants)),
    )
    s1_calls, s2_calls, _structures = _install_variant_fakes(monkeypatch)
    output = tmp_path / "run"
    _orchestrator().run(StubRecord(), output, stop_after="s1")
    before = json.loads((output / "pipeline.state").read_text(encoding="utf-8"))
    recorded_s1_signature = before["stages"]["s1"]["scopes"]["product_major"]["signature"]

    resumed = _orchestrator()
    resumed.config["step1"]["censo_lite"] = {
        "xtb_thermo": {"gfn_level": 2}
    }
    resumed.resume_policy = "use-existing-upstream"
    resumed.start_from = "s2"
    resumed.recompute_from = None
    resumed.run(StubRecord(), output, stop_after="s2")

    after = json.loads((output / "pipeline.state").read_text(encoding="utf-8"))
    assert s1_calls == ["product_major"]
    assert s2_calls == ["product_major"]
    assert after["stages"]["s1"]["scopes"]["product_major"]["signature"] == recorded_s1_signature


def test_explicit_recompute_from_s1_reruns_s1(
    monkeypatch: MonkeyPatch, tmp_path: Path
):
    variants = [{
        "variant_id": "product_major",
        "directory_name": "product_major",
        "branch_product_smiles": "C1CCC1",
    }]
    monkeypatch.setattr(
        v4_orchestrator.V4Orchestrator,
        "_write_s0_from_record",
        staticmethod(_fake_s0_writer(variants)),
    )
    s1_calls, _s2_calls, _structures = _install_variant_fakes(monkeypatch)
    output = tmp_path / "run"
    _orchestrator().run(StubRecord(), output, stop_after="s1")

    recompute = _orchestrator()
    recompute.resume_policy = "recompute"
    recompute.start_from = None
    recompute.recompute_from = "s1"
    recompute.run(StubRecord(), output, stop_after="s1")

    assert s1_calls == ["product_major", "product_major"]


def test_selective_refresh_recomputes_s0_s2_but_reuses_s1(
    monkeypatch: MonkeyPatch, tmp_path: Path
):
    variants = [{
        "variant_id": "product_major",
        "directory_name": "product_major",
        "branch_product_smiles": "C1CCC1",
    }]
    monkeypatch.setattr(
        v4_orchestrator.V4Orchestrator,
        "_write_s0_from_record",
        staticmethod(_fake_s0_writer(variants)),
    )
    s1_calls, s2_calls, _structures = _install_variant_fakes(monkeypatch)
    output = tmp_path / "run"
    _orchestrator().run(StubRecord(), output, stop_after="s2")

    refreshed = _orchestrator()
    refreshed.refresh_stages = {"s0", "s2"}
    refreshed.run(StubRecord(), output, stop_after="s2")

    assert s1_calls == ["product_major"]
    assert s2_calls == ["product_major", "product_major"]


def test_s1_scheduling_changes_do_not_change_scientific_signature():
    first = _orchestrator()
    first.config["step1"]["censo_lite"] = {
        "b97_3c": {"parallel_jobs": 1, "cores_per_job": 16},
        "xtb_thermo": {"gfn_level": 2},
    }
    second = _orchestrator()
    second.config["step1"]["censo_lite"] = {
        "b97_3c": {"parallel_jobs": 16, "cores_per_job": 1},
        "xtb_thermo": {"gfn_level": 2},
    }

    assert first._s1_signature_payload("C=C", "product_major", "product_major") == second._s1_signature_payload(
        "C=C", "product_major", "product_major"
    )
