import json
from pathlib import Path

from scripts.parse_s1_progress import aggregate_progress, parse_events


def test_parse_and_aggregate_s1_progress(tmp_path: Path) -> None:
    log_path = tmp_path / "rph.log"
    lines = [
        "2026-03-13 10:00:00 - rph_core.steps - INFO - S1_PROGRESS|"
        + json.dumps(
            {
                "schema": "s1_progress_v1",
                "event": "run_started",
                "molecule": "product",
                "state_file": "/tmp/product/conformer_state.json",
                "summary": {"total_conformers": 2, "completed": 0, "failed": 0, "running": 0},
                "payload": {"smiles": "C=C"},
            }
        ),
        "2026-03-13 10:01:00 - rph_core.steps - INFO - S1_PROGRESS|"
        + json.dumps(
            {
                "schema": "s1_progress_v1",
                "event": "conformer_started",
                "molecule": "product",
                "state_file": "/tmp/product/conformer_state.json",
                "summary": {"total_conformers": 2, "completed": 0, "failed": 0, "running": 1},
                "payload": {"conformer": "conf_000", "attempt": 0},
            }
        ),
        "2026-03-13 10:02:00 - rph_core.steps - INFO - S1_PROGRESS|"
        + json.dumps(
            {
                "schema": "s1_progress_v1",
                "event": "conformer_completed",
                "molecule": "product",
                "state_file": "/tmp/product/conformer_state.json",
                "summary": {"total_conformers": 2, "completed": 1, "failed": 0, "running": 0},
                "payload": {"conformer": "conf_000", "energy_hartree": -10.123},
            }
        ),
        "2026-03-13 10:02:10 - rph_core.steps - INFO - S1_PROGRESS|"
        + json.dumps(
            {
                "schema": "s1_progress_v1",
                "event": "best_conformer_selected",
                "molecule": "product",
                "state_file": "/tmp/product/conformer_state.json",
                "summary": {"total_conformers": 2, "completed": 1, "failed": 0, "running": 0},
                "payload": {"conformer": "conf_000", "weight": 0.8, "energy_hartree": -10.123},
            }
        ),
        "2026-03-13 10:03:00 - rph_core.steps - INFO - S1_PROGRESS|"
        + json.dumps(
            {
                "schema": "s1_progress_v1",
                "event": "run_completed",
                "molecule": "product",
                "state_file": "/tmp/product/conformer_state.json",
                "summary": {"total_conformers": 2, "completed": 1, "failed": 0, "running": 0},
                "payload": {"global_min_xyz": "/tmp/product/product_global_min.xyz", "energy_hartree": -10.123},
            }
        ),
        "2026-03-13 10:03:30 - rph_core.steps - INFO - S1_ANCHOR_PROGRESS|"
        + json.dumps(
            {
                "schema": "s1_anchor_progress_v1",
                "event": "molecule_completed",
                "payload": {"molecule": "product", "index": 1, "total": 1, "energy_hartree": -10.123},
            }
        ),
        "2026-03-13 10:04:00 - rph_core.steps - INFO - S1_ANCHOR_PROGRESS|"
        + json.dumps(
            {
                "schema": "s1_anchor_progress_v1",
                "event": "anchor_finished",
                "payload": {"status": "completed", "total": 1, "completed": 1, "failed": 0},
            }
        ),
    ]
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    s1_events, anchor_events = parse_events(log_path)
    aggregate = aggregate_progress(s1_events, anchor_events, log_file=log_path)

    assert aggregate["counts"]["s1_progress_events"] == 5
    assert aggregate["counts"]["anchor_events"] == 2
    assert aggregate["counts"]["molecules"] == 1
    assert aggregate["molecules"]["product"]["run_status"] == "completed"
    assert aggregate["molecules"]["product"]["conformers"]["conf_000"]["status"] == "completed"
    assert aggregate["molecules"]["product"]["conformers"]["conf_000"]["best_selected"] is True
    assert aggregate["anchor"]["status"] == "completed"


def test_parse_events_ignores_invalid_json_payload(tmp_path: Path) -> None:
    log_path = tmp_path / "rph.log"
    log_path.write_text(
        "2026-03-13 - INFO - S1_PROGRESS|{bad json}\n"
        "2026-03-13 - INFO - S1_PROGRESS|[]\n"
        "2026-03-13 - INFO - normal line\n",
        encoding="utf-8",
    )

    s1_events, anchor_events = parse_events(log_path)
    assert s1_events == []
    assert anchor_events == []


def test_conformer_skipped_marked_cached(tmp_path: Path) -> None:
    log_path = tmp_path / "rph.log"
    log_path.write_text(
        "2026-03-13 - INFO - S1_PROGRESS|"
        + json.dumps(
            {
                "schema": "s1_progress_v1",
                "event": "conformer_skipped",
                "molecule": "product",
                "state_file": "/tmp/product/conformer_state.json",
                "summary": {"total_conformers": 1, "completed": 1, "failed": 0, "running": 0},
                "payload": {"conformer": "conf_000"},
            }
        )
        + "\n",
        encoding="utf-8",
    )

    s1_events, anchor_events = parse_events(log_path)
    aggregate = aggregate_progress(s1_events, anchor_events, log_file=log_path)
    assert aggregate["molecules"]["product"]["conformers"]["conf_000"]["status"] == "cached"
