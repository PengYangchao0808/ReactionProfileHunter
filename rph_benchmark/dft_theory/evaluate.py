#!/usr/bin/env python3
# pyright: reportDeprecated=false, reportUnusedCallResult=false, reportAny=false, reportExplicitAny=false, reportUnknownMemberType=false, reportUnknownVariableType=false, reportUnknownArgumentType=false, reportUnknownLambdaType=false, reportReturnType=false
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from benchmark.dft_theory.lib.state_manager import read_benchmark_manifest, read_rx_manifest, write_benchmark_manifest
from benchmark.dft_theory.lib.state_manager import require_phase1_canonical_results, require_phase2_canonical_results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Aggregate final DFT benchmark report.")
    parser.add_argument("--session-dir", required=True, help="Benchmark session directory")
    parser.add_argument(
        "--phase",
        choices=["final"],
        required=True,
        help="Evaluation phase",
    )
    parser.add_argument("--rx-id", nargs="*", default=None, help="Optional subset of rx ids")
    return parser.parse_args()


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def resolve_rx_ids(session_dir: Path, requested_rx: Optional[Sequence[str]]) -> List[str]:
    benchmark_manifest = read_benchmark_manifest(session_dir)
    configured = [str(item) for item in benchmark_manifest.get("cases", []) if str(item).strip()]
    if not requested_rx:
        return configured
    requested = {str(item).strip() for item in requested_rx if str(item).strip()}
    return [rx_id for rx_id in configured if rx_id in requested]


def _winner_from_payload(payload: Any) -> Optional[str]:
    if isinstance(payload, str) and payload.strip():
        return payload.strip()
    if not isinstance(payload, Mapping):
        return None
    for key in ("method_id", "winner_method_id", "id", "method", "winner"):
        value = payload.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()
    return None


def resolve_phase_winner(rx_manifest: Mapping[str, Any], phase_key: str, current_key: str) -> Optional[str]:
    phase_payload = rx_manifest.get(phase_key)
    if isinstance(phase_payload, dict):
        if "winner_artifact" in phase_payload:
            raise ValueError(f"Legacy {phase_key} winner_artifact is not supported by benchmark v2")
        for key in ("winner", "winner_method_id"):
            winner = _winner_from_payload(phase_payload.get(key))
            if winner:
                return winner
    current = rx_manifest.get(current_key)
    if isinstance(current, str) and current.strip():
        return current.strip()
    return None


def tally_winner(votes: Dict[str, int], winner: Optional[str]) -> None:
    if winner:
        votes[winner] = votes.get(winner, 0) + 1


def render_final_markdown(final_report: Mapping[str, Any]) -> str:
    phase1_votes = final_report.get("phase1_winner_votes", {})
    phase2_votes = final_report.get("phase2_winner_votes", {})
    lines = [
        "# Final DFT Benchmark Report",
        "",
        f"- generated_at: {final_report.get('generated_at')}",
        f"- recommended_l2_sp_method: {final_report.get('recommended_l2_sp_method') or 'NA'}",
        f"- recommended_l1_geo_method: {final_report.get('recommended_l1_geo_method') or 'NA'}",
        "",
        "## Phase1 Winner Votes",
        "",
        "| method | votes |",
        "|---|---:|",
    ]
    if isinstance(phase1_votes, dict):
        for method, votes in sorted(phase1_votes.items(), key=lambda item: item[1], reverse=True):
            lines.append(f"| {method} | {votes} |")
    lines.extend(["", "## Phase2 Winner Votes", "", "| method | votes |", "|---|---:|"])
    if isinstance(phase2_votes, dict):
        for method, votes in sorted(phase2_votes.items(), key=lambda item: item[1], reverse=True):
            lines.append(f"| {method} | {votes} |")
    return "\n".join(lines) + "\n"


def evaluate_phase_final(session_dir: Path, rx_ids: Sequence[str]) -> int:
    benchmark_manifest = read_benchmark_manifest(session_dir)

    phase1_winner_votes: Dict[str, int] = {}
    phase2_winner_votes: Dict[str, int] = {}
    for rx_id in rx_ids:
        rx_manifest = read_rx_manifest(session_dir, rx_id)
        phase1_results = rx_manifest.get("phase1", {}).get("results", {})
        if phase1_results:
            if not isinstance(phase1_results, Mapping):
                raise ValueError(f"Final evaluation requires canonical phase1 results for rx{rx_id}")
            for method_id, row in phase1_results.items():
                if isinstance(row, Mapping):
                    require_phase1_canonical_results(row, context=f"Final evaluation phase1 cache rx{rx_id}/{method_id}")

        phase2_results = rx_manifest.get("phase2", {}).get("results", {})
        if phase2_results:
            if not isinstance(phase2_results, Mapping):
                raise ValueError(f"Final evaluation requires canonical phase2 results for rx{rx_id}")
            for method_id, row in phase2_results.items():
                if isinstance(row, Mapping):
                    require_phase2_canonical_results(row, context=f"Final evaluation phase2 cache rx{rx_id}/{method_id}")
        tally_winner(phase1_winner_votes, resolve_phase_winner(rx_manifest, "phase1", "current_phase1_winner"))
        tally_winner(phase2_winner_votes, resolve_phase_winner(rx_manifest, "phase2", "current_phase2_winner"))

    recommended_l2_sp_method = (
        max(phase1_winner_votes.items(), key=lambda item: item[1])[0] if phase1_winner_votes else None
    )
    recommended_l1_geo_method = (
        max(phase2_winner_votes.items(), key=lambda item: item[1])[0] if phase2_winner_votes else None
    )

    final_report = {
        "generated_at": now_iso(),
        "session_dir": str(session_dir),
        "cases": [f"rx{rx_id}" for rx_id in rx_ids],
        "phase1_winner_votes": phase1_winner_votes,
        "phase2_winner_votes": phase2_winner_votes,
        "recommended_l2_sp_method_id": recommended_l2_sp_method,
        "recommended_l2_sp_method": recommended_l2_sp_method,
        "recommended_l1_geo_method": recommended_l1_geo_method,
    }

    final_report_json = (session_dir / "reports" / "final_report.json").resolve()
    final_report_md = (session_dir / "reports" / "final_report.md").resolve()
    write_json(final_report_json, final_report)
    write_text(final_report_md, render_final_markdown(final_report))

    benchmark_manifest.setdefault("stages", {})
    benchmark_manifest["stages"]["final_report"] = "done"
    benchmark_manifest.setdefault("reports", {})
    benchmark_manifest["reports"]["final_report_json"] = str(final_report_json)
    benchmark_manifest["reports"]["final_report_md"] = str(final_report_md)
    write_benchmark_manifest(session_dir, benchmark_manifest)
    return 0


def main() -> int:
    args = parse_args()
    session_dir = Path(str(args.session_dir)).expanduser().resolve()
    if not session_dir.exists():
        raise FileNotFoundError(f"Session dir does not exist: {session_dir}")

    rx_ids = resolve_rx_ids(session_dir, args.rx_id)
    if not rx_ids:
        raise RuntimeError("No rx ids resolved for evaluation")

    return evaluate_phase_final(session_dir, rx_ids)


if __name__ == "__main__":
    raise SystemExit(main())
