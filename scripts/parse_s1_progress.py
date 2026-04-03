import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


S1_PROGRESS_PREFIX = "S1_PROGRESS|"
S1_ANCHOR_PROGRESS_PREFIX = "S1_ANCHOR_PROGRESS|"


def _extract_prefixed_json(line: str, prefix: str) -> Optional[Dict[str, Any]]:
    marker_idx = line.find(prefix)
    if marker_idx < 0:
        return None

    payload = line[marker_idx + len(prefix):].strip()
    if not payload:
        return None

    try:
        parsed = json.loads(payload)
    except json.JSONDecodeError:
        return None

    if not isinstance(parsed, dict):
        return None
    return parsed


def parse_events(log_file: Path) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    s1_events: List[Dict[str, Any]] = []
    anchor_events: List[Dict[str, Any]] = []

    with open(log_file, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            s1_payload = _extract_prefixed_json(line, S1_PROGRESS_PREFIX)
            if s1_payload is not None:
                s1_events.append(s1_payload)
                continue

            anchor_payload = _extract_prefixed_json(line, S1_ANCHOR_PROGRESS_PREFIX)
            if anchor_payload is not None:
                anchor_events.append(anchor_payload)

    return s1_events, anchor_events


def _ensure_molecule(bucket: Dict[str, Any], molecule: str, state_file: str) -> Dict[str, Any]:
    molecules = bucket.setdefault("molecules", {})
    if molecule not in molecules:
        molecules[molecule] = {
            "event_count": 0,
            "state_file": state_file,
            "run_status": "unknown",
            "last_event": "",
            "summary": {},
            "conformers": {},
        }
    return molecules[molecule]


def _map_conformer_status(event_name: str) -> str:
    if event_name in {"conformer_completed", "sp_completed"}:
        return "completed"
    if event_name in {"conformer_failed", "sp_failed", "opt_failed"}:
        return "failed"
    if event_name in {"conformer_started", "opt_converged"}:
        return "running"
    if event_name == "conformer_skipped":
        return "cached"
    return "unknown"


def aggregate_progress(
    s1_events: List[Dict[str, Any]],
    anchor_events: List[Dict[str, Any]],
    log_file: Optional[Path] = None,
) -> Dict[str, Any]:
    result: Dict[str, Any] = {
        "schema": "s1_progress_aggregate_v1",
        "log_file": "" if log_file is None else str(log_file),
        "counts": {
            "s1_progress_events": len(s1_events),
            "anchor_events": len(anchor_events),
        },
        "anchor": {
            "status": "unknown",
            "summary": {},
            "molecules": {},
        },
        "molecules": {},
    }

    for event in s1_events:
        molecule = str(event.get("molecule", ""))
        if not molecule:
            continue

        state_file = str(event.get("state_file", ""))
        entry = _ensure_molecule(result, molecule, state_file)
        entry["event_count"] = int(entry.get("event_count", 0)) + 1
        entry["last_event"] = str(event.get("event", ""))
        entry["summary"] = event.get("summary", {}) if isinstance(event.get("summary"), dict) else {}
        if state_file:
            entry["state_file"] = state_file

        event_name = str(event.get("event", ""))
        if event_name == "run_started":
            entry["run_status"] = "running"
        elif event_name in {"run_completed", "run_reused"}:
            entry["run_status"] = "completed"
        elif event_name == "run_failed":
            entry["run_status"] = "failed"

        payload = event.get("payload", {})
        if not isinstance(payload, dict):
            continue

        conf_name = payload.get("conformer")
        if conf_name is None:
            continue

        conf_key = str(conf_name)
        conformers = entry.setdefault("conformers", {})
        conf_entry = conformers.setdefault(
            conf_key,
            {
                "status": "unknown",
                "last_event": "",
                "attempt": None,
                "energy_hartree": None,
                "best_selected": False,
            },
        )
        mapped_status = _map_conformer_status(event_name)
        if mapped_status != "unknown":
            conf_entry["status"] = mapped_status
        conf_entry["last_event"] = event_name
        if "attempt" in payload:
            conf_entry["attempt"] = payload.get("attempt")
        if "energy_hartree" in payload:
            conf_entry["energy_hartree"] = payload.get("energy_hartree")
        if event_name == "best_conformer_selected":
            conf_entry["best_selected"] = True
            if "weight" in payload:
                conf_entry["weight"] = payload.get("weight")

    anchor_status = result["anchor"]
    for event in anchor_events:
        event_name = str(event.get("event", ""))
        payload = event.get("payload", {})
        if not isinstance(payload, dict):
            payload = {}

        if event_name == "anchor_finished":
            anchor_status["status"] = str(payload.get("status", "unknown"))
            anchor_status["summary"] = payload

        molecule = payload.get("molecule")
        if molecule is not None:
            key = str(molecule)
            anchor_molecules = anchor_status.setdefault("molecules", {})
            current = anchor_molecules.setdefault(key, {})
            current.update(payload)
            if event_name == "molecule_started":
                current["status"] = "running"
            elif event_name == "molecule_completed":
                current["status"] = "completed"
            elif event_name == "molecule_failed":
                current["status"] = "failed"

    result["counts"]["molecules"] = len(result.get("molecules", {}))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate S1 structured progress events from rph.log")
    parser.add_argument("--log", required=True, help="Path to rph.log")
    parser.add_argument("--output", help="Optional output JSON path")
    args = parser.parse_args()

    log_path = Path(args.log)
    if not log_path.exists():
        raise FileNotFoundError(f"Log file not found: {log_path}")

    s1_events, anchor_events = parse_events(log_path)
    aggregate = aggregate_progress(s1_events, anchor_events, log_file=log_path)

    text = json.dumps(aggregate, indent=2, ensure_ascii=False, sort_keys=True)
    if args.output:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
