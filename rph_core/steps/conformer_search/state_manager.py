import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional


logger = logging.getLogger(__name__)


class ConformerStateManager:

    STATE_FILE = "conformer_state.json"

    def __init__(self, molecule_dir: Path, molecule_name: str):
        self.molecule_dir = Path(molecule_dir)
        self.molecule_dir.mkdir(parents=True, exist_ok=True)
        self.molecule_name = molecule_name
        self.state_path = self.molecule_dir / self.STATE_FILE
        self.state = self._load_or_create()

    def _now(self) -> str:
        return datetime.now().isoformat()

    def _base_state(self) -> Dict[str, Any]:
        now = self._now()
        return {
            "version": "1.0",
            "molecule_name": self.molecule_name,
            "smiles": "",
            "created_at": now,
            "updated_at": now,
            "run": {
                "status": "initialized",
                "two_stage_enabled": False,
                "protocol": "",
                "funnel_signature": {},
                "handoff_signature": {},
                "stages_executed": [],
                "last_error": "",
            },
            "crest": {},
            "conformers": {},
            "summary": {
                "total_conformers": 0,
                "completed": 0,
                "failed": 0,
                "running": 0,
                "best_conformer": "",
                "global_min_energy": None,
                "global_min_xyz": "",
                "ranking_basis": "",
                "candidate_total": 0,
                "selected_candidate_count": 0,
                "fallback_triggered": False,
                "fallback_trigger": False,
                "window_count": 0,
                "gap_rank1_rank2": None,
                "mode_effective": "",
            },
        }

    def _load_or_create(self) -> Dict[str, Any]:
        if self.state_path.exists():
            try:
                with open(self.state_path, "r", encoding="utf-8") as handle:
                    data = json.load(handle)
                if isinstance(data, dict):
                    data.setdefault("crest", {})
                    data.setdefault("conformers", {})
                    data.setdefault("summary", {})
                    data.setdefault("run", {})
                    return data
            except (OSError, json.JSONDecodeError) as exc:
                logger.warning("Failed to load conformer state %s: %s", self.state_path, exc)

        data = self._base_state()
        self.state = data
        self.save()
        return data

    def save(self) -> None:
        self.state["updated_at"] = self._now()
        with open(self.state_path, "w", encoding="utf-8") as handle:
            json.dump(self.state, handle, indent=2, ensure_ascii=False)

    def start_run(self, smiles: str, two_stage_enabled: bool) -> None:
        run = self.state.setdefault("run", {})
        run["status"] = "running"
        run["two_stage_enabled"] = bool(two_stage_enabled)
        run["last_error"] = ""
        self.state["smiles"] = smiles
        self.save()

    def set_protocol_signature(
        self,
        *,
        protocol: str,
        funnel_signature: Dict[str, Any],
        handoff_signature: Dict[str, Any],
        stages_executed: Optional[List[str]] = None,
    ) -> None:
        run = self.state.setdefault("run", {})
        run["protocol"] = str(protocol)
        run["funnel_signature"] = dict(funnel_signature)
        run["handoff_signature"] = dict(handoff_signature)
        if stages_executed is not None:
            run["stages_executed"] = list(stages_executed)
        self.save()

    def is_resume_compatible(
        self,
        *,
        protocol: str,
        funnel_signature: Dict[str, Any],
        handoff_signature: Dict[str, Any],
    ) -> bool:
        run = self.state.get("run", {})
        if not isinstance(run, dict):
            return False
        return (
            run.get("protocol") == str(protocol)
            and run.get("funnel_signature") == dict(funnel_signature)
            and run.get("handoff_signature") == dict(handoff_signature)
        )

    def mark_run_complete(self) -> None:
        self.state.setdefault("run", {})["status"] = "completed"
        self.save()

    def mark_run_failed(self, error_message: str) -> None:
        run = self.state.setdefault("run", {})
        run["status"] = "failed"
        run["last_error"] = str(error_message)
        self.save()

    def mark_crest_stage(
        self,
        stage_name: str,
        status: str,
        output_file: Optional[Path] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        crest = self.state.setdefault("crest", {})
        payload: Dict[str, Any] = {
            "status": status,
            "updated_at": self._now(),
        }
        if output_file is not None:
            try:
                payload["output"] = str(Path(output_file).resolve())
            except OSError:
                payload["output"] = str(output_file)
        if metadata:
            payload.update(metadata)
        crest[stage_name] = payload
        self.save()

    def get_crest_output(self, stage_name: str) -> Optional[Path]:
        stage = self.state.get("crest", {}).get(stage_name, {})
        output = stage.get("output")
        if not output:
            return None
        return Path(output)

    def upsert_conformer(self, conf_name: str, source_xyz: Path, source_index: int) -> None:
        conformers = self.state.setdefault("conformers", {})
        entry = conformers.setdefault(conf_name, {})
        entry.setdefault("name", conf_name)
        entry["source_xyz"] = str(source_xyz)
        entry["source_index"] = source_index
        entry.setdefault("attempts", [])
        entry.setdefault("status", "pending")
        entry.setdefault("record", None)
        self._refresh_summary()
        self.save()

    def mark_conformer_running(self, conf_name: str) -> None:
        entry = self._get_conf(conf_name)
        entry["status"] = "running"
        entry["started_at"] = self._now()
        self._refresh_summary()
        self.save()

    def mark_opt_attempt(self, conf_name: str, attempt: int, status: str, log_file: Path, note: str = "") -> None:
        entry = self._get_conf(conf_name)
        attempts = entry.setdefault("attempts", [])
        attempts.append(
            {
                "stage": "opt",
                "attempt": int(attempt),
                "status": status,
                "log_file": str(log_file),
                "note": note,
                "timestamp": self._now(),
            }
        )
        self.save()

    def mark_sp_result(
        self,
        conf_name: str,
        status: str,
        output_file: Path,
        sp_energy: Optional[float] = None,
        note: str = "",
    ) -> None:
        entry = self._get_conf(conf_name)
        attempts = entry.setdefault("attempts", [])
        payload: Dict[str, Any] = {
            "stage": "sp",
            "status": status,
            "output_file": str(output_file),
            "note": note,
            "timestamp": self._now(),
        }
        if sp_energy is not None:
            payload["sp_energy"] = float(sp_energy)
        attempts.append(payload)
        self.save()

    def mark_conformer_completed(self, conf_name: str, record: Dict[str, Any]) -> None:
        entry = self._get_conf(conf_name)
        entry["status"] = "completed"
        entry["record"] = record
        entry["completed_at"] = self._now()
        self._refresh_summary()
        self.save()

    def mark_conformer_failed(self, conf_name: str, note: str) -> None:
        entry = self._get_conf(conf_name)
        entry["status"] = "failed"
        entry["last_error"] = note
        entry["completed_at"] = self._now()
        self._refresh_summary()
        self.save()

    def is_conformer_completed(self, conf_name: str) -> bool:
        entry = self.state.get("conformers", {}).get(conf_name)
        if not isinstance(entry, dict):
            return False
        if entry.get("status") != "completed":
            return False
        record = entry.get("record")
        return isinstance(record, dict) and record.get("sp_energy") is not None

    def get_conformer_record(self, conf_name: str) -> Optional[Dict[str, Any]]:
        entry = self.state.get("conformers", {}).get(conf_name)
        if not isinstance(entry, dict):
            return None
        record = entry.get("record")
        if not isinstance(record, dict):
            return None
        return record

    def set_global_min(self, best_conf_name: str, energy: float, global_min_xyz: Path) -> None:
        summary = self.state.setdefault("summary", {})
        summary["best_conformer"] = best_conf_name
        summary["global_min_energy"] = float(energy)
        summary["global_min_xyz"] = str(global_min_xyz)
        self.save()

    def set_funnel_summary(
        self,
        *,
        ranking_basis: str,
        candidate_total: int,
        selected_candidate_count: int,
        fallback_triggered: bool,
        window_count: Optional[int] = None,
        gap_rank1_rank2: Optional[float] = None,
        mode_effective: str = "",
        stages_executed: Optional[List[str]] = None,
    ) -> None:
        summary = self.state.setdefault("summary", {})
        summary["ranking_basis"] = str(ranking_basis)
        summary["candidate_total"] = int(candidate_total)
        summary["selected_candidate_count"] = int(selected_candidate_count)
        summary["fallback_triggered"] = bool(fallback_triggered)
        summary["fallback_trigger"] = bool(fallback_triggered)
        if window_count is not None:
            summary["window_count"] = int(window_count)
        summary["gap_rank1_rank2"] = gap_rank1_rank2
        summary["mode_effective"] = str(mode_effective)
        if stages_executed is not None:
            self.state.setdefault("run", {})["stages_executed"] = list(stages_executed)
        self.save()

    def get_summary(self) -> Dict[str, Any]:
        summary = self.state.get("summary", {})
        if isinstance(summary, dict):
            return dict(summary)
        return {}

    def _get_conf(self, conf_name: str) -> Dict[str, Any]:
        conformers = self.state.setdefault("conformers", {})
        if conf_name not in conformers:
            conformers[conf_name] = {
                "name": conf_name,
                "status": "pending",
                "attempts": [],
                "record": None,
            }
        entry = conformers[conf_name]
        if not isinstance(entry, dict):
            conformers[conf_name] = {
                "name": conf_name,
                "status": "pending",
                "attempts": [],
                "record": None,
            }
            entry = conformers[conf_name]
        return entry

    def _refresh_summary(self) -> None:
        conformers = self.state.get("conformers", {})
        if not isinstance(conformers, dict):
            return
        total = len(conformers)
        completed = 0
        failed = 0
        running = 0
        for value in conformers.values():
            if not isinstance(value, dict):
                continue
            status = value.get("status")
            if status == "completed":
                completed += 1
            elif status == "failed":
                failed += 1
            elif status == "running":
                running += 1
        summary = self.state.setdefault("summary", {})
        summary["total_conformers"] = total
        summary["completed"] = completed
        summary["failed"] = failed
        summary["running"] = running
