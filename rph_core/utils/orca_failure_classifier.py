"""Structured ORCA failure classification helpers."""

from __future__ import annotations

from dataclasses import dataclass
import enum
from functools import lru_cache
import logging
from pathlib import Path
import re
from typing import Any, Mapping, Optional

import yaml

logger = logging.getLogger(__name__)

_RE_FLAGS = re.IGNORECASE | re.MULTILINE
_DEFAULT_ENABLED = True
_DEFAULT_SCAN_TAIL_CHARS = 50000
_DEFAULT_EVIDENCE_CONTEXT_LINES = 1

_NORMAL_TERMINATION_PATTERN = r"ORCA TERMINATED NORMALLY"
_SIGNAL_RUNTIME_PATTERN = r"(?:SIGSEGV|Segmentation fault|forrtl: severe)"
_SIGNAL_MPI_PATTERN = r"(?:MPI|mpi|pml|ob1)"
_SIGNAL_MEMORY_PATTERN = r"(?:out of memory|Cannot allocate memory|MemoryError|mmap failed)"
_SCF_SUBPROCESS_QCMSG_PATTERN = r"qcmsg\.cpp:\s*\d+"
# "aborting" alone is not enough: ORCA prints the NON-fatal "* ABORTING THE RUN *"
# SOSCF banner and then recovers via TRAH. The negative lookahead keeps
# "error termination in SCF gradient" owned by SCF_GRADIENT_ABORTED.
_SCF_SUBPROCESS_ABORTING_PATTERN = r"error\s+termination\s+in\s+SCF(?!\s+gradient)"
_SCF_NOT_CONVERGED_PATTERN = r"SCF NOT CONVERGED"
_SCF_CONVERGENCE_HALTED_PATTERN = r"SCF convergence halted"
_SCF_NOT_CONVERGED_LOWER_PATTERN = r"scf not converged"
_SCF_GRADIENT_ABORT_PATTERN = r"(?:error\s+termination\s+in\s+SCF\s+gradient|\bgradient\b.{0,40}\babort)"
_SCF_GRADIENT_ABNORMAL_PATTERN = r"ORCA finished computing gradients abnormally"
_MPI_RUNTIME_PATTERN = r"MPI.{0,80}(?:abort|fail|error)"
_OPENMPI_RUNTIME_PATTERN = r"(?:ompi|openmpi|mpiexec).{0,80}(?:error|abort)"
_MEMORY_PATTERN = r"(?:out of memory|cannot allocate memory|mmap failed|memoryerror)"
_MP2_ABORT_PATTERN = r"mp2.{0,40}(?:abort|fail)"
_RI_ABORT_PATTERN = r"resolution-of-identity.{0,40}(?:abort|fail)"
_INTERRUPTED_ABNORMAL_PATTERN = r"ORCA terminated abnormally"
_INTERRUPTED_ERROR_PATTERN = r"ORCA finished with error"
# ORCA 5/6 can wrap the final MaxIter warning across lines, e.g.::
#
#   The optimization did not converge but reached the maximum number of
#   optimization cycles.
#
# Keep this distinct from the per-cycle "has not yet converged" progress
# message: only the terminal MaxIter wording is a valid rescue trigger.
_OPT_NOT_CONVERGED_PATTERN = (
    r"THE\s+OPTIMIZATION\s+DID\s+NOT\s+CONVERGE"
    r"(?:\s+BUT\s+REACHED\s+THE\s+MAXIMUM\s+NUMBER\s+OF"
    r"\s+OPTIMIZATION\s+CYCLES)?"
)
_OPT_MAX_CYCLES_PATTERN = r"maximum\s+number\s+of\s+optimization\s+cycles"
# ORCA's GFN-xTB looks for xtb only inside the ORCA binary directory (never PATH).
_XTB_BINARY_UNAVAILABLE_PATTERN = (
    r"(?:otool_xtb\s+not\s+available|"
    r"Please\s+provide\s+the\s+xtb\s+executables\s+in\s+the\s+same\s+path)"
)
# ORCA 6.x failure signatures (verified against official 6.1 manual and a real
# 6.0.1 crash log).  Keep the qcmsg-based SCF patterns owned by the SCF types.
_ORCA_ERROR_TERMINATION_PATTERN = r"ORCA\s+finished\s+by\s+error\s+termination\s+in\s+(?!SCF\b)\S+"
_QCMSG_ABORT_PATTERN = r"\[file\s+orca_tools/qcmsg\.cpp,\s*line\s+\d+\].*?aborting\s+the\s+run"
_CALLING_COMMAND_PATTERN = r"Calling\s+Command:\s+mpirun"
_GBW_VERSION_MISMATCH_PATTERN = (
    r"Your\s+GBWFile\s+is\s+either\s+corrupt\s+or\s+from\s+a\s+different\s+ORCA\s+version"
)
_INPUT_NON_ASCII_PATTERN = r"non-ASCII\s+character\(s\)\s+found\s+in\s+words\s+on\s+line"
_EXTERNAL_KILLED_PATTERN = (
    r"(?:slurmstepd:\s*error:.*CANCELLED|srun:\s*forcing\s+job\s+termination|"
    r"Job\s+terminated\s+from\s+outer|\bmake:\s*\*\*\*.*\bTerminated\b)"
)
_THREAD_SPAWN_FAILED_PATTERN = (
    r"libgomp:\s*Thread\s+creation\s+failed:\s*Resource\s+temporarily\s+unavailable"
)
_MPI_RANK_KILLED_PATTERN = r"mpirun\s+noticed\s+that\s+process\s+rank"
_MPI_PRIMARY_ABORT_PATTERN = (
    r"Primary\s+job\s+terminated\s+normally,\s*but\s+\d+\s+process(?:es)?\s+returned"
)
_MPI_PROCS_KILLED_PATTERN = r"total\s+processes\s+killed"
_OOM_PROCESS_PATTERN = r"\[file\s+orca_tools/qcmem\.cpp.*Process\s+\d+\].*OUT\s+OF\s+MEMORY"
_MAXCORE_PATTERN = r"Please\s+increase\s+MaxCore"

_NORMAL_TERMINATION_RE = re.compile(_NORMAL_TERMINATION_PATTERN, _RE_FLAGS)
_SIGNAL_RUNTIME_RE = re.compile(_SIGNAL_RUNTIME_PATTERN, _RE_FLAGS)
_SIGNAL_MPI_RE = re.compile(_SIGNAL_MPI_PATTERN, _RE_FLAGS)
_SIGNAL_MEMORY_RE = re.compile(_SIGNAL_MEMORY_PATTERN, _RE_FLAGS)
_SCF_SUBPROCESS_QCMSG_RE = re.compile(_SCF_SUBPROCESS_QCMSG_PATTERN, _RE_FLAGS)
_SCF_SUBPROCESS_ABORTING_RE = re.compile(_SCF_SUBPROCESS_ABORTING_PATTERN, _RE_FLAGS)
_SCF_NOT_CONVERGED_RE = re.compile(_SCF_NOT_CONVERGED_PATTERN, _RE_FLAGS)
_SCF_CONVERGENCE_HALTED_RE = re.compile(_SCF_CONVERGENCE_HALTED_PATTERN, _RE_FLAGS)
_SCF_NOT_CONVERGED_LOWER_RE = re.compile(_SCF_NOT_CONVERGED_LOWER_PATTERN, _RE_FLAGS)
_SCF_GRADIENT_ABORT_RE = re.compile(_SCF_GRADIENT_ABORT_PATTERN, _RE_FLAGS)
_SCF_GRADIENT_ABNORMAL_RE = re.compile(_SCF_GRADIENT_ABNORMAL_PATTERN, _RE_FLAGS)
_MPI_RUNTIME_RE = re.compile(_MPI_RUNTIME_PATTERN, _RE_FLAGS)
_OPENMPI_RUNTIME_RE = re.compile(_OPENMPI_RUNTIME_PATTERN, _RE_FLAGS)
_MEMORY_RE = re.compile(_MEMORY_PATTERN, _RE_FLAGS)
_MP2_ABORT_RE = re.compile(_MP2_ABORT_PATTERN, _RE_FLAGS)
_RI_ABORT_RE = re.compile(_RI_ABORT_PATTERN, _RE_FLAGS)
_INTERRUPTED_ABNORMAL_RE = re.compile(_INTERRUPTED_ABNORMAL_PATTERN, _RE_FLAGS)
_INTERRUPTED_ERROR_RE = re.compile(_INTERRUPTED_ERROR_PATTERN, _RE_FLAGS)
_OPT_NOT_CONVERGED_RE = re.compile(_OPT_NOT_CONVERGED_PATTERN, _RE_FLAGS)
_OPT_MAX_CYCLES_RE = re.compile(_OPT_MAX_CYCLES_PATTERN, _RE_FLAGS)
_XTB_BINARY_UNAVAILABLE_RE = re.compile(_XTB_BINARY_UNAVAILABLE_PATTERN, _RE_FLAGS)
_ORCA_ERROR_TERMINATION_RE = re.compile(_ORCA_ERROR_TERMINATION_PATTERN, _RE_FLAGS)
_QCMSG_ABORT_RE = re.compile(_QCMSG_ABORT_PATTERN, _RE_FLAGS)
_CALLING_COMMAND_RE = re.compile(_CALLING_COMMAND_PATTERN, _RE_FLAGS)
_GBW_VERSION_MISMATCH_RE = re.compile(_GBW_VERSION_MISMATCH_PATTERN, _RE_FLAGS)
_INPUT_NON_ASCII_RE = re.compile(_INPUT_NON_ASCII_PATTERN, _RE_FLAGS)
_EXTERNAL_KILLED_RE = re.compile(_EXTERNAL_KILLED_PATTERN, _RE_FLAGS)
_THREAD_SPAWN_FAILED_RE = re.compile(_THREAD_SPAWN_FAILED_PATTERN, _RE_FLAGS)
_MPI_RANK_KILLED_RE = re.compile(_MPI_RANK_KILLED_PATTERN, _RE_FLAGS)
_MPI_PRIMARY_ABORT_RE = re.compile(_MPI_PRIMARY_ABORT_PATTERN, _RE_FLAGS)
_MPI_PROCS_KILLED_RE = re.compile(_MPI_PROCS_KILLED_PATTERN, _RE_FLAGS)
_OOM_PROCESS_RE = re.compile(_OOM_PROCESS_PATTERN, _RE_FLAGS)
_MAXCORE_RE = re.compile(_MAXCORE_PATTERN, _RE_FLAGS)


class OrcaFailureType(str, enum.Enum):
    SCF_SUBPROCESS_ABORTED = "scf_subprocess_aborted"
    SCF_NOT_CONVERGED = "scf_not_converged"
    SCF_GRADIENT_ABORTED = "scf_gradient_aborted"
    MPI_RUNTIME_FAILED = "mpi_runtime_failed"
    MEMORY_PRESSURE = "memory_pressure"
    TIMEOUT = "timeout"
    GEOMETRY_OPTIMIZATION_NOT_CONVERGED = "geometry_optimization_not_converged"
    INTERRUPTED = "interrupted"
    XTB_BINARY_UNAVAILABLE = "xtb_binary_unavailable"
    ORCA_ERROR_TERMINATION = "orca_error_termination"
    GBW_VERSION_MISMATCH = "gbw_version_mismatch"
    INPUT_NON_ASCII = "input_non_ascii"
    EXTERNAL_KILLED = "external_killed"
    THREAD_SPAWN_FAILED = "thread_spawn_failed"


@dataclass(frozen=True)
class OrcaFailureClassification:
    failure_type: OrcaFailureType
    evidence_lines: tuple[str, ...]
    returncode: int | None
    retry_hint: str
    summary: str


@dataclass(frozen=True)
class OrcaFailureClassifierSettings:
    enabled: bool = _DEFAULT_ENABLED
    scan_tail_chars: int = _DEFAULT_SCAN_TAIL_CHARS
    evidence_context_lines: int = _DEFAULT_EVIDENCE_CONTEXT_LINES


def extract_evidence(output_text: str, pattern: str, context_lines: int = 1) -> tuple[str, ...]:
    """Return matching lines plus surrounding context for a regex pattern."""
    if not output_text:
        return ()

    compiled = re.compile(pattern, _RE_FLAGS)
    lines = output_text.splitlines()
    if not lines:
        return ()

    starts: list[int] = []
    cursor = 0
    for line in lines:
        starts.append(cursor)
        cursor += len(line) + 1

    selected_indexes: list[int] = []
    seen_indexes: set[int] = set()
    for match in compiled.finditer(output_text):
        line_index = 0
        for index, line_start in enumerate(starts):
            next_start = starts[index + 1] if index + 1 < len(starts) else cursor + 1
            if line_start <= match.start() < next_start:
                line_index = index
                break
        start = max(0, line_index - context_lines)
        stop = min(len(lines), line_index + context_lines + 1)
        for idx in range(start, stop):
            if idx not in seen_indexes:
                seen_indexes.add(idx)
                selected_indexes.append(idx)

    return tuple(_trim_line(lines[idx]) for idx in selected_indexes)


def has_normal_termination(output_text: str) -> bool:
    """Return True when ORCA reports normal termination."""
    return bool(_NORMAL_TERMINATION_RE.search(output_text or ""))


def resolve_orca_failure_classifier_settings(
    config: Optional[Mapping[str, Any]] = None,
) -> OrcaFailureClassifierSettings:
    """Resolve classifier settings from runtime config or repo defaults."""
    if config is None:
        return _load_default_settings()

    section = config.get("orca_failure_classification", {}) if isinstance(config, Mapping) else {}
    if not isinstance(section, Mapping):
        section = {}
    defaults = OrcaFailureClassifierSettings()

    enabled = section.get("enabled", defaults.enabled)
    scan_tail_chars = section.get("scan_tail_chars", defaults.scan_tail_chars)
    evidence_context_lines = section.get(
        "evidence_context_lines",
        defaults.evidence_context_lines,
    )
    return OrcaFailureClassifierSettings(
        enabled=bool(enabled),
        scan_tail_chars=_coerce_positive_int(scan_tail_chars, defaults.scan_tail_chars),
        evidence_context_lines=_coerce_nonnegative_int(
            evidence_context_lines,
            defaults.evidence_context_lines,
        ),
    )


def classify_orca_failure(
    *,
    output_text: str,
    stderr_text: str = "",
    returncode: int | None = None,
    timed_out: bool = False,
    max_cycles_opt: int | None = None,
    scf_converged: bool | None = None,
) -> OrcaFailureClassification | None:
    return _classify_orca_failure_impl(
        output_text=output_text,
        stderr_text=stderr_text,
        returncode=returncode,
        timed_out=timed_out,
        max_cycles_opt=max_cycles_opt,
        scf_converged=scf_converged,
        settings=_load_default_settings(),
    )


def classify_orca_failure_for_config(
    *,
    output_text: str,
    stderr_text: str = "",
    returncode: int | None = None,
    timed_out: bool = False,
    max_cycles_opt: int | None = None,
    scf_converged: bool | None = None,
    config: Optional[Mapping[str, Any]] = None,
) -> OrcaFailureClassification | None:
    return _classify_orca_failure_impl(
        output_text=output_text,
        stderr_text=stderr_text,
        returncode=returncode,
        timed_out=timed_out,
        max_cycles_opt=max_cycles_opt,
        scf_converged=scf_converged,
        settings=resolve_orca_failure_classifier_settings(config),
    )


# Failure types that legitimately coexist with a normal-termination banner and a
# zero/None exit code:
#   - GEOMETRY_OPTIMIZATION_NOT_CONVERGED: ORCA prints "****ORCA TERMINATED
#     NORMALLY****" even when the optimization reaches MaxIter.
#   - SCF_NOT_CONVERGED: an SP whose SCF did not converge may still print the
#     normal banner and exit 0 while reporting a (flagged) final energy.
_NORMAL_TERMINATION_COMPATIBLE_FAILURES = frozenset(
    {
        OrcaFailureType.GEOMETRY_OPTIMIZATION_NOT_CONVERGED,
        OrcaFailureType.SCF_NOT_CONVERGED,
    }
)


def _pattern_match_is_false_positive(
    output_text: str,
    returncode: int | None,
    classification: OrcaFailureClassification,
) -> bool:
    """A pattern match on a normally-terminating zero-rc run is a false positive
    unless the matched failure type legitimately coexists with that banner."""
    if not has_normal_termination(output_text or ""):
        return False
    if returncode not in (None, 0):
        return False
    return classification.failure_type not in _NORMAL_TERMINATION_COMPATIBLE_FAILURES


def _classify_orca_failure_impl(
    *,
    output_text: str,
    stderr_text: str,
    returncode: int | None,
    timed_out: bool,
    max_cycles_opt: int | None,
    scf_converged: bool | None,
    settings: OrcaFailureClassifierSettings,
) -> OrcaFailureClassification | None:
    if not settings.enabled:
        return None

    scanned_output = _tail_text(output_text or "", settings.scan_tail_chars)
    scanned_stderr = _tail_text(stderr_text or "", settings.scan_tail_chars)

    if timed_out:
        return _emit_classification(
            OrcaFailureClassification(
                failure_type=OrcaFailureType.TIMEOUT,
                evidence_lines=("process timed out",),
                returncode=returncode,
                retry_hint="increase_timeout",
                summary="ORCA process timed out",
            )
        )

    if returncode is not None and returncode < 0 and abs(returncode) in (6, 9, 11, 15):
        evidence_lines = _signal_evidence(scanned_stderr, returncode)
        if _SIGNAL_MEMORY_RE.search(scanned_stderr):
            return _emit_classification(
                OrcaFailureClassification(
                    failure_type=OrcaFailureType.MEMORY_PRESSURE,
                    evidence_lines=evidence_lines,
                    returncode=returncode,
                    retry_hint="increase_memory",
                    summary="ORCA failed due to memory pressure",
                )
            )
        if _SIGNAL_MPI_RE.search(scanned_stderr) or _SIGNAL_RUNTIME_RE.search(scanned_stderr):
            return _emit_classification(
                OrcaFailureClassification(
                    failure_type=OrcaFailureType.MPI_RUNTIME_FAILED,
                    evidence_lines=evidence_lines,
                    returncode=returncode,
                    retry_hint="reduce_ranks",
                    summary="ORCA MPI runtime failed",
                )
            )
        return _emit_classification(
            OrcaFailureClassification(
                failure_type=OrcaFailureType.MPI_RUNTIME_FAILED,
                evidence_lines=evidence_lines,
                returncode=returncode,
                retry_hint="reduce_ranks",
                summary="ORCA MPI runtime failed",
            )
        )

    pattern_classification = _match_output_patterns(
        scanned_output,
        returncode=returncode,
        max_cycles_opt=max_cycles_opt,
        scf_converged=scf_converged,
        evidence_context_lines=settings.evidence_context_lines,
    )
    if pattern_classification is not None:
        if _pattern_match_is_false_positive(
            output_text, returncode, pattern_classification
        ):
            return None
        return _emit_classification(pattern_classification)

    if has_normal_termination(output_text or "") and returncode in (None, 0):
        return None

    if returncode not in (None, 0):
        return _emit_classification(
            OrcaFailureClassification(
                failure_type=OrcaFailureType.INTERRUPTED,
                evidence_lines=_fallback_evidence(scanned_output, scanned_stderr, returncode),
                returncode=returncode,
                retry_hint="no_retry",
                summary=f"Unrecognized ORCA failure (rc={returncode})",
            )
        )

    return None


def _match_output_patterns(
    output_text: str,
    *,
    returncode: int | None,
    max_cycles_opt: int | None,
    scf_converged: bool | None,
    evidence_context_lines: int,
) -> OrcaFailureClassification | None:
    if not output_text:
        return None

    ordered_patterns = (
        (
            OrcaFailureType.XTB_BINARY_UNAVAILABLE,
            (_XTB_BINARY_UNAVAILABLE_RE,),
            "provide_xtb_next_to_orca",
            "ORCA GFN-xTB requires the xtb binary inside the ORCA directory",
        ),
        (
            OrcaFailureType.SCF_SUBPROCESS_ABORTED,
            (_SCF_SUBPROCESS_QCMSG_RE, _SCF_SUBPROCESS_ABORTING_RE),
            "restart_geometry_fresh_scratch",
            "ORCA SCF subprocess aborted",
        ),
        (
            OrcaFailureType.SCF_NOT_CONVERGED,
            (_SCF_NOT_CONVERGED_RE, _SCF_CONVERGENCE_HALTED_RE, _SCF_NOT_CONVERGED_LOWER_RE),
            "slowconv_high_scfiter",
            "ORCA SCF did not converge" if scf_converged is False else "ORCA SCF did not converge",
        ),
        (
            OrcaFailureType.SCF_GRADIENT_ABORTED,
            (_SCF_GRADIENT_ABORT_RE, _SCF_GRADIENT_ABNORMAL_RE),
            "restart_geometry_fresh_scratch",
            "ORCA gradient evaluation aborted",
        ),
        (
            OrcaFailureType.MPI_RUNTIME_FAILED,
            (_MPI_RUNTIME_RE, _OPENMPI_RUNTIME_RE, _MPI_RANK_KILLED_RE, _MPI_PRIMARY_ABORT_RE, _MPI_PROCS_KILLED_RE),
            "reduce_ranks",
            "ORCA MPI runtime failed",
        ),
        (
            OrcaFailureType.MEMORY_PRESSURE,
            (_MEMORY_RE, _MP2_ABORT_RE, _RI_ABORT_RE, _OOM_PROCESS_RE, _MAXCORE_RE),
            "increase_memory",
            "ORCA failed due to memory pressure",
        ),
        (
            OrcaFailureType.GBW_VERSION_MISMATCH,
            (_GBW_VERSION_MISMATCH_RE,),
            "regenerate_gbw_for_orca_version",
            "ORCA GBW file is from a different ORCA version",
        ),
        (
            OrcaFailureType.INPUT_NON_ASCII,
            (_INPUT_NON_ASCII_RE,),
            "fix_input_encoding",
            "ORCA input contains non-ASCII characters",
        ),
        (
            OrcaFailureType.ORCA_ERROR_TERMINATION,
            (_ORCA_ERROR_TERMINATION_RE,),
            "no_retry",
            "ORCA finished by error termination",
        ),
        (
            OrcaFailureType.EXTERNAL_KILLED,
            (_EXTERNAL_KILLED_RE,),
            "check_scheduler",
            "ORCA job killed by scheduler/shell",
        ),
        (
            OrcaFailureType.THREAD_SPAWN_FAILED,
            (_THREAD_SPAWN_FAILED_RE,),
            "reduce_threads",
            "ORCA OpenMP thread creation failed",
        ),
        (
            OrcaFailureType.INTERRUPTED,
            (_INTERRUPTED_ABNORMAL_RE, _INTERRUPTED_ERROR_RE),
            "no_retry",
            "ORCA terminated abnormally",
        ),
        (
            OrcaFailureType.GEOMETRY_OPTIMIZATION_NOT_CONVERGED,
            (_OPT_NOT_CONVERGED_RE, _OPT_MAX_CYCLES_RE),
            "restart_from_last_geometry",
            _geometry_summary(max_cycles_opt),
        ),
    )

    for failure_type, patterns, retry_hint, summary in ordered_patterns:
        for compiled in patterns:
            if compiled.search(output_text):
                evidence_lines = extract_evidence(
                    output_text,
                    compiled.pattern,
                    context_lines=evidence_context_lines,
                )
                return OrcaFailureClassification(
                    failure_type=failure_type,
                    evidence_lines=evidence_lines or (_trim_line(compiled.pattern),),
                    returncode=returncode,
                    retry_hint=retry_hint,
                    summary=summary,
                )
    return None


def _signal_evidence(stderr_text: str, returncode: int | None) -> tuple[str, ...]:
    evidence: list[str] = []
    for pattern in (_SIGNAL_MEMORY_PATTERN, _SIGNAL_MPI_PATTERN, _SIGNAL_RUNTIME_PATTERN):
        evidence.extend(extract_evidence(stderr_text, pattern))
    if returncode is not None:
        evidence.append(_trim_line(f"process exited with signal {abs(returncode)}"))
    if not evidence and stderr_text.strip():
        evidence.append(_trim_line(stderr_text.strip().splitlines()[-1]))
    return tuple(dict.fromkeys(evidence)) or ("process terminated by signal",)


def _fallback_evidence(output_text: str, stderr_text: str, returncode: int | None) -> tuple[str, ...]:
    evidence: list[str] = []
    if stderr_text.strip():
        evidence.append(_trim_line(stderr_text.strip().splitlines()[-1]))
    if output_text.strip():
        evidence.append(_trim_line(output_text.strip().splitlines()[-1]))
    if returncode is not None:
        evidence.append(_trim_line(f"process returncode={returncode}"))
    return tuple(dict.fromkeys(evidence))


def _geometry_summary(max_cycles_opt: int | None) -> str:
    if max_cycles_opt is None:
        return "ORCA geometry optimization did not converge"
    return f"ORCA geometry optimization did not converge within {max_cycles_opt} cycles"


def _emit_classification(
    classification: OrcaFailureClassification,
) -> OrcaFailureClassification:
    logger.info(
        "Classified ORCA failure as %s (rc=%s): %s",
        classification.failure_type.value,
        classification.returncode,
        classification.summary,
    )
    return classification


@lru_cache(maxsize=1)
def _load_default_settings() -> OrcaFailureClassifierSettings:
    config_path = Path(__file__).resolve().parents[2] / "config" / "defaults.yaml"
    try:
        with config_path.open("r", encoding="utf-8") as handle:
            config = yaml.safe_load(handle) or {}
    except OSError as exc:
        logger.debug("Falling back to built-in ORCA classifier defaults: %s", exc)
        return OrcaFailureClassifierSettings()
    section = config.get("orca_failure_classification", {}) if isinstance(config, Mapping) else {}
    if not isinstance(section, Mapping):
        section = {}
    defaults = OrcaFailureClassifierSettings()
    return OrcaFailureClassifierSettings(
        enabled=bool(section.get("enabled", defaults.enabled)),
        scan_tail_chars=_coerce_positive_int(section.get("scan_tail_chars", defaults.scan_tail_chars), defaults.scan_tail_chars),
        evidence_context_lines=_coerce_nonnegative_int(
            section.get("evidence_context_lines", defaults.evidence_context_lines),
            defaults.evidence_context_lines,
        ),
    )


def _tail_text(text: str, max_chars: int) -> str:
    if max_chars <= 0 or len(text) <= max_chars:
        return text
    return text[-max_chars:]


def _trim_line(line: str) -> str:
    return line[:300]


def _coerce_positive_int(value: Any, default: int) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return default
    return parsed if parsed > 0 else default


def _coerce_nonnegative_int(value: Any, default: int) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return default
    return parsed if parsed >= 0 else default
