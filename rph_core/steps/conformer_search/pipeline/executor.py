from pathlib import Path
from typing import Callable, List, Tuple

from rph_core.steps.conformer_search.candidates import CandidateSet, ConformerCandidate, candidate_set_from_paths
from rph_core.steps.conformer_search.protocols import ProtocolSpec
from rph_core.steps.conformer_search.pipeline.stages.final_opt_sp import FinalOptSPResult, execute_final_opt_sp


class PipelineExecutor:
    def __init__(self, protocol_spec: ProtocolSpec):
        self.protocol_spec = protocol_spec

    def stage_sequence(self) -> List[str]:
        sequence = ["final_opt_sp"] if self.protocol_spec.final_opt_sp_enabled else []
        return sequence

    def execute_final_opt_sp(
        self,
        *,
        candidates: List[Path],
        runner: Callable[[List[ConformerCandidate]], Tuple[Path, float]],
    ) -> FinalOptSPResult:
        candidate_set = candidate_set_from_paths(
            candidates,
            source_protocol=self.protocol_spec.name,
        )
        result = execute_final_opt_sp(
            candidate_set=candidate_set,
            runner=runner,
            protocol_spec=self.protocol_spec,
        )
        return result

    def execute_handoff(
        self,
        *,
        candidate_set: CandidateSet,
        runner: Callable[[List[ConformerCandidate]], Tuple[Path, float]],
    ) -> FinalOptSPResult:
        result = execute_final_opt_sp(
            candidate_set=candidate_set,
            runner=runner,
            protocol_spec=self.protocol_spec,
        )
        return result
