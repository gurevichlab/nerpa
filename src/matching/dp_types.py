from typing import NamedTuple, Optional
import numpy.typing as npt
from src.matching.matching_types_alignment_step import AlignmentStepType
from src.data_types import LogProb


class DP_State(NamedTuple):
    module_pos: int
    monomer_pos: int
    gene_reps: int
    module_reps: int


class DP_Value(NamedTuple):
    score: LogProb
    parent: DP_State
    step_type: AlignmentStepType

DP_Table = npt.NDArray[Optional[DP_Value]]

def valid_state(dp_table: DP_Table, state: DP_State) -> bool:
    return all(0 <= idx < dp_table.shape[i]
               for i, idx in enumerate(state)) and \
        dp_table[state] is not None

START_STATE = DP_State(0,0,0,0)
