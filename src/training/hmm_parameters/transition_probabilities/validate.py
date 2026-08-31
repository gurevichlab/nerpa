import math
from functools import partial
from logging import Logger
from typing import (
    Dict,
    NamedTuple,
    List,
    Tuple
)
from src.general_type_aliases import LogProb, Prob
    
from src.hmm.hmm_constructor.hmm_constructor_state_edge_context_relations import (
    STATE_TYPE_TO_EDGE_TYPES,
    RELEVANT_GENOMIC_CONTEXT,
    NORMAL_FLOW_EDGE_TYPES,
)
from src.hmm.hmm_constructor.hmm_constructor_state_edge_context_relations import all_possible_genomic_contexts
from src.training.hmm_parameters.short_aliases import *


def transition_weights_errors(
        transition_weights: Dict[ST, Dict[GC, Dict[ET, Prob]]],
) -> List[str]:
    error_msgs = []
    states_with_transitions = (
            set(STATE_TYPE_TO_EDGE_TYPES.keys())
            - {ST.FINAL}
    )

    if set(transition_weights.keys()) != states_with_transitions:
        error_msgs.append(
            'Estimated transition probabilities do not cover all state types. '
            f'Missing:\n'
            f'{states_with_transitions - set(transition_weights.keys())}'
        )


    for state_type, ctx_dict in transition_weights.items():
        missing_ctxts = (
            set(all_possible_genomic_contexts(state_type)) - set(ctx_dict.keys())
        )
        if missing_ctxts:
            error_msgs.append(
                f'Log probabilities for state type {state_type.name} '
                'do not cover all genomic contexts. Missing:\n'
                f'{missing_ctxts}'
            )

        for gc, edge_type_log_probs in ctx_dict.items():
            missing_edge_types = (
                STATE_TYPE_TO_EDGE_TYPES[state_type] - set(edge_type_log_probs.keys())
            )
            if missing_edge_types:
                error_msgs.append(
                    f'Log probabilities for state type {state_type.name} '
                    f'and genomic context {[f.name for f in gc]} '
                    f'do not cover all edge types. '
                    f'Missing:\n'
                 f'{missing_edge_types}'
                )

            total_prob = sum(pow(math.e, lp) for lp in edge_type_log_probs.values())
            if not math.isclose(total_prob, 1.0, rel_tol=1e-6):
                error_msgs.append(
                    f'Log probabilities for state type {state_type.name} and genomic context {[f.name for f in gc]} '
                    f'do not sum to 1 (sum={total_prob}).'
                )

    return error_msgs 
