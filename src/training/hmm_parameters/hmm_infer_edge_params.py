import math
from functools import partial
from logging import Logger
from typing import (
    Dict
)
from src.antismash_parsing.genomic_context import ModuleGenomicContextFeature
from src.general_type_aliases import LogProb
from src.hmm.hmm_auxiliary_types import (
    DetailedHMMEdgeType, GenomicContext, DetailedHMMStateType
)
from src.training.hmm_parameters.training_types import EdgeChoices_df, EdgeChoicesSchema
from src.hmm.hmm_constructor.hmm_constructor_state_edge_context_relations import (
    STATE_TYPE_TO_EDGE_TYPES,
    RELEVANT_GENOMIC_CONTEXT,
    NORMAL_FLOW_EDGE_TYPES,
)
from src.hmm.hmm_constructor.hmm_constructor_state_edge_context_relations import all_possible_genomic_contexts as all_possible_ctxts
from math import log

MCF = ModuleGenomicContextFeature

def get_log_probs_from_cnts(edge_type_cnts: Dict[DetailedHMMEdgeType, int],
                            min_default_probability: float) -> Dict[DetailedHMMEdgeType, LogProb]:
    log_probs = {}
    total_cnts = sum(edge_type_cnts.values())
    # print(f'Total counts: {total_cnts}')
    num_edge_types = len(edge_type_cnts)
    for edge_type, cnt in edge_type_cnts.items():
        pseudo_cnt_weight = (
                min_default_probability * total_cnts /
                (1 - num_edge_types * min_default_probability)
        )  # to ensure that each edge type has at least min_default_probability
        # print(f'Edge type: {edge_type.name}, count: {cnt}, pseudocount weight: {pseudo_cnt_weight}')
        prob = (cnt + pseudo_cnt_weight) / (total_cnts + num_edge_types * pseudo_cnt_weight)  # pseudocounts
        log_probs[edge_type] = log(prob)

    return log_probs


def check_log_probs(log_probs: Dict[DetailedHMMStateType, Dict[GenomicContext, Dict[DetailedHMMEdgeType, LogProb]]]):
    print('Sanity check of inferred log probabilities...')
    states_with_transitions = (
            set(STATE_TYPE_TO_EDGE_TYPES.keys())
            - {ST.FINAL}
    )
    assert set(log_probs.keys()) == states_with_transitions, \
        ('Log probabilities do not cover all state types.'
         f' Missing: '
         f'{states_with_transitions - set(log_probs.keys())}')

    for state_type, ctx_dict in log_probs.items():
        assert set(all_possible_ctxts(state_type)) == set(ctx_dict.keys()), \
            (f'Log probabilities for state type {state_type.name} do not cover all genomic contexts.'
             f' Missing: '
             f'{set(all_possible_ctxts(state_type)).difference(set(ctx_dict.keys()))}')

        for gc, edge_type_log_probs in ctx_dict.items():
            assert set(edge_type_log_probs.keys()) == STATE_TYPE_TO_EDGE_TYPES[state_type], \
                (f'Log probabilities for state type {state_type.name} and genomic context {[f.name for f in gc]} '
                 f'do not cover all edge types.'
                 f' Missing: '
                 f'{set(STATE_TYPE_TO_EDGE_TYPES[state_type]).difference(set(edge_type_log_probs.keys()))}')

            total_prob = sum(pow(math.e, lp) for lp in edge_type_log_probs.values())

            if state_type != ST.SKIPPING_MODULES_AT_END:
                assert abs(total_prob - 1.0) < 1e-6, \
                    (f'Log probabilities for state type {state_type.name} and genomic context {[f.name for f in gc]} '
                     f'do not sum to 1 (sum={total_prob}).')


ST = DetailedHMMStateType
ET = DetailedHMMEdgeType


def get_log_probs(edge_choices_per_gc: Dict[GenomicContext, Dict[DetailedHMMEdgeType, int]],
                  state_type: DetailedHMMStateType,
                  gc: GenomicContext,
                  min_num_cnts: int,
                  min_default_probability: float,
                  logger: Logger) -> Dict[DetailedHMMEdgeType, LogProb]:
    edge_type_cnts = edge_choices_per_gc.get(gc, {})
    if len(edge_type_cnts) == 1:
        return {
            next(iter(edge_type_cnts.keys())): 0.0,
        }

    try:
        normal_flow_edge = next(et for et in edge_type_cnts
                                if et in NORMAL_FLOW_EDGE_TYPES)
    except StopIteration:  # debug
        print(f'No normal flow edge found for {state_type.name} and {[f.name for f in gc]}\n'
              f'Edge counts: {edge_type_cnts}')
        raise
    other_edge = next(et for et in edge_type_cnts
                      if et != normal_flow_edge)

    if all(_edge_cnts.get(other_edge, 0) == 0
           for _gc, _edge_cnts in edge_choices_per_gc.items()):
        logger.info(f'{other_edge.name} never chosen for {state_type.name} -> '
                    f'asssigning min_default_probability.')
        return {
            normal_flow_edge: log(1 - min_default_probability),
            other_edge: log(min_default_probability)
        }

    if sum(edge_type_cnts.values()) >= min_num_cnts:
        return get_log_probs_from_cnts(
            edge_type_cnts,
            min_default_probability
        )

    logger.warning(f'Too few data for {state_type.name} and {[f.name for f in gc]} ')

    # Find a close genomic context with enough data
    try:
        gc_wo_feature = next(_gc
                             for _gc, _edge_cnts in edge_choices_per_gc.items()
                             if len(_gc) == len(gc) - 1 and sum(_edge_cnts.values()) >= min_num_cnts)
    except StopIteration:
        logger.error(f'No close genomic context with enough data for {state_type.name} and {[f.name for f in gc]} ')
        raise

    edge_type_cnts_wo_feature = edge_choices_per_gc.get(gc_wo_feature, {})

    log_probs_wo_feature = get_log_probs_from_cnts(
        edge_type_cnts_wo_feature,
        min_default_probability
    )

    log_probs_fict_cnts = get_log_probs_from_cnts(
        {
            normal_flow_edge: min_num_cnts - edge_type_cnts[other_edge],
            other_edge: edge_type_cnts[other_edge]
        },
        min_default_probability
    )

    if edge_type_cnts[other_edge] == 0:
        logger.info(f'{other_edge.name} never chosen -> '
                    f'asssigning probability from {tuple(f.name for f in gc_wo_feature)}.')
        return log_probs_wo_feature

    if log_probs_fict_cnts[other_edge] > log_probs_wo_feature[other_edge]:
        logger.info(f'Enough {other_edge.name} events -> '
                    f'adding fictional {normal_flow_edge.name}.')
        return log_probs_fict_cnts
    else:
        logger.info(f'Too few {other_edge.name} events -> '
                    f'using probabilities from {tuple(f.name for f in gc_wo_feature)}')
        return log_probs_wo_feature


ECS = EdgeChoicesSchema

def get_num_skip_at_start_runs(edge_choices_df: EdgeChoices_df) -> int:
    # return the number unique NRP_IDs with at least one SKIP_UNTIL_NEXT_TENTATIVE_ASSEMBLY_LINE_START edge type
    mask = (edge_choices_df[ECS.EDGE_TYPE] == ET.SKIP_UNTIL_NEXT_TENTATIVE_ASSEMBLY_LINE_START)
    nrp_ids_with_skip = edge_choices_df.loc[mask, ECS.NRP_ID].unique()
    return len(nrp_ids_with_skip)

def get_num_skip_at_end_runs(edge_choices_df: EdgeChoices_df) -> int:
    # return the number unique NRP_IDs with at least one SKIP_UNTIL_NEXT_TENTATIVE_ASSEMBLY_LINE_START edge type
    mask = (edge_choices_df[ECS.EDGE_TYPE] == ET.SKIP_MODULES_AT_END)
    nrp_ids_with_skip = edge_choices_df.loc[mask, ECS.NRP_ID].unique()
    return len(nrp_ids_with_skip)


def get_log_probs_all_gcs(edge_cnts_per_gc: Dict[GenomicContext, Dict[DetailedHMMEdgeType, int]],
                          state_type: DetailedHMMStateType,
                          min_num_cnts: int,
                          min_default_probability: float,
                          logger: Logger) -> Dict[GenomicContext, Dict[DetailedHMMEdgeType, LogProb]]:
    return {
        gc: get_log_probs(edge_cnts_per_gc,
                          state_type,
                          gc,
                          min_num_cnts,
                          min_default_probability,
                          logger)
        for gc in edge_cnts_per_gc.keys()
    }


def assign_skips_at_start_probs(edge_choices_per_state_cntxt: Dict[DetailedHMMStateType, Dict[GenomicContext, Dict[DetailedHMMEdgeType, int]]],
                                edge_choices_df: EdgeChoices_df,
                                log_probs: Dict[DetailedHMMStateType, Dict[GenomicContext, Dict[DetailedHMMEdgeType, LogProb]]],
                                min_num_cnts: int,
                                min_default_probability: float,
                                logger: Logger = None):
    # Idea: I first pretend that ASSEMBLY_LINE_START_MARKER is not important,
    # and pull up counts from contexts with and without it.
    # Then, I use the pulled-up counts to set the probability for CHOOSE_IF_START_MATCHING
    # when ASSEMBLY_LINE_START_MARKER is present,
    # and set the probability of START_MATCHING to -inf
    # when ASSEMBLY_LINE_START_MARKER is not present.

    _get_log_probs_all_gcs = partial(get_log_probs_all_gcs,  # for brevity
                                     state_type=ST.CHOOSE_IF_START_MATCHING,
                                     min_num_cnts=min_num_cnts,
                                     min_default_probability=min_default_probability,
                                     logger=logger)

    log_probs[ST.CHOOSE_IF_START_MATCHING] = {}

    # 1. Probability to start skipping at start
    edge_cnts_bgc_start = {
        gc_bgc_start: edge_choices_per_state_cntxt[ST.CHOOSE_IF_START_MATCHING].get(gc_bgc_start, {})
        for gc_bgc_start in all_possible_ctxts(ST.CHOOSE_IF_START_MATCHING)
        if MCF.START_OF_BGC in gc_bgc_start
    }
    log_probs[ST.CHOOSE_IF_START_MATCHING].update(
        _get_log_probs_all_gcs(edge_cnts_per_gc=edge_cnts_bgc_start)
    )

    # 2. Probability to continue skipping at start
    pulled_cnts_start_marker = {}
    for gc_with_start_marker in filter(lambda gc: MCF.START_OF_BGC not in gc and MCF.ASSEMBLY_LINE_START_MARKER in gc,
                                       all_possible_ctxts(ST.CHOOSE_IF_START_MATCHING)):
        gc_wo_start_marker = tuple(f for f in gc_with_start_marker if f != MCF.ASSEMBLY_LINE_START_MARKER)
        edge_cnts_with_start_marker = edge_choices_per_state_cntxt[ST.CHOOSE_IF_START_MATCHING].get(gc_with_start_marker, {})
        edge_cnts_wo_start_marker = edge_choices_per_state_cntxt[ST.CHOOSE_IF_START_MATCHING].get(gc_wo_start_marker, {})
        pulled_cnts_start_marker[gc_with_start_marker] = {
            edge_type: edge_cnts_with_start_marker.get(edge_type, 0) +
                       edge_cnts_wo_start_marker.get(edge_type, 0)
            for edge_type in STATE_TYPE_TO_EDGE_TYPES[ST.CHOOSE_IF_START_MATCHING]
        }

    log_probs[ST.CHOOSE_IF_START_MATCHING].update(
        _get_log_probs_all_gcs(edge_cnts_per_gc=pulled_cnts_start_marker)
    )


    # 3. Disable START_MATCHING when ASSEMBLY_LINE_START_MARKER is not present
    log_probs[ST.CHOOSE_IF_START_MATCHING].update({
        gc_wo_start_marker:
            {
                ET.START_MATCHING: float('-inf'),
                ET.SKIP_UNTIL_NEXT_TENTATIVE_ASSEMBLY_LINE_START: 0
            }
        for gc_wo_start_marker in all_possible_ctxts(ST.CHOOSE_IF_START_MATCHING)
        if MCF.START_OF_BGC not in gc_wo_start_marker
           and MCF.ASSEMBLY_LINE_START_MARKER not in gc_wo_start_marker
    }
    )


def assign_skips_at_end_probs(edge_choices_per_state_cntxt: Dict[DetailedHMMStateType, Dict[GenomicContext, Dict[DetailedHMMEdgeType, int]]],
                              edge_choices_wo_filtering: Dict[DetailedHMMStateType, Dict[GenomicContext, Dict[DetailedHMMEdgeType, int]]],
                              edge_choices_df: EdgeChoices_df,
                              log_probs: Dict[DetailedHMMStateType, Dict[GenomicContext, Dict[DetailedHMMEdgeType, LogProb]]],
                              min_num_cnts: int,
                              min_default_probability: float,
                              logger: Logger = None):
    # Idea: I pretend that HMM goes in reverse.
    # Starting from the FINAL state, it decides whether to skip modules at end or "finish" matching.
    # Then skips_continued/skips_started is the probability of skipping an individual module at the end.
    # When assigninng weights, I will sum up the log probabilities for each skipped module and assign that
    # to the ET.FINISH_MATCHING

    _get_log_probs_all_gcs = partial(get_log_probs_all_gcs,  # for brevity
                                     min_num_cnts=min_num_cnts,
                                     min_default_probability=min_default_probability,
                                     logger=logger)

    log_probs[ST.CHOOSE_IF_FINISH_MATCHING] = {}
    log_probs[ST.SKIPPING_MODULES_AT_END] = {}

    # 1. Probability to start skipping at end
    edge_cnts_finish_marker = {
        gc_finish_marker: edge_choices_per_state_cntxt[ST.CHOOSE_IF_FINISH_MATCHING].get(gc_finish_marker, {})
        for gc_finish_marker in all_possible_ctxts(ST.CHOOSE_IF_FINISH_MATCHING)
        if MCF.ASSEMBLY_LINE_FINISH_MARKER in gc_finish_marker
    }

    log_probs[ST.CHOOSE_IF_FINISH_MATCHING].update(
        _get_log_probs_all_gcs(edge_cnts_per_gc=edge_cnts_finish_marker,
                               state_type=ST.CHOOSE_IF_FINISH_MATCHING)
    )

    # 2. Disable FINISH_MATCHING when ASSEMBLY_LINE_FINISH_MARKER is not present
    log_probs[ST.CHOOSE_IF_FINISH_MATCHING].update({
        gc: {
                ET.FINISH_MATCHING: float('-inf'),
                ET.CONTINUE_MATCHING: 0
        }
        for gc in all_possible_ctxts(ST.CHOOSE_IF_FINISH_MATCHING)
        if MCF.ASSEMBLY_LINE_FINISH_MARKER not in gc
    })

    # 3. Probability to continue skipping at end
    # if we view a skip run at the end from the end to the start,
    # it comprises a number of SKIP_MODULES_AT_END edges followed by a FINISH_MATCHING edge
    # I view each SKIP_MODULES_AT_END edge as a decision to continue skipping,
    # and the FINISH_MATCHING edge as a decision to stop skipping
    edge_cnts_skipping = {}
    for gc in all_possible_ctxts(ST.SKIPPING_MODULES_AT_END):
        skips_finished = 0
        # ST.CHOOSE_IF_FINISH_MATCHING has different relevant contexts than ST.SKIPPING_MODULES_AT_END,
        # so I need to sum up counts from unfiltered data
        for gc_ in edge_choices_wo_filtering[ST.CHOOSE_IF_FINISH_MATCHING].keys():
            gc_rel = tuple(feature
                           for feature in gc_
                           if feature in RELEVANT_GENOMIC_CONTEXT[ST.SKIPPING_MODULES_AT_END])
            if gc_rel == gc:
                skips_finished += (
                    edge_choices_wo_filtering[ST.CHOOSE_IF_FINISH_MATCHING][gc_][ET.FINISH_MATCHING]
                )

        skips_continued = edge_choices_per_state_cntxt[ST.SKIPPING_MODULES_AT_END][gc][ET.SKIP_MODULES_AT_END]

        edge_cnts_skipping[gc] = {
                ET.START_MATCHING: skips_finished,  # ET.START_MATCHING is a placeholder. This edge doesn't exist in SKIPPING_MODULES_AT_END
                ET.SKIP_MODULES_AT_END: skips_continued,
            }

    # Note that ST.SKIPPING_MODULES_AT_END has a single edge,
    # and its weight is the log probability to continue skipping (not 0)
    # during HMM construction, I will move these probabilities to ST.CHOOSE_IF_FINISH_MATCHING
    # and set the weight of SKIPPING_MODULES_AT_END to 0
    skipping_log_probs = _get_log_probs_all_gcs(edge_cnts_per_gc=edge_cnts_skipping,
                                                state_type=ST.SKIPPING_MODULES_AT_END)
    log_probs[ST.SKIPPING_MODULES_AT_END] = {
        gc: {ET.SKIP_MODULES_AT_END: gc_log_probs[ET.SKIP_MODULES_AT_END]}  # remove placeholder edge
        for gc, gc_log_probs in skipping_log_probs.items()
    }



def infer_edge_params(edge_choices_per_state_cntxt: Dict[DetailedHMMStateType, Dict[GenomicContext, Dict[DetailedHMMEdgeType, int]]],
                      edge_choices_wo_filtering: Dict[DetailedHMMStateType, Dict[GenomicContext, Dict[DetailedHMMEdgeType, int]]],
                      edge_choices_df: EdgeChoices_df,
                      min_default_probability: float = 0.01,
                      min_num_cnts: int = 5,
                      logger: Logger = None) \
        -> Dict[DetailedHMMStateType, Dict[GenomicContext, Dict[DetailedHMMEdgeType, LogProb]]]:

    # enforce that all edge types are present

    for state_type in STATE_TYPE_TO_EDGE_TYPES.keys():
        for gc in all_possible_ctxts(state_type):
            if gc not in edge_choices_per_state_cntxt[state_type]:
                edge_choices_per_state_cntxt[state_type][gc] = {}
            for edge_type in STATE_TYPE_TO_EDGE_TYPES[state_type]:
                if edge_type not in edge_choices_per_state_cntxt[state_type][gc]:
                    edge_choices_per_state_cntxt[state_type][gc][edge_type] = 0

    log_probs = {}

    assign_skips_at_start_probs(edge_choices_per_state_cntxt,
                                edge_choices_df,
                                log_probs,
                                min_num_cnts,
                                min_default_probability,
                                logger)

    assign_skips_at_end_probs(edge_choices_per_state_cntxt,
                              edge_choices_wo_filtering,
                              edge_choices_df,
                              log_probs,
                              min_num_cnts,
                              min_default_probability,
                              logger)

    for state_type in STATE_TYPE_TO_EDGE_TYPES.keys():
        if (state_type in log_probs  # already assigned. Don't overwrite
                or state_type == ST.FINAL):  # no outgoing edges
            continue

        log_probs[state_type] = {}
        for gc in all_possible_ctxts(state_type):
            log_probs[state_type][gc] = get_log_probs(edge_choices_per_state_cntxt[state_type],
                                                      state_type,
                                                      gc,
                                                      min_num_cnts,
                                                      min_default_probability,
                                                      logger)

    check_log_probs(log_probs)

    return log_probs



'''
def infer_edge_params_old(edge_choices_cnts: Dict[DetailedHMMEdgeType,
                                              Dict[Optional[GenomicContext], ChoicesCnts]],
                      min_allowed_log_prob: LogProb=-5) \
        -> Dict[DetailedHMMEdgeType, Dict[GenomicContext, LogProb]]:
    print('Inferring edge parameters from counts...')
    ET = DetailedHMMEdgeType
    MGF = ModuleGenomicContextFeature
    edge_weight: Dict[DetailedHMMEdgeType, Dict[GenomicContext, LogProb]] = {}

    def from_cnts(cnts_dict: Dict[GenomicContext, ChoicesCnts]) -> Dict[GenomicContext, LogProb]:
        return {gc: log(cnts.to_prob(pseudocounts=True))
                for gc, cnts in cnts_dict.items()}

    ##### INSERTS
    edge_weight[ET.START_INSERTING_AT_START] = from_cnts(edge_choices_cnts[ET.START_INSERTING_AT_START])

    print('WARNING(ok): not enough data for INSERT_AT_START,'
          ' so I use the same as for INSERT')
    edge_weight[ET.INSERT_AT_START] = from_cnts(edge_choices_cnts[ET.INSERT])

    edge_weight[ET.START_INSERTING] = from_cnts(edge_choices_cnts[ET.START_INSERTING])

    print('WARNING(ok): not enough data for influence of PKS_DOWNSTREAM'
          ' on START_INSERTING, so I ignore it')
    for gc in edge_weight[ET.START_INSERTING]:
        if MGF.PKS_DOWNSTREAM in gc:
            upd_gc = tuple(feature for feature in gc if feature != MGF.PKS_DOWNSTREAM)
            edge_weight[ET.START_INSERTING][gc] = edge_weight[ET.START_INSERTING][upd_gc]

    edge_weight[ET.INSERT] = from_cnts(edge_choices_cnts[ET.INSERT])
    edge_weight[ET.START_INSERTING_AT_END] = from_cnts(edge_choices_cnts[ET.START_INSERTING_AT_END])

    print('WARNING(ok): not enough data for influence of PKS_DOWNSTREAM'
          ' on START_INSERTING_AT_END, so I ignore it')
    edge_weight[ET.START_INSERTING_AT_END][(MGF.PKS_DOWNSTREAM,)] = edge_weight[ET.START_INSERTING_AT_END][()]

    edge_weight[ET.INSERT_AT_END] = from_cnts(edge_choices_cnts[ET.INSERT_AT_END])

    ##### SKIPS
    edge_weight[ET.SKIP_MODULE] = from_cnts(edge_choices_cnts[ET.SKIP_MODULE])
    edge_weight[ET.START_SKIPPING_MODULES_AT_START] = from_cnts(edge_choices_cnts[ET.START_SKIPPING_MODULES_AT_START])

    skip_at_start_runs = edge_choices_cnts[ET.START_SKIPPING_MODULES_AT_START][()].CHOSEN
    skip_at_start_total = edge_choices_cnts[ET.SKIP_MODULE_AT_START][()].CHOSEN  # ONLY_A_DOMAIN skips not counted
    if skip_at_start_runs > 0:
        avg = skip_at_start_total / skip_at_start_runs
    else:
        print('WARNING(should happen only in debug): no data for SKIP_MODULE_AT_START')
        avg = 0
    skip_at_start_prob = 1 / (1 + avg)  # geometric distribution

    # all occurences of MGF.ONLY_A_DOMAIN *inside* BGC
    only_a_domains_middle = sum(edge_choices_cnts[ET.SKIP_MODULE][(MGF.ONLY_A_DOMAIN,)])  # CHOSEN + NOT_CHOSEN
    only_a_domains_start = edge_choices_cnts[ET.SKIP_MODULE_AT_START][(MGF.ONLY_A_DOMAIN,)].CHOSEN
    only_a_domains_end = edge_choices_cnts[ET.SKIP_MODULE_AT_END][(MGF.ONLY_A_DOMAIN,)].CHOSEN
    if (only_a_domains_middle + only_a_domains_start + only_a_domains_end) > 0:
        only_a_domain_skip_prob = (only_a_domains_start + only_a_domains_end) / (only_a_domains_middle + only_a_domains_start + only_a_domains_end)
    else:
        print('WARNING(should happen only in debug): no data for SKIP_MODULE_AT_START, ONLY_A_DOMAIN')
        only_a_domain_skip_prob = 0.5


    edge_weight[ET.SKIP_MODULE_AT_START] = {(): log(skip_at_start_prob),
                                            (MGF.ONLY_A_DOMAIN,): log(only_a_domain_skip_prob)}

    edge_weight[ET.END_MATCHING] = from_cnts(edge_choices_cnts[ET.END_MATCHING])
    if () not in edge_weight[ET.END_MATCHING]:
        print('WARNING(should only happen in debug): no data for END_MATCHING, setting to min allowed value')
        edge_weight[ET.END_MATCHING][()] = min_allowed_log_prob

    edge_weight[ET.SKIP_MODULE_END_MATCHING] = {
        gc: skip_log_prob + edge_weight[ET.END_MATCHING][()]
        for gc, skip_log_prob in edge_weight[ET.SKIP_MODULE].items()
    }

    skip_at_end_runs = (sum(cnts.CHOSEN for cnts in edge_choices_cnts[ET.SKIP_MODULE_END_MATCHING].values())
                        + edge_choices_cnts[ET.END_MATCHING][()].CHOSEN)
    skip_at_end_total = edge_choices_cnts[ET.SKIP_MODULE_AT_END][()].CHOSEN  # ONLY_A_DOMAIN skips not counted
    if skip_at_end_runs > 0:
        avg = skip_at_end_total / skip_at_end_runs
    else:
        print('WARNING(should happen only in debug): no data for SKIP_MODULE_AT_END')
        avg = 0
    skip_at_end_prob = 1 / (1 + avg)  # geometric distribution
    edge_weight[ET.SKIP_MODULE_AT_END] = {(): log(skip_at_end_prob),
                                            (MGF.ONLY_A_DOMAIN,): log(only_a_domain_skip_prob)}

    edge_weight[ET.SKIP_GENE] = from_cnts(edge_choices_cnts[ET.SKIP_GENE])

    ##### ITERATIONS
    edge_weight[ET.ITERATE_MODULE] = from_cnts(edge_choices_cnts[ET.ITERATE_MODULE])
    edge_weight[ET.ITERATE_GENE] = from_cnts(edge_choices_cnts[ET.ITERATE_GENE])

    ##### Check if all edge types and genomic contexts are present (for debugging)
    #### and ensure that the weights are above the min_allowed_log_prob
    for edge_type, relevant_mgfs in EDGE_TYPE_DEPENDENCIES.items():
        if edge_type not in edge_weight:
            raise Exception(f'Weights for edge type {edge_type} not computed')
        for mgf_set in powerset(relevant_mgfs):
            gc = tuple(sorted(mgf_set, key=lambda mgf: mgf.name))
            if gc not in edge_weight[edge_type]:
                print(f'WARNING: Weights for edge type {edge_type} and genomic context {gc} not computed'
                      f' - setting to min allowed value: {min_allowed_log_prob}')
            edge_weight[edge_type][gc] = max(edge_weight[edge_type].get(gc, min_allowed_log_prob),
                                             min_allowed_log_prob)


    return edge_weight
    
    
def handle_states_with_too_little_data(edge_choices_per_state_cntxt: Dict[DetailedHMMStateType, Dict[GenomicContext, Dict[DetailedHMMEdgeType, int]]],
                                       log_probs: Dict[DetailedHMMStateType, Dict[GenomicContext, Dict[DetailedHMMEdgeType, LogProb]]],
                                       min_num_cnts: int,
                                       min_default_probability: float,
                                       logger: Logger = None):
    # These are ad-hoc rules to handle states with too little data.
    for state_type, ctxt_dict in edge_choices_per_state_cntxt.items():
        for gc in all_possible_ctxts(state_type):
            edge_type_cnts = ctxt_dict.get(gc, {})
            if log_probs[state_type].get(gc, None) is None:
                total_cnts = sum(edge_type_cnts.values())
                logger.warning(f'WARNING: Not enough data for state type {state_type.name}'
                            f' and genomic context {[f.name for f in gc]} '
                            f'({total_cnts} data points).')


                if state_type == ST.INSERT:
                    # It seems it never happens, so I just assign min_default_probability
                    if any(_edge_type_cnts.get(ET.INSERT, 0) > 0  # in case it ever happens
                           for gc2, _edge_type_cnts in ctxt_dict.items()):
                        raise NotImplementedError()
                    log_probs[ST.INSERT][gc] = {
                        ET.INSERT: log(min_default_probability),
                        ET.FINISH_INSERTING: log(1 - min_default_probability)
                    }
                    continue

                if len(gc) <= 2 and MCF.PKS_DOWNSTREAM in gc:
                    # logger.info('Pulling up probabilities from context without PKS_DOWNSTREAM')
                    normal_flow_edge = next(et for et in STATE_TYPE_TO_EDGE_TYPES[state_type]
                                            if et in NORMAL_FLOW_EDGE_TYPES)
                    other_edge = next(et for et in STATE_TYPE_TO_EDGE_TYPES[state_type]
                                      if et != normal_flow_edge)
                    # Add enough fictional counts to normal_flow_edge to reach min_num_cnts
                    new_cnts = {normal_flow_edge: min_num_cnts - edge_type_cnts[other_edge],
                                other_edge: edge_type_cnts[other_edge]}

                    new_log_probs = get_log_probs_from_cnts(new_cnts,
                                                            min_default_probability)

                    # Choose whether to use new_cnts or cnts for context without PKS_DOWNSTREAM
                    # Pick the one with higher probability for other_edge
                    ctxt_wo_pks = tuple(feature for feature in gc if feature != MCF.PKS_DOWNSTREAM)
                    if sum(ctxt_dict.get(ctxt_wo_pks, {}).values()) < min_num_cnts:
                        raise NotImplementedError()

                    if new_log_probs[other_edge] > log_probs[state_type][ctxt_wo_pks][other_edge]:
                        logger.info(f'Adding fictional {normal_flow_edge.name} so that'
                                    f'the total counts reach {min_num_cnts}')
                        log_probs[state_type][gc] = new_log_probs
                    else:
                        logger.info(f'Using probabilities from context {ctxt_wo_pks}')
                        log_probs[state_type][gc] = log_probs[state_type][ctxt_wo_pks]
                    continue

                logger.warning(f'Not enough data to estimate probabilities '
                               f'for state type {state_type.name} '
                               f' and genomic context {[f.name for f in gc]}.')
                normal_flow_edge = next(et for et in STATE_TYPE_TO_EDGE_TYPES[state_type]
                                        if et in NORMAL_FLOW_EDGE_TYPES)
                other_edge = next(et for et in STATE_TYPE_TO_EDGE_TYPES[state_type]
                                  if et != normal_flow_edge)
                # Add enough fictional counts to normal_flow_edge to reach min_num_cnts
                new_cnts = {normal_flow_edge: min_num_cnts - edge_type_cnts[other_edge],
                            other_edge: edge_type_cnts[other_edge]}
                logger.info(f'Adding fictional {normal_flow_edge.name} so that'
                            f'the total counts reach {min_num_cnts}')
                log_probs[state_type][gc] = get_log_probs_from_cnts(new_cnts,
                                                                    min_default_probability)
'''