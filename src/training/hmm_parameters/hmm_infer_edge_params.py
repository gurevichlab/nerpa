from typing import (
    Dict,
    List,
    Tuple,
    Optional,
    Union
)
from src.antismash_parsing.genomic_context import ModuleGenomicContextFeature
from src.data_types import LogProb
from src.matching.hmm_auxiliary_types import (
    DetailedHMMEdgeType, GenomicContext
)
from src.training.hmm_parameters.training_types import ChoicesCnts
from src.training.hmm_parameters.filter_edge_data import EDGE_TYPE_DEPENDENCIES
from math import log
from more_itertools import powerset


def infer_edge_params(edge_choices_cnts: Dict[DetailedHMMEdgeType,
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

    print('WARNING: not enough data for INSERT_AT_START,'
          ' so I use the same as for INSERT')
    edge_weight[ET.INSERT_AT_START] = from_cnts(edge_choices_cnts[ET.INSERT])

    edge_weight[ET.START_INSERTING] = from_cnts(edge_choices_cnts[ET.START_INSERTING])

    print('WARNING: not enough data for influence of PKS_DOWNSTREAM'
          ' on START_INSERTING, so I ignore it')
    for gc in edge_weight[ET.START_INSERTING]:
        if MGF.PKS_DOWNSTREAM in gc:
            upd_gc = tuple(feature for feature in gc if feature != MGF.PKS_DOWNSTREAM)
            edge_weight[ET.START_INSERTING][gc] = edge_weight[ET.START_INSERTING][upd_gc]

    edge_weight[ET.INSERT] = from_cnts(edge_choices_cnts[ET.INSERT])
    edge_weight[ET.START_INSERTING_AT_END] = from_cnts(edge_choices_cnts[ET.START_INSERTING_AT_END])

    print('WARNING: not enough data for influence of PKS_DOWNSTREAM'
          ' on START_INSERTING_AT_END, so I ignore it')
    edge_weight[ET.START_INSERTING_AT_END][(MGF.PKS_DOWNSTREAM,)] = edge_weight[ET.START_INSERTING_AT_END][()]

    edge_weight[ET.INSERT_AT_END] = from_cnts(edge_choices_cnts[ET.INSERT_AT_END])

    ##### SKIPS
    edge_weight[ET.SKIP_MODULE] = from_cnts(edge_choices_cnts[ET.SKIP_MODULE])
    edge_weight[ET.START_SKIP_MODULES_AT_START] = from_cnts(edge_choices_cnts[ET.START_SKIP_MODULES_AT_START])

    skip_at_start_runs = edge_choices_cnts[ET.START_SKIP_MODULES_AT_START][()].CHOSEN
    skip_at_start_total = edge_choices_cnts[ET.SKIP_MODULE_AT_START][()].CHOSEN  # ONLY_A_DOMAIN skips not counted
    avg = skip_at_start_total / skip_at_start_runs if skip_at_start_runs > 0 else 0
    skip_at_start_prob = 1 / (1 + avg)  # geometric distribution

    # all occurences of MGF.ONLY_A_DOMAIN *inside* BGC
    only_a_domains_middle = sum(edge_choices_cnts[ET.SKIP_MODULE][(MGF.ONLY_A_DOMAIN,)])  # CHOSEN + NOT_CHOSEN
    only_a_domains_start = edge_choices_cnts[ET.SKIP_MODULE_AT_START][(MGF.ONLY_A_DOMAIN,)].CHOSEN
    only_a_domains_end = edge_choices_cnts[ET.SKIP_MODULE_AT_END][(MGF.ONLY_A_DOMAIN,)].CHOSEN
    if (only_a_domains_middle + only_a_domains_start + only_a_domains_end) > 0:
        only_a_domain_skip_prob = (only_a_domains_start + only_a_domains_end) / (only_a_domains_middle + only_a_domains_start + only_a_domains_end)
    else:
        only_a_domain_skip_prob = 0.5


    edge_weight[ET.SKIP_MODULE_AT_START] = {(): log(skip_at_start_prob),
                                            (MGF.ONLY_A_DOMAIN,): log(only_a_domain_skip_prob)}

    edge_weight[ET.END_MATCHING] = from_cnts(edge_choices_cnts[ET.END_MATCHING])
    edge_weight[ET.SKIP_MODULE_END_MATCHING] = {
        gc: skip_log_prob + edge_weight[ET.END_MATCHING][()]
        for gc, skip_log_prob in edge_weight[ET.SKIP_MODULE].items()
    }

    skip_at_end_runs = (sum(cnts.CHOSEN for cnts in edge_choices_cnts[ET.SKIP_MODULE_END_MATCHING].values())
                        + edge_choices_cnts[ET.END_MATCHING][()].CHOSEN)
    skip_at_end_total = edge_choices_cnts[ET.SKIP_MODULE_AT_END][()].CHOSEN  # ONLY_A_DOMAIN skips not counted
    avg = skip_at_end_total / skip_at_end_runs if skip_at_end_runs > 0 else 0
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