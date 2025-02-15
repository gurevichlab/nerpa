from typing import (
    Dict,
    List,
    Tuple,
    Optional,
    Union
)
from src.antismash_parsing.location_features import ModuleLocFeature, GeneLocFeature, BGC_Fragment_Loc_Feature
from src.matching.matcher_viterbi_types import (
    DetailedHMMEdgeType
)
from src.training.statistical_inference import *

SingleFeatureContext = Union[ModuleLocFeature, GeneLocFeature, BGC_Fragment_Loc_Feature, None]


def infer_edge_params(data: Dict[DetailedHMMEdgeType, Dict[Tuple[SingleFeatureContext, ...], Tuple[int, int]]]):
    with open('edge_data.txt', 'w') as f:
        f.write(str(data))
    default_probability: float = 0.01  # near impossible event
    params: Dict[DetailedHMMEdgeType, Dict[SingleFeatureContext, float]] = {}

    # START_INSERTING_AT_START
    data_ins_st = data[DetailedHMMEdgeType.START_INSERTING_AT_START]
    estimated_probabilities: Dict[Union[SingleFeatureContext, str], float] = {}

    # Context 0: []
    estimated_probabilities[None] = default_probability
    base_prob = estimated_probabilities[None]  # to keep things concise

    # Context 1: ['START_OF_BGC']
    failures, successes = data_ins_st[(ModuleLocFeature.START_OF_BGC,)]
    estimated_probabilities[ModuleLocFeature.START_OF_BGC] = (
        (successes + 1) / (failures + successes + 2) - base_prob
    )

    # Context 2: ['START_OF_FRAGMENT']
    failures, successes = data_ins_st[(ModuleLocFeature.START_OF_BGC,)]
    total = 10  # prior pseudocounts / skepticism level
    alpha = total * (successes / (failures + successes))
    beta = total - alpha

    if (ModuleLocFeature.START_OF_FRAGMENT,) in data_ins_st:
        failures, successes = data_ins_st[(ModuleLocFeature.START_OF_FRAGMENT,)]
        alpha += successes
        beta += failures

        estimated_probabilities[ModuleLocFeature.START_OF_FRAGMENT] = (
            alpha / (alpha + beta)
        ) - base_prob
    else:
        estimated_probabilities[ModuleLocFeature.START_OF_FRAGMENT] = default_probability

    # Context 3: [ModuleLocFeature.PKS_UPSTREAM_PREV_GENE]
    prior = UNIFORM

    failures, successes = data_ins_st[(ModuleLocFeature.START_OF_BGC, ModuleLocFeature.PKS_UPSTREAM_PREV_GENE)]
    a = (1 - base_prob) * (1 - estimated_probabilities[ModuleLocFeature.START_OF_BGC])
    b = 1 - a
    f = (a, b)
    post = bayesian_inference(prior, f, failures, successes)

    failures, successes = data_ins_st[(ModuleLocFeature.START_OF_FRAGMENT, ModuleLocFeature.PKS_UPSTREAM_PREV_GENE)]
    a = (1 - base_prob) * (1 - estimated_probabilities[ModuleLocFeature.START_OF_FRAGMENT])
    b = 1 - a
    f = (a, b)
    post = bayesian_inference(post, f, failures, successes)

    estimated_probabilities[ModuleLocFeature.PKS_UPSTREAM_PREV_GENE] = expectation(post)

    # Context 4: [ModuleLocFeature.PKS_UPSTREAM_SAME_GENE]
    prior = UNIFORM

    failures, successes = data_ins_st[(ModuleLocFeature.START_OF_BGC, ModuleLocFeature.PKS_UPSTREAM_SAME_GENE)]
    a = (1 - base_prob) * (1 - estimated_probabilities[ModuleLocFeature.START_OF_BGC])
    b = 1 - a
    f = (a, b)
    post = bayesian_inference(prior, f, failures, successes)

    estimated_probabilities[ModuleLocFeature.PKS_UPSTREAM_SAME_GENE] = expectation(post)

    params[DetailedHMMEdgeType.START_INSERTING_AT_START] = estimated_probabilities

    # START_INSERTING
    data_ins = data[DetailedHMMEdgeType.START_INSERTING]
    estimated_probabilities = {}

    estimated_probabilities[None] = default_probability
    base_prob = estimated_probabilities[None]

    failures, successes = data_ins[(ModuleLocFeature.END_OF_BGC,)]
    estimated_probabilities[ModuleLocFeature.END_OF_BGC] = (
        (successes + 1) / (failures + successes + 2) - base_prob
    )

    failures, successes = data_ins[(ModuleLocFeature.END_OF_FRAGMENT,)]
    estimated_probabilities[ModuleLocFeature.END_OF_FRAGMENT] = (
        (successes + 1) / (failures + successes + 2) - base_prob
    )

    failures, successes = data_ins[(ModuleLocFeature.END_OF_GENE,)]
    estimated_probabilities[ModuleLocFeature.END_OF_GENE] = (
        (successes + 1) / (failures + successes + 2) - base_prob
    )

    prior = UNIFORM
    for loc_feature in [
        ModuleLocFeature.END_OF_BGC,
        ModuleLocFeature.END_OF_FRAGMENT,
        ModuleLocFeature.END_OF_GENE,
    ]:
        if (loc_feature, ModuleLocFeature.PKS_DOWNSTREAM_NEXT_GENE) not in data_ins:
            continue
        failures, successes = data_ins[(loc_feature, ModuleLocFeature.PKS_DOWNSTREAM_NEXT_GENE)]
        a = (1 - base_prob) * (1 - estimated_probabilities[loc_feature])
        b = 1 - a
        f = (a, b)
        post = bayesian_inference(prior, f, failures, successes)
        prior = post

    estimated_probabilities[ModuleLocFeature.PKS_DOWNSTREAM_NEXT_GENE] = expectation(post)

    prior = UNIFORM
    for loc_feature in [
        ModuleLocFeature.END_OF_BGC,
        ModuleLocFeature.END_OF_FRAGMENT,
        ModuleLocFeature.END_OF_GENE,
    ]:
        if (loc_feature, ModuleLocFeature.PKS_DOWNSTREAM_SAME_GENE) not in data_ins:
            continue
        failures, successes = data_ins[(loc_feature, ModuleLocFeature.PKS_DOWNSTREAM_SAME_GENE)]
        a = (1 - base_prob) * (1 - estimated_probabilities[loc_feature])
        b = 1 - a
        f = (a, b)
        post = bayesian_inference(prior, f, failures, successes)
        prior = post

    estimated_probabilities[ModuleLocFeature.PKS_DOWNSTREAM_SAME_GENE] = expectation(post)

    params[DetailedHMMEdgeType.START_INSERTING] = estimated_probabilities

    for edge_type, edge_data in data.items():
        if edge_type in [
            DetailedHMMEdgeType.START_INSERTING_AT_START,
            DetailedHMMEdgeType.START_INSERTING,
        ]:
            continue
        params[edge_type] = {}
        failures, successes = edge_data[()]
        prob_est = (successes + 1) / (failures + successes + 2)
        params[edge_type][None] = max(prob_est, default_probability)

    return params