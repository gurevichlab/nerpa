from pathlib import Path
from typing import (
    Dict,
    Iterable,
    List,
    Optional,
    Tuple, NamedTuple
)

from src import write_results
from src.data_types import (
    BGC_ID,
    BGC_Variant,
    BGC_Variant_ID,
    NRP_Variant,
    NRP_Fragment,
    LogProb, Prob
)
from src.config import MatchingConfig
from src.matching.detailed_hmm import DetailedHMM
from src.matching.hmm_auxiliary_types import HMM
from src.matching.hmm_match import HMM_Match
from src.matching.hmm_scoring_helper import HMMHelper
from src.matching.match_type import Match, NRP_Variant_ID
from src.matching.alignment_type import Alignment, alignment_score, combined_alignments_score
from src.matching.p_values_estimation import PValueEstimator
from src.matching.viterbi_algorithm import get_opt_path_with_score
from src.matching.hmm_checkpoints_heuristic import get_checkpoints
from src.pipeline.logger import NerpaLogger
from src.rban_parsing.get_linearizations import NRP_Linearizations, Linearization, num_linearizations
from src.rban_parsing.rban_monomer import rBAN_Monomer
from src.generic.combinatorics import split_sequence_subseqs
from src.monomer_names_helper import MonomerNamesHelper, MonomerResidue
from collections import defaultdict, Counter
from itertools import chain, islice
from joblib import delayed, Parallel
import math


def join_hmm_matches(hmm_matches: List[HMM_Match]) -> HMM_Match:
    score = sum(hmm_match.score for hmm_match in hmm_matches)
    bgc_variant_id = hmm_matches[0].bgc_variant_id
    nrp_id = hmm_matches[0].nrp_id

    nrp_linearizations = [linearization
                          for hmm_match in hmm_matches
                          for linearization in hmm_match.nrp_linearizations]
    optimal_paths = [opt_path
                     for hmm_match in hmm_matches
                     for opt_path in hmm_match.optimal_paths]
    return HMM_Match(score=score,
                     bgc_variant_id=bgc_variant_id,
                     nrp_id=nrp_id,
                     nrp_linearizations=nrp_linearizations,
                     optimal_paths=optimal_paths)


def get_best_match_for_nrp(hmm: HMM,
                           nrp_linearizations: NRP_Linearizations,
                           monomer_names_helper: MonomerNamesHelper,
                           detailed_hmm: DetailedHMM) -> HMM_Match:

    def match_for_linearization(linearization: List[rBAN_Monomer]) -> HMM_Match:
        mon_codes = [monomer_names_helper.mon_to_int[mon.to_base_mon()]
                     for mon in linearization]
        score, opt_path = get_opt_path_with_score(hmm, mon_codes)
        return HMM_Match(score=score,
                         bgc_variant_id=detailed_hmm.bgc_variant.bgc_variant_id,
                         nrp_id=nrp_linearizations.nrp_id,
                         nrp_linearizations=[[mon.rban_idx for mon in linearization]],
                         optimal_paths=[opt_path])

    best_noniterative_match = max((match_for_linearization(non_iterative_linearization)
                                   for non_iterative_linearization in nrp_linearizations.non_iterative),
                                  key=lambda match: match.score)

    best_iterative_score = float('-inf')
    best_iterative_match = []
    for groups_linearizations in nrp_linearizations.iterative:
        split_score = 0.0
        split_matches = []
        for group_linearizations in groups_linearizations:
            group_match= max((match_for_linearization(linearization)
                              for linearization in group_linearizations),
                             key=lambda match: match.score)
            split_matches.append(group_match)
            split_score += group_match.score
        if split_score > best_iterative_score:
            best_iterative_score = split_score
            best_iterative_match = join_hmm_matches(split_matches)

    if best_iterative_score > best_noniterative_match.score:
        return best_iterative_match
    else:
        return best_noniterative_match


def get_matches_for_hmm(detailed_hmm: DetailedHMM,
                        nrp_linearizations_all: List[NRP_Linearizations],
                        matching_cfg: MatchingConfig,
                        log=None) -> List[HMM_Match]:
    if log is not None:
        log.info(f'Processing BGC {detailed_hmm.bgc_variant.bgc_variant_id.bgc_id.genome_id} '
                 f'variant {detailed_hmm.bgc_variant.bgc_variant_id.variant_idx}')

    bgc_variant_id = detailed_hmm.bgc_variant.bgc_variant_id
    max_num_matches_per_bgc_variant = matching_cfg.max_num_matches_per_bgc \
        if matching_cfg.max_num_matches_per_bgc != 0 else None

    hmm = detailed_hmm.to_hmm()
    matches: List[HMM_Match] = [get_best_match_for_nrp(hmm, nrp_linearizations,
                                                       detailed_hmm.hmm_helper.monomer_names_helper,
                                                       detailed_hmm)
                                for nrp_linearizations in nrp_linearizations_all]


    return list(islice(sorted(matches,
                              key=lambda match: match.score,
                              reverse=True),
                       max_num_matches_per_bgc_variant))

NRP_ID = str

def filter_and_sort_matches(matches: List[HMM_Match],
                            hmms: List[DetailedHMM],
                            nrp_linearizations: List[NRP_Linearizations],
                            config: MatchingConfig) -> List[HMM_Match]:
    nrps_null_hypothesis_scores = {
        linearizations.nrp_id: null_hypothesis_score(linearizations, hmms[0].hmm_helper)
        for linearizations in nrp_linearizations
    }

    def normalized_match_score(match: HMM_Match) -> LogProb:
        """Calculate the normalized score of a match."""
        nrp_id = match.nrp_id
        null_hypothesis_score = nrps_null_hypothesis_scores[nrp_id]
        return match.score - null_hypothesis_score

    # Step 0: Keep only the top match per (BGC, NRP) pair
    matches_map: Dict[Tuple[BGC_ID, NRP_ID], List[HMM_Match]] = defaultdict(list)

    for match in matches:
        nrp_id = match.nrp_id
        bgc_id = match.bgc_variant_id.bgc_id
        matches_map[(bgc_id, nrp_id)].append(match)

    # Collect the best matches
    unique_matches = [max(matches_for_bgc_nrp, key=lambda m: normalized_match_score(m))
                      for matches_for_bgc_nrp in matches_map.values()]

    # Step 1: Sort all matches by descending score
    sorted_matches = sorted(unique_matches, key=lambda m: normalized_match_score(m), reverse=True)

    # Step 2: Filter by BGC
    bgc_counts: Dict[BGC_ID, int] = Counter()
    bgc_filtered = []

    for match in sorted_matches:
        bgc_id = match.bgc_variant_id.bgc_id
        if bgc_counts[bgc_id] < config.max_num_matches_per_bgc or config.max_num_matches_per_bgc == 0:
            bgc_filtered.append(match)
            bgc_counts[bgc_id] += 1

    # Step 3: Filter by NRP_ID
    nrp_counts: Dict[NRP_ID, int] = Counter()
    nrp_filtered = []
    for match in bgc_filtered:
        nrp_id = match.nrp_id
        if nrp_counts[nrp_id] < config.max_num_matches_per_nrp or config.max_num_matches_per_nrp == 0:
            nrp_filtered.append(match)
            nrp_counts[nrp_id] += 1

    # Step 4: Limit total number of matches
    if config.max_num_matches > 0:
        nrp_filtered = nrp_filtered[:config.max_num_matches]

    return nrp_filtered


def get_all_hmm_matches(hmms: List[DetailedHMM],
                    nrp_linearizations: List[NRP_Linearizations],
                    matching_cfg: MatchingConfig,
                    num_threads: int = 1,
                    log: Optional[NerpaLogger] = None) -> List[HMM_Match]:
    if log is not None:
        total_linearizations = sum(num_linearizations(nrp_linearization)
                                   for nrp_linearization in nrp_linearizations)
        log.info(f'Matching {len(hmms)} BGC variants against {total_linearizations} NRP linearizations')

    return list(
        chain(*Parallel(n_jobs=num_threads,  backend="loky")(
            delayed(get_matches_for_hmm)(hmm, nrp_linearizations,
                                         matching_cfg,
                                         log=None)
            for hmm in hmms))
    )

def null_hypothesis_score(nrp_linearizations: NRP_Linearizations,
                          hmm_scoring_helper: HMMHelper) -> LogProb:
        monomers = nrp_linearizations.non_iterative[0]  # any linearization will do, I just need the set of monomers
        return sum(hmm_scoring_helper.monomer_default_score(mon.to_base_mon())
                   for mon in monomers)

def get_hmm_matches(hmms: List[DetailedHMM],
                    nrp_linearizations: List[NRP_Linearizations],
                    matching_cfg: MatchingConfig,
                    num_threads: int = 1,
                    log: Optional[NerpaLogger] = None) -> List[HMM_Match]:
    debug = False
    matches = get_all_hmm_matches(hmms, nrp_linearizations, matching_cfg, num_threads, log)

    if log is not None:
        log.info('Matches obtained.')
        
    if debug:
        write_results.draw_hmms_with_optimal_paths(hmms, matches,
                                                   Path('/home/ilianolhin/git/nerpa2/hmms'))


    if log is not None:
        log.info('Estimating p-values for matches...')

    DetailedHMM.set_p_value_estimators_for_hmms(hmms, num_threads)

    if log is not None:
        log.info('P-values obtained. Filtering and sorting matches...')

    matches = filter_and_sort_matches(matches,
                                      hmms,
                                      nrp_linearizations,
                                      matching_cfg)

    if log is not None:
        log.info(f'Obtained {len(matches)} matches after filtering and sorting.')

    return matches
