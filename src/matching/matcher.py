from typing import (
    Dict,
    Iterable,
    List,
    Optional,
    Tuple
)
from src.data_types import BGC_Variant, NRP_Variant, NRP_Fragment, LogProb
from src.config import MatchingConfig
from src.matching.detailed_hmm import DetailedHMM
from src.matching.hmm_auxiliary_types import HMM
from src.matching.match_type import Match, Match_BGC_Variant_Info, Match_NRP_Variant_Info
from src.matching.alignment_type import Alignment, alignment_score, combined_alignments_score
from src.matching.viterbi_algorithm import get_hmm_score, get_opt_path_with_emissions
from src.rban_parsing.get_linearizations import NRP_Linearizations, Linearization
from src.rban_parsing.rban_monomer import rBAN_Monomer
from src.generic.combinatorics import split_sequence_subseqs
from src.monomer_names_helper import MonomerNamesHelper
from collections import defaultdict
from itertools import chain, islice
from joblib import delayed, Parallel


def get_best_linearizations_for_nrp(hmm: HMM,
                                    nrp_linearizations: NRP_Linearizations,
                                    monomer_names_helper: MonomerNamesHelper,
                                    detailed_hmm: DetailedHMM,
                                    use_anchors_heuristic: bool = False) \
    -> Tuple[float, List[Linearization]]:
    if use_anchors_heuristic:
        _get_hmm_score = lambda mons: alignment_score(detailed_hmm.get_alignment(mons))
    else:
        _get_hmm_score = lambda mons: get_hmm_score(hmm,
                                                    [monomer_names_helper.mon_to_int[mon.to_base_mon()]
                                                     for mon in mons])
    def linearization_score(linearization: List[rBAN_Monomer]) -> LogProb:
        #mons = [monomer_names_helper.mon_to_int[mon.to_base_mon()]
        #        for mon in linearization]
        return _get_hmm_score(linearization)

    best_non_iterative_score, best_non_iterative_linearization = \
        max(((linearization_score(non_iterative_linearization), non_iterative_linearization)
            for non_iterative_linearization in nrp_linearizations.non_iterative),
            key=lambda x: x[0])

    best_iterative_score = float('-inf')
    best_iterative_linearizations = []
    for groups_linearizations in nrp_linearizations.iterative:
        split_score = 0
        split_linearizations = []
        for group_linearizations in groups_linearizations:
            group_score, group_linearization = max(((linearization_score(linearization), linearization)
                                                   for linearization in group_linearizations),
                                                    key=lambda x: x[0])
            split_linearizations.append(group_linearization)
            split_score += group_score
        if split_score > best_iterative_score:
            best_iterative_score = split_score
            best_iterative_linearizations = split_linearizations

    if best_non_iterative_score > best_iterative_score:
        return best_non_iterative_score, [best_non_iterative_linearization]
    else:
        return best_iterative_score, best_iterative_linearizations


def get_matches_for_hmm(detailed_hmm: DetailedHMM,
                        nrp_linearizations: Dict[Match_NRP_Variant_Info, NRP_Linearizations],
                        max_num_matches_per_bgc_variant: int,
                        log=None) -> List[Match]:
    max_num_matches_per_bgc_variant = max_num_matches_per_bgc_variant \
        if max_num_matches_per_bgc_variant != 0 else None
    if log is not None:
        log.info(f'Processing BGC {detailed_hmm.bgc_variant.genome_id} variant {detailed_hmm.bgc_variant.variant_idx}')

    hmm = detailed_hmm.to_hmm()
    matched_nrps_with_linearizations = []
    for nrp_info, nrp_linearizations in nrp_linearizations.items():
        score, linearizations = get_best_linearizations_for_nrp(hmm, nrp_linearizations,
                                                                detailed_hmm.hmm_helper.monomer_names_helper,
                                                                detailed_hmm,
                                                                use_anchors_heuristic=True)
        matched_nrps_with_linearizations.append((score, nrp_info, linearizations))

    best_matched_linearizations = list(islice(sorted(matched_nrps_with_linearizations,
                                                     key=lambda x: x[0],
                                                     reverse=True),
                                              max_num_matches_per_bgc_variant))

    def get_alignments_for_linearizations(linearizations: List[Linearization]) -> List[Alignment]:
        return [detailed_hmm.get_alignment(linearization) for linearization in linearizations]

    bgc_variant_info = Match_BGC_Variant_Info.from_bgc_variant(detailed_hmm.bgc_variant)
    best_matches = [Match(bgc_variant_info=bgc_variant_info,
                          nrp_variant_info=nrp_info,
                          alignments=get_alignments_for_linearizations(linearizations),
                          normalized_score=score)
                    for score, nrp_info, linearizations in best_matched_linearizations]
    return best_matches


def remove_bad_matches_for_nrps(matches_: Iterable[Match],
                                max_num_matches_per_nrp: int):
    matches = sorted(matches_,
                     key=lambda m: m.normalized_score, reverse=True)
    matches_per_nrp = defaultdict(list)
    for match in matches:
        if max_num_matches_per_nrp == 0 or \
            len(matches_per_nrp[match.nrp_variant_info]) < max_num_matches_per_nrp:
            matches_per_nrp[match.nrp_variant_info].append(match)

    return sorted(chain(*matches_per_nrp.values()),
                  key=lambda m: m.normalized_score, reverse=True)


def get_matches(hmms: List[DetailedHMM],
                nrp_linearizations: Dict[Match_NRP_Variant_Info, NRP_Linearizations],
                matching_cfg: MatchingConfig,
                num_threads: int = 1,
                log=None) -> List[Match]:
    total_linearizations = 0
    for nrp_linearization in nrp_linearizations.values():
        total_linearizations += len(nrp_linearization.non_iterative)
        total_linearizations += sum(len(group) for groups in nrp_linearization.iterative
                                    for group in groups)
    if log is not None:
        log.info(f'Matching {len(hmms)} BGC variants against {total_linearizations} NRP linearizations')
    matches = chain(*Parallel(n_jobs=num_threads)(delayed(get_matches_for_hmm)(hmm, nrp_linearizations,
                                                                               matching_cfg.max_num_matches_per_bgc,
                                                                               log)
                                                  for hmm in hmms))
    matches = remove_bad_matches_for_nrps(matches, matching_cfg.max_num_matches_per_nrp)
    # q: sort matches by normalized score
    max_num_matches = matching_cfg.max_num_matches if matching_cfg.max_num_matches != 0 else None
    return list(islice(sorted(matches, key=lambda m: m.normalized_score, reverse=True),
                       max_num_matches))
