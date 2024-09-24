from typing import Iterable, List, Dict, NamedTuple
from itertools import chain, pairwise, product
from collections import defaultdict
from src.data_types import MonomerResidue, LogProb, BGC_Variant, NRP_Variant, NRP_Fragment
from src.generic.combinatorics import intersection_with_repeats
from src.matching.scoring_config import HeuristicMatchingConfig


def top_n_predictions(residue_scores: Dict[MonomerResidue, LogProb],
                      num_top_predictions: int) -> List[MonomerResidue]:
    return sorted(residue_scores, key=residue_scores.get, reverse=True)[:num_top_predictions]


def bgc_residues(bgc_variant: BGC_Variant,
                 num_top_predictions: int = 3) -> List[MonomerResidue]:
    return list(chain(*[top_n_predictions(bgc_module.residue_score, num_top_predictions)
                        for bgc_fragment in bgc_variant.fragments
                        for bgc_module in bgc_fragment]))


def nrp_residues(nrp_variant: NRP_Variant) -> List[MonomerResidue]:
    return [mon.residue for fragment in nrp_variant.fragments
            for mon in fragment.monomers]


def get_common_cons_pairs_cnt_nonit(bgc_variant: BGC_Variant, nrp_variant: NRP_Variant) -> int:
    nrp_pairs = chain(*(pairwise([mon.residue for mon in fragment.monomers])
                        for fragment in nrp_variant.fragments))

    bgc_set_pairs = pairwise(top_n_predictions(bgc_module.residue_score, 3)
                             for bgc_fragment in bgc_variant.fragments
                             for bgc_module in bgc_fragment)
    bgc_pairs = chain(*(product(s1, s2) for s1, s2 in bgc_set_pairs))
    return len(intersection_with_repeats(nrp_pairs, bgc_pairs))


def get_common_cons_pairs_cnt_nit(bgc_variant: BGC_Variant, nrp_variant: NRP_Variant) -> int:
    bgc_set_pairs = pairwise(top_n_predictions(bgc_module.residue_score, 3)
                             for bgc_fragment in bgc_variant.fragments
                             for bgc_module in bgc_fragment)
    bgc_pairs = chain(*(product(s1, s2) for s1, s2 in bgc_set_pairs))
    return sum(len(intersection_with_repeats(pairwise([mon.residue for mon in fragment.monomers]),
                                             bgc_pairs))
               for fragment in nrp_variant.fragments)


def heuristic_match_discard(bgc_variant: BGC_Variant,
                            nrp_fragments: List[NRP_Fragment],
                            heuristic_matching_cfg: HeuristicMatchingConfig) -> bool:
    nrp_len = sum(len(nrp_fragment.monomers) for nrp_fragment in nrp_fragments)
    bgc_len = sum(len(bgc_fragment) for bgc_fragment in bgc_variant.fragments)

    bgc_res = bgc_residues(bgc_variant, heuristic_matching_cfg.NUM_TOP_PREDICTIONS)
    nrp_mons = [mon.residue for fragment in nrp_fragments for mon in fragment.monomers]

    bgc_pairs = chain(*(product(top_n_predictions(bgc_module1.residue_score, heuristic_matching_cfg.NUM_TOP_PREDICTIONS),
                                top_n_predictions(bgc_module2.residue_score, heuristic_matching_cfg.NUM_TOP_PREDICTIONS))
                        for bgc_fragment in bgc_variant.fragments
                        for bgc_module1, bgc_module2 in pairwise(bgc_fragment)))
    nrp_pairs = list(chain(*(pairwise([mon.residue for mon in fragment.monomers])
                             for fragment in nrp_fragments)))

    # check that lengths are close
    slope, intercept, margin = heuristic_matching_cfg.LINEAR_DISCARD_PARAMS_LENGTHS
    if abs(slope * nrp_len + intercept - bgc_len) > margin:
        return False

    # check that number of common residues is not too low
    num_common_residues = len(intersection_with_repeats(bgc_res, nrp_mons))
    slope, intercept, margin = heuristic_matching_cfg.LINEAR_DISCARD_PARAMS_AA_CONTENTS
    if num_common_residues < slope * nrp_len + intercept - margin:
        return False

    # check that number of common consecutive pairs is not too low
    num_common_pairs = len(intersection_with_repeats(nrp_pairs, bgc_pairs))
    slope, intercept, margin = heuristic_matching_cfg.LINEAR_DISCARD_PARAMS_CONS_PAIRS
    if num_common_pairs < slope * nrp_len + intercept - margin:
        return False
    return True


class DiscardVerdict(NamedTuple):
    joined: bool
    iterative: bool


def get_heuristic_discard(bgc_variant: BGC_Variant,
                          nrp_variant: NRP_Variant,
                          heuristic_matching_cfg: HeuristicMatchingConfig) -> DiscardVerdict:
        # num_cons_pairs = get_common_cons_pairs_cnt(match, bgc_variant)
        discard_joined = heuristic_match_discard(bgc_variant, nrp_variant.fragments, heuristic_matching_cfg)
        discard_it = any(heuristic_match_discard(bgc_variant,
                                                 [nrp_fragment],
                                                 heuristic_matching_cfg)
                         for nrp_fragment in nrp_variant.fragments)
        return DiscardVerdict(discard_joined, discard_it)


def filter_out_bad_nrp_candidates(bgc_variant: BGC_Variant,
                                  nrp_variants: Iterable[NRP_Variant],
                                  heuristic_matching_cfg: HeuristicMatchingConfig) -> Iterable[NRP_Variant]:
    bgc_res = bgc_residues(bgc_variant, heuristic_matching_cfg.NUM_TOP_PREDICTIONS)
    bgc_len = sum(len(bgc_fragment) for bgc_fragment in bgc_variant.fragments)

    def discard_variant(nrp_variant: NRP_Variant) -> bool:
        nrp_res = nrp_residues(nrp_variant)
        nrp_len = len(nrp_res)
        discard_nonit = heuristic_match_discard(bgc_len, nrp_len,
                                                len(intersection_with_repeats(bgc_res, nrp_res)),
                                                heuristic_matching_cfg)
        discard_it = any(heuristic_match_discard(bgc_len,
                                                 len(nrp_fragment.monomers),
                                                 len(intersection_with_repeats(bgc_res,
                                                                               [mon.residue for mon in nrp_fragment.monomers])),
                                                 heuristic_matching_cfg)
                            for nrp_fragment in nrp_variant.fragments)
        return discard_nonit and discard_it

    return filter(lambda nrp_variant: not discard_variant(nrp_variant), nrp_variants)
