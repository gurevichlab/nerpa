from typing import Iterable, List, Dict
from itertools import chain
from collections import defaultdict
from src.data_types import MonomerResidue, LogProb, BGC_Variant, NRP_Variant
from src.generic.combinatorics import intersection_with_repeats
from src.matching.scoring_config import HeuristicMatchingConfig


def bgc_residues(bgc_variant: BGC_Variant,
                 num_top_predictions: int = 3) -> List[MonomerResidue]:
    def top_n_predictions(residue_scores: Dict[MonomerResidue, LogProb]) -> List[MonomerResidue]:
        return sorted(residue_scores, key=residue_scores.get, reverse=True)[:num_top_predictions]


    return list(chain(*[top_n_predictions(bgc_module.residue_score)
                        for bgc_fragment in bgc_variant.fragments
                        for bgc_module in bgc_fragment]))


def nrp_residues(nrp_variant: NRP_Variant) -> List[MonomerResidue]:
    return [mon.residue for fragment in nrp_variant.fragments
            for mon in fragment.monomers]


def heuristic_match_discard(bgc_len: int,
                            nrp_len: int,
                            num_common_residues: int,
                            heuristic_matching_cfg: HeuristicMatchingConfig) -> bool:
    slope, intercept, margin = heuristic_matching_cfg.LINEAR_DISCARD_PARAMS_LENGTHS
    if abs(slope * nrp_len + intercept  - bgc_len) > margin:
        return False

    slope, intercept, margin = heuristic_matching_cfg.LINEAR_DISCARD_PARAMS_AA_CONTENTS
    if abs(slope * nrp_len + intercept - num_common_residues) > margin:
        return False

    return True


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
        result = discard_nonit and discard_it
        if result:  # for debugging
            print(f'Discarding {nrp_variant.nrp_id} due to heuristic matching')
        return result

    return filter(lambda nrp_variant: not discard_variant(nrp_variant), nrp_variants)
