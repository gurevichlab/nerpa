from __future__ import annotations
from typing import List
from dataclasses import dataclass
from src.data_types import LogProb
from src.matching.match_type import Match_BGC_Variant_Info, Match_NRP_Variant_Info, Match
from src.rban_parsing.get_linearizations import LinearizationLight
from src.matching.detailed_hmm import DetailedHMM
from src.data_types import NRP_Variant
from src.rban_parsing.rban_monomer import rBAN_idx
from collections import defaultdict

@dataclass
class HMM_Match:
    score: LogProb
    bgc_variant_info: Match_BGC_Variant_Info
    nrp_id: str
    nrp_linearizations: List[List[rBAN_idx]]
    optimal_paths: List[List[int]]

    @classmethod
    def from_json(cls, json_data: dict) -> HMM_Match:
        return cls(score=json_data['score'],
                   bgc_variant_info=Match_BGC_Variant_Info(**json_data['bgc_variant_info']),
                   nrp_id=json_data['nrp_id'],
                   nrp_linearizations=json_data['nrp_linearizations'],
                   optimal_paths=json_data['optimal_paths'])


def convert_to_detailed_matches(hmms: List[DetailedHMM],
                                nrp_variants: List[NRP_Variant],
                                hmm_matches: List[HMM_Match]) -> List[Match]:
    hmm_by_bgc_info = {}
    for hmm in hmms:
        bgc_variant_info = Match_BGC_Variant_Info.from_bgc_variant(hmm.bgc_variant)
        hmm_by_bgc_info[bgc_variant_info] = hmm

    rban_idx_to_rban_mon = defaultdict(dict)
    for nrp_variant in nrp_variants:
        for fragment in nrp_variant.fragments:
            for rban_mon in fragment.monomers:
                rban_idx_to_rban_mon[nrp_variant.nrp_id][rban_mon.rban_idx] = rban_mon

    matches = []
    for hmm_match in hmm_matches:
        hmm = hmm_by_bgc_info[hmm_match.bgc_variant_info]
        alignments = []
        for nrp_linearization, optimal_path in zip(hmm_match.nrp_linearizations, hmm_match.optimal_paths):
            rban_monomers = [rban_idx_to_rban_mon[hmm_match.nrp_id][rban_idx]
                             for rban_idx in nrp_linearization]
            alignments.append(hmm.path_to_alignment(optimal_path, rban_monomers))

        matches.append(Match(bgc_variant_info=hmm_match.bgc_variant_info,
                             nrp_variant_info=Match_NRP_Variant_Info(nrp_id=hmm_match.nrp_id, variant_idx=0),
                             normalized_score=hmm_match.score,
                             alignments=alignments))

    return matches