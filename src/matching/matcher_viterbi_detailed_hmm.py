from __future__ import annotations
from typing import (
    ClassVar,
    Dict,
    List,
    NamedTuple,
    Tuple,
    Optional
)
from dataclasses import dataclass
from enum import Enum, auto

from src.antismash_parsing.location_features import (
    ModuleLocFeature,
    GeneLocFeature,
    BGC_Fragment_Loc_Feature
)
from src.data_types import NRP_Monomer, BGC_Module, BGC_Variant, GeneId, LogProb
from src.matching.matcher_viterbi_types import (
    DetailedHMMEdgeType,
    DetailedHMMStateType,
    DetailedHMMState,
    DetailedHMMEdge,
    HMM
)
from src.matching.matcher_viterbi_algorithm import get_opt_path_with_emissions
from src.matching.alignment_to_path_in_hmm import alignment_to_hmm_path
from src.matching.hmm_edge_weights import EdgeWeightsParams
from src.build_output.draw_hmm import draw_hmm
from src.monomer_names_helper import MonomerNamesHelper
from src.rban_parsing.rban_monomer import rBAN_Monomer
from collections import defaultdict

from src.matching.bgc_to_hmm import bgc_variant_to_detailed_hmm
from src.matching.hmm_to_alignment import hmm_path_to_alignment
from src.matching.matching_types_alignment import Alignment
from src.matching.hmm_anchors_heuristic import heuristic_opt_path
from src.matching.hmm_scoring_helper import HMMHelper
from itertools import pairwise
from graphviz import Digraph
from pathlib import Path


@dataclass
class DetailedHMM:
    states: List[DetailedHMMState]
    adj_list: List[Dict[int, DetailedHMMEdge]]
    start_state_idx: int
    final_state_idx: int
    bgc_variant: BGC_Variant
    state_idx_to_module_idx: Dict[int, int]
    _module_idx_to_state_idx: List[int]  # points to MODULE_START state for each module. Used for building hmm from alignment

    hmm_helper: ClassVar[Optional[HMMHelper]] = None

    @classmethod
    def from_bgc_variant(cls, bgc_variant: BGC_Variant) -> DetailedHMM:
        if cls.hmm_helper is None:
            raise ValueError("HMM helper must be set before calling DetailedHMM.from_bgc_variant")
        return bgc_variant_to_detailed_hmm(DetailedHMM, bgc_variant)

    def to_hmm(self) -> HMM:
        num_states = len(self.states)

        adj_list = [[(edge_to, edge_data.log_prob)
                     for edge_to, edge_data in self.adj_list[edge_from].items()]
                    for edge_from in range(num_states)]

        emission_log_probs = [[state.emissions[mon]
                               for mon in sorted(state.emissions,
                                                 key=lambda m: self.hmm_helper.monomer_names_helper.mon_to_int[m])]  # TODO: define int(m) or use monomer_names_helper
                              for state in self.states]
        return HMM(adj_list=adj_list,
                   emission_log_probs=emission_log_probs)

    def get_opt_path_with_emissions(self,
                                    start: int,
                                    finish: int,
                                    emitted_monomers: List[rBAN_Monomer]) -> List[Tuple[int, Optional[rBAN_Monomer]]]:
        monomer_codes = [self.hmm_helper.monomer_names_helper.mon_to_int[mon.to_base_mon()] for mon in emitted_monomers]
        score, path_with_emissions = get_opt_path_with_emissions(self.to_hmm(),
                                                                 start, finish,
                                                                 monomer_codes)
        emitted_monomers_iter = iter(emitted_monomers)
        return [(u, None if mon_int is None else next(emitted_monomers_iter))
                for u, mon_int in path_with_emissions]

    def path_to_alignment(self,
                          path: List[int],
                          nrp_monomers: List[rBAN_Monomer]) -> Alignment:
        return hmm_path_to_alignment(self, path, nrp_monomers)

    def alignment_to_path_with_emisions(self, alignment: Alignment) -> List[Tuple[int, Optional[NRP_Monomer]]]:
        return alignment_to_hmm_path(self, alignment)

    def get_alignment(self, emitted_monomers: List[rBAN_Monomer],
                      heuristic_matching: bool = True) -> Alignment:
        if heuristic_matching:
            bgc_predictions = [self.hmm_helper.get_emissions(module, self.bgc_variant.has_pks_domains())
                               for module in self.bgc_variant.modules]
            opt_path = heuristic_opt_path(self, bgc_predictions, emitted_monomers)
        else:
            opt_path, _ = self.get_opt_path_with_emissions(self.start_state_idx,
                                                           self.final_state_idx,
                                                           emitted_monomers)
        #self.draw(Path(f'{self.bgc_variant.genome_id}.png'), opt_path)  # for debugging
        #print('opt_path', opt_path)
        return self.path_to_alignment(opt_path, emitted_monomers)


    def draw(self,
             filename: Path,
             highlight_path: Optional[List[int]] = None) -> Digraph:
        return draw_hmm(self, filename, highlight_path)


# captures all context relevant for parameter training
class EdgeFingerprint(NamedTuple):
    edge_type: DetailedHMMEdgeType
    bgc_module_context: Optional[ModuleLocFeature] = None
    bgc_gene_context: Optional[GeneLocFeature] = None
    bgc_fragment_context: Optional[BGC_Fragment_Loc_Feature] = None

def edge_fingerprint(hmm: DetailedHMM,
                     path: List[Tuple[int, NRP_Monomer]],
                     edge_idx: int) -> EdgeFingerprint:
    u, v = path[edge_idx][0], path[edge_idx + 1][0]
    edge = hmm.adj_list[u][v]

    # find corresponding bgc module
    match edge.edge_type:
        # genomic info not relevant
        case (
        # before matching
            DetailedHMMEdgeType.START_INSERTING_AT_START |
            DetailedHMMEdgeType.INSERT_AT_START |
            DetailedHMMEdgeType.START_SKIP_MODULES_AT_START |
            DetailedHMMEdgeType.START_SKIP_GENES_AT_START |
            DetailedHMMEdgeType.START_SKIP_FRAGMENTS_AT_START |
        # at matching
            DetailedHMMEdgeType.INSERT |
            DetailedHMMEdgeType.END_INSERTING
        ):
            return EdgeFingerprint(edge_type=edge.edge_type)
        # only fragment context relevant
        case (
            DetailedHMMEdgeType.SKIP_FRAGMENT_AT_START |
            DetailedHMMEdgeType.SKIP_FRAGMENT |
            DetailedHMMEdgeType.SKIP_FRAGMENT_AT_END
        ):
            return EdgeFingerprint(edge_type=edge.edge_type,
                                   bgc_fragment_context=hmm.bgc_variant.fragments[hmm.state_idx_to_module_idx[u]].loc)


