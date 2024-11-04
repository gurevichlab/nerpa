import math

from typing import (
    Dict,
    List,
    Tuple,
    Optional
)
from src.data_types import (
    BGC_Module,
    BGC_Module_Modification,
    BGC_Fragment,
    Chirality,
    GeneId,
    NRP_Monomer,
    NRP_Fragment,
    LogProb,
    NRP_Monomer_Modification,
)
from src.antismash_parsing.location_features import (
    module_features_to_fragment_features,
    module_features_to_gene_features,
    ModuleLocFeatures
)
from src.rban_parsing.monomer_features import MonomerFeature, MonomerFeatures
from src.monomer_names_helper import UNKNOWN_RESIDUE, MonomerResidue

from src.matching.dp_types import DP_State
from src.matching.scoring_config import ScoringConfig, ModMatch, ChiralityMatch
from src.matching.matching_types_alignment_step import AlignmentStepType, MatchDetailedScore
from src.generic.other import get_score
from dataclasses import dataclass
from itertools import groupby
from enum import Enum, auto
from math import log

FUNCTION_NAME_TO_STEP_TYPE = {
    'bgc_module_skip': AlignmentStepType.BGC_MODULE_SKIP,
    'gene_skip': AlignmentStepType.GENE_SKIP,
    'bgc_fragment_skip': AlignmentStepType.BGC_FRAGMENT_SKIP,
    'nrp_mon_insert': AlignmentStepType.NRP_MONOMER_INSERT,
    'match': AlignmentStepType.MATCH,
    'iterate_module': AlignmentStepType.ITERATE_MODULE,
    'iterate_gene': AlignmentStepType.ITERATE_GENE,
    'skip_from_start': AlignmentStepType.SKIP_FROM_START,
    'skip_to_end': AlignmentStepType.SKIP_TO_END
}

@dataclass
class ScoringHelper:
    scoring_config: ScoringConfig
    heuristic_discard_on: bool = False  # immediately discard bad matches based on a heuristic
    bgc_modules: List[BGC_Module] = None
    nrp_monomers: List[NRP_Monomer] = None
    pks_domains_in_bgc: bool = False
    gene_pos_in_fragment: Dict[GeneId, int] = None
    genes_in_fragment: Dict[int, int] = None
    modules_in_gene: Dict[GeneId, int] = None

    def set_bgc_modules(self, bgc_modules: List[BGC_Module]):
        self.bgc_modules = bgc_modules
        self.gene_pos_in_fragment = {}
        self.genes_in_fragment = {}
        self.modules_in_gene = {}
        for fragment_idx, fragment_modules in groupby(bgc_modules, key=lambda module: module.fragment_idx):
            fragment_modules = list(fragment_modules)
            self.genes_in_fragment[fragment_idx] = len(set(module.gene_id for module in fragment_modules))
            for i, (gene_id, gene_modules) in enumerate(groupby(fragment_modules, key=lambda module: module.gene_id)):
                self.modules_in_gene[gene_id] = len(list(gene_modules))
                self.gene_pos_in_fragment[gene_id] = i




    def set_nrp_monomers(self, nrp_monomers: List[NRP_Monomer]):
        self.nrp_monomers = nrp_monomers

    def set_pks_domains_in_bgc(self, pks_domains_in_bgc: bool):
        self.pks_domains_in_bgc = pks_domains_in_bgc

    def bgc_module_skip(self, bgc_module_idx: int, dp_state: Optional[DP_State] = None) -> LogProb:  # argument is kept for better alignment backtracking
        module_features = self.bgc_modules[bgc_module_idx].module_loc
        return get_score(self.scoring_config.bgc_module_skip_score, module_features)

    def nrp_monomer_insert_score(self, module_features: ModuleLocFeatures,
                                 monomer_features: MonomerFeatures,
                                 start_of_bgc: bool) -> LogProb:
        if start_of_bgc:
            return get_score(self.scoring_config.nrp_monomer_insert_at_start_score, module_features)
        else:
            return get_score(self.scoring_config.nrp_monomer_insert_score, module_features)

    def nrp_mon_insert(self, mon_idx: int, dp_state: DP_State) -> LogProb:
        mon = self.nrp_monomers[mon_idx]
        start_of_bgc = dp_state.module_pos == 0
        bgc_module = self.bgc_modules[dp_state.module_pos - 1] if not start_of_bgc else self.bgc_modules[0]
        return self.nrp_monomer_insert_score(module_features=bgc_module.module_loc,
                                             monomer_features=mon.monomer_features,
                                             start_of_bgc=start_of_bgc)

    def match_residue_score(self, residue_scores: Dict[MonomerResidue, LogProb],
                            module_features: ModuleLocFeatures,
                            nrp_residue: MonomerResidue,
                            monomer_features: MonomerFeatures,
                            pks_domains_in_bgc: bool) -> LogProb:
        if MonomerFeature.PKS in monomer_features:
            return -math.inf  # PKS residues are never matched to BGC modules
        if MonomerFeature.PKS_HYBRID in monomer_features and not pks_domains_in_bgc:
             nrp_residue = UNKNOWN_RESIDUE  # if there are no PKS domains in the BGC, PKS hybrids are treated as unknown residues
        return residue_scores[nrp_residue]


    def match_modifications_score(self, bgc_modifications: Tuple[BGC_Module_Modification, ...],
                                  module_features: ModuleLocFeatures,
                                  nrp_modifications: Tuple[NRP_Monomer_Modification, ...],
                                  monomer_features: MonomerFeatures) -> LogProb:
        mod_match = ModMatch(NRP_Monomer_Modification.METHYLATION,
                             BGC_Module_Modification.METHYLATION in bgc_modifications,
                             NRP_Monomer_Modification.METHYLATION in nrp_modifications)
        return self.scoring_config.mod_score[mod_match]

    def match_chirality_score(self, bgc_modifications: Tuple[BGC_Module_Modification, ...],
                                module_features: ModuleLocFeatures,
                                nrp_chirality: Chirality,
                                monomer_features: MonomerFeatures) -> LogProb:
        if nrp_chirality == Chirality.UNKNOWN:
            return 0
        chirality_match = ChiralityMatch(BGC_Module_Modification.EPIMERIZATION in bgc_modifications,
                                         nrp_chirality)
        return self.scoring_config.chirality_score[chirality_match]

    def match_steptype_score(self, module_features: ModuleLocFeatures) -> LogProb:
        return 0  # placeholder for future scoring

    def match_detailed_score(self, bgc_module: BGC_Module,
                             nrp_monomer: NRP_Monomer) -> MatchDetailedScore:
        residue_score = self.match_residue_score(bgc_module.residue_score,
                                                 bgc_module.module_loc,
                                                 nrp_monomer.residue,
                                                 nrp_monomer.monomer_features,
                                                 self.pks_domains_in_bgc)
        modifications_score = self.match_modifications_score(bgc_module.modifications,
                                                             bgc_module.module_loc,
                                                             nrp_monomer.modifications,
                                                             nrp_monomer.monomer_features)
        chirality_score = self.match_chirality_score(bgc_module.modifications,
                                                     bgc_module.module_loc,
                                                     nrp_monomer.chirality,
                                                     nrp_monomer.monomer_features)

        match_steptype_score = self.match_steptype_score(bgc_module.module_loc)
        return MatchDetailedScore(residue_score,
                                  modifications_score,
                                  chirality_score)

    def match(self, bgc_module_idx: int, nrp_monomer_idx: int, dp_state: Optional[DP_State] = None) -> LogProb:
        return sum(self.match_detailed_score(self.bgc_modules[bgc_module_idx],self.nrp_monomers[nrp_monomer_idx]))

    def iterate_module(self, dp_state: Optional[DP_State] = None) -> LogProb:
        return -0.5  # placeholder for future scoring

    def iterate_gene(self, dp_state: Optional[DP_State] = None) -> LogProb:
        return -0.5

    def null_hypothesis_score(self, nrp_monomer: NRP_Monomer) -> LogProb:
        return self.scoring_config.null_hypothesis_residue_score[nrp_monomer.residue] \
            + self.scoring_config.null_hypothesis_chirality_score[nrp_monomer.chirality] \
            + self.scoring_config.null_hypothesis_mod_score[(NRP_Monomer_Modification.METHYLATION,
                                                             NRP_Monomer_Modification.METHYLATION in nrp_monomer.modifications)]

    def gene_skip(self, fst_module_idx: int, last_module_idx: int,
                  dp_state: Optional[DP_State] = None) -> LogProb:
        gene_features = module_features_to_gene_features(self.bgc_modules[fst_module_idx].module_loc,
                                                         self.bgc_modules[last_module_idx].module_loc)
        return get_score(self.scoring_config.gene_skip_score, gene_features)


    def bgc_fragment_skip(self, fst_module_idx: int, last_module_idx: int,
                          dp_state: Optional[DP_State] = None) -> LogProb:
            fragment_features = module_features_to_fragment_features(self.bgc_modules[fst_module_idx].module_loc,
                                                                    self.bgc_modules[last_module_idx].module_loc)
            return get_score(self.scoring_config.bgc_fragment_skip_score, fragment_features)

    # skip all modules [0:i) and monomers [0:j) from the start
    # TODO: account for context
    def skip_from_start(self,
                        module_idx: int,
                        monomer_idx: int,
                        dp_state: DP_State) -> LogProb:
        module = self.bgc_modules[module_idx]
        fragment_skips = module.fragment_idx
        gene_skips = self.gene_pos_in_fragment[module.gene_id]
        module_skips = module.a_domain_idx
        return -3 + monomer_idx * get_score(self.scoring_config.nrp_monomer_insert_at_start_score, module.module_loc)
        '''
        return (fragment_skips * log(0.1)  # TODO: put in config
                + gene_skips * log(0.05)
                + module_skips * log(0.03)
                + monomer_idx * get_score(self.scoring_config.nrp_monomer_insert_at_start_score, module.module_loc)
                )
        '''

    # skip modules [i:] and monomers [j:]
    def skip_to_end(self,
                    pi: int,
                    pj: int,
                    dp_state: Optional[DP_State] = None) -> LogProb:
        module = self.bgc_modules[pi - 1]  # last matched module
        fragment_skips = self.bgc_modules[-1].fragment_idx - module.fragment_idx
        gene_skips = self.genes_in_fragment[module.fragment_idx] - self.gene_pos_in_fragment[module.gene_id]
        module_skips = self.modules_in_gene[module.gene_id] - module.a_domain_idx
        monomer_skips = len(self.nrp_monomers) - pj
        return -3 + monomer_skips * get_score(self.scoring_config.nrp_monomer_insert_at_start_score, module.module_loc)  # should be insert at end score
        '''
        return (fragment_skips * log(0.1)  # TODO: put in config
                + gene_skips * log(0.05)
                + module_skips * log(0.03)
                + monomer_skips * get_score(self.scoring_config.nrp_monomer_insert_at_start_score, module.module_loc)  # should be insert at end score
                )
        '''