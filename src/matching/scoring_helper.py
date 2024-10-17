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
from src.matching.alignment_types import AlignmentStepType
from src.generic.other import get_score
from dataclasses import dataclass
from enum import Enum, auto

FUNCTION_NAME_TO_STEP_TYPE = {'bgc_module_skip': AlignmentStepType.BGC_MODULE_SKIP,
                              'nrp_mon_skip': AlignmentStepType.NRP_MONOMER_SKIP,
                              'match': AlignmentStepType.MATCH,
                              'iterate_module': AlignmentStepType.ITERATE_MODULE,
                              'iterate_gene': AlignmentStepType.ITERATE_GENE}

@dataclass
class ScoringHelper:
    scoring_config: ScoringConfig
    heuristic_discard_on: bool = False  # immediately discard bad matches based on a heuristic
    bgc_modules: List[BGC_Module] = None
    nrp_monomers: List[NRP_Monomer] = None
    pks_domains_in_bgc: bool = False

    def set_bgc_modules(self, bgc_modules: List[BGC_Module]):
        self.bgc_modules = bgc_modules

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
        chirality_match = ChiralityMatch(BGC_Module_Modification.EPIMERIZATION in bgc_modifications,
                                        nrp_chirality)
        return self.scoring_config.chirality_score[chirality_match]

    def match_steptype_score(self, module_features: ModuleLocFeatures) -> LogProb:
        return 0  # placeholder for future scoring

    def match_detailed_score(self, bgc_module: BGC_Module,
                             nrp_monomer: NRP_Monomer) -> Tuple[LogProb, LogProb, LogProb, LogProb]:
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
        return (residue_score,
                modifications_score,
                chirality_score,
                match_steptype_score)

    def match(self, bgc_module_idx: int, nrp_monomer_idx: int, dp_state: Optional[DP_State] = None) -> LogProb:
        return sum(self.match_detailed_score(self.bgc_modules[bgc_module_idx],self.nrp_monomers[nrp_monomer_idx]))

    def iterate_module(self, dp_state: Optional[DP_State] = None) -> LogProb:
        return 0

    def iterate_gene(self, dp_state: Optional[DP_State] = None) -> LogProb:
        return 0

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


    def fragment_skip(self, fst_module_idx: int, last_module_idx: int,
                      dp_state: Optional[DP_State] = None) -> LogProb:
            fragment_features = module_features_to_fragment_features(self.bgc_modules[fst_module_idx].module_loc,
                                                                    self.bgc_modules[last_module_idx].module_loc)
            return get_score(self.scoring_config.fragment_skip_score, fragment_features)

    def skip_bgc_fragment_score(self, bgc_fragment: BGC_Fragment) -> LogProb:
        return self.scoring_config.bgc_fragment_skip_penalty * len(bgc_fragment)

    def skip_nrp_fragment_score(self, nrp_fragment: NRP_Fragment) -> LogProb:
        return self.scoring_config.nrp_fragment_skip_penalty * len(nrp_fragment.monomers)
