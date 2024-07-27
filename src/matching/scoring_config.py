from typing import (
    Callable,
    Dict,
    Literal,
    List,
    NamedTuple,
    Tuple
)
from src.data_types import (
    Chirality,
    LogProb,
    MonomerResidue,
    NRP_Monomer_Modification,
    UNKNOWN_RESIDUE
)
from src.matching.alignment_types import AlignmentStepType
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import yaml
from math import log


class ChiralityMatch(NamedTuple):
    bgc_epim: bool
    nrp_chr: Chirality


class ModMatch(NamedTuple):
    mod: NRP_Monomer_Modification
    bgc_mod: bool
    nrp_mod: bool


@dataclass
class ScoringConfig:
    match_score: LogProb
    bgc_module_skip_score: LogProb
    nrp_monomer_skip_score: Dict[MonomerResidue, LogProb]

    num_unknown_residues: int
    mod_score: Dict[ModMatch, LogProb]
    chirality_score: Dict[ChiralityMatch, LogProb]

    null_hypothesis_residue_score: Dict[MonomerResidue, LogProb]
    null_hypothesis_mod_score: Dict[Tuple[NRP_Monomer_Modification, bool], LogProb]
    null_hypothesis_chirality_score: Dict[Chirality, LogProb]

    max_gene_reps: int
    max_module_reps: int

    normalization: Literal['AVERAGING', 'RANDOM_NULL_MODEL', 'NONE']

    pks_residues: List[str]


def load_modificatons_score(cfg: dict) -> Dict[ModMatch, LogProb]:
    return {ModMatch(mod=mod, bgc_mod=bgc_mod, nrp_mod=nrp_mod):
                cfg['modification_score'][f'Mod_{mod.name}_BGC_{bgc_mod}_NRP_{nrp_mod}']
            for mod in NRP_Monomer_Modification
            for bgc_mod in (False, True)
            for nrp_mod in (False, True)
            if mod != NRP_Monomer_Modification.UNKNOWN}

def load_chirality_score(cfg: dict) -> Dict[ChiralityMatch, LogProb]:
    return {ChiralityMatch(bgc_epim=bgc_epim, nrp_chr=nrp_chr):
                cfg['chirality_score'][f'BGC_{bgc_epim}_NRP_{nrp_chr.name}']
            for bgc_epim in (False, True)
            for nrp_chr in Chirality}

def load_nrp_monomer_skip_score(cfg: dict) -> Dict[MonomerResidue, LogProb]:
    nrp_mon_skip_log_prob = log(cfg['alignment_step_frequency']['NRP_MONOMER_SKIP'])
    '''
    return defaultdict(lambda: nrp_mon_skip_log_prob + cfg['nrp_monomer_frequencies'][UNKNOWN_RESIDUE],
                       {residue: nrp_mon_skip_log_prob + residue_frequency
                        for residue, residue_frequency in cfg['nrp_monomer_frequencies'].items()})
    '''
    return defaultdict(lambda: nrp_mon_skip_log_prob)

def load_null_hypothesis_residue_score(cfg: dict) -> Dict[MonomerResidue, LogProb]:
    return defaultdict(lambda: cfg['null_hypothesis_residue_score'][UNKNOWN_RESIDUE],
                       cfg['null_hypothesis_residue_score'])

def load_null_hypothesis_mod_score(cfg: dict) -> Dict[Tuple[NRP_Monomer_Modification, bool], LogProb]:
    mod_score = load_modificatons_score(cfg)

    def get_null_hyp_mod_score(mod: NRP_Monomer_Modification,
                               mod_in_nrp: bool) -> LogProb:
        mod_freq = cfg['module_modification_frequency'][mod.name]
        return mod_freq * mod_score[ModMatch(mod=mod, bgc_mod=True, nrp_mod=mod_in_nrp)] \
            + (1-mod_freq) * mod_score[ModMatch(mod=mod, bgc_mod=False, nrp_mod=mod_in_nrp)]

    return {(mod, mod_in_nrp): get_null_hyp_mod_score(mod, mod_in_nrp)
            for mod in NRP_Monomer_Modification
            for mod_in_nrp in (False, True)
            if mod != NRP_Monomer_Modification.UNKNOWN}

def load_null_hypothesis_chirality_score(cfg: dict) -> Dict[Chirality, LogProb]:
    chirality_score = load_chirality_score(cfg)

    def get_null_hyp_chr_score(nrp_chr: Chirality) -> LogProb:
        epim_freq = cfg['module_epimerization_frequency']
        return epim_freq * chirality_score[ChiralityMatch(bgc_epim=True, nrp_chr=nrp_chr)] \
            + (1-epim_freq) * chirality_score[ChiralityMatch(bgc_epim=False, nrp_chr=nrp_chr)]

    return {chr: get_null_hyp_chr_score(chr)
            for chr in Chirality}


def load_scoring_config(path_to_config: Path) -> ScoringConfig:
    cfg = yaml.safe_load(path_to_config.open('r'))

    nrp_monomer_skip_score = load_nrp_monomer_skip_score(cfg)
    mod_score = load_modificatons_score(cfg)
    chirality_score = load_chirality_score(cfg)

    null_hypothesis_residue_score = load_null_hypothesis_residue_score(cfg)
    null_hypothesis_mod_score = load_null_hypothesis_mod_score(cfg)
    null_hypothesis_chirality_score = load_null_hypothesis_chirality_score(cfg)

    return ScoringConfig(match_score=log(cfg['alignment_step_frequency']['MATCH']),
                         bgc_module_skip_score=log(cfg['alignment_step_frequency']['BGC_MODULE_SKIP']),
                         nrp_monomer_skip_score=nrp_monomer_skip_score,
                         num_unknown_residues=cfg['num_unknown_residues'],
                         mod_score=mod_score,
                         chirality_score=chirality_score,
                         null_hypothesis_residue_score=null_hypothesis_residue_score,
                         null_hypothesis_mod_score=null_hypothesis_mod_score,
                         null_hypothesis_chirality_score=null_hypothesis_chirality_score,
                         max_gene_reps=cfg['max_gene_repetitions'],
                         max_module_reps=cfg['max_module_repetitions'],
                         normalization=cfg['normalization'],
                         pks_residues=cfg['pks_residues'])
