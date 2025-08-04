import math
from collections import defaultdict
from math import log
from typing import (
    Any,
    Dict,
    Literal,
    List,
    Optional,
    NamedTuple,
    NewType,
    Tuple
)
import pandas as pd
from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path

from scripts.extract_norine_monomers_info import norine_graphs_path
from src.antismash_parsing.antismash_name_mappings import KNOWN_SUBSTRATES
from src.generic.functional import cached_by_key
import yaml


AA10 = NewType('AA10', str)  # 10-letter Stachelhaus code
AA34 = NewType('AA34', str)  # 34-letter code
antiSMASH_MonomerName = NewType('antiSMASH_MonomerName', str)  # antiSMASH Short
NorineMonomerName = NewType('NorineMonomerName', str)
MonomerResidue = NewType('MonomerResidue', str)
UNKNOWN_RESIDUE = MonomerResidue('unknown')
PKS = MonomerResidue('PKS')
MonCode = NewType('MonCode', int)

def enum_representer(dumper, e: Enum):
    return dumper.represent_scalar(f'!{e.__class__.__name__}', e.name)

class Chirality(Enum):
    L = auto()
    D = auto()
    UNKNOWN = auto()


class NRP_Monomer_Modification(Enum):  # post-translational modification
    METHYLATION = auto()
    UNKNOWN = auto()

yaml.add_representer(Chirality, enum_representer)
yaml.add_representer(NRP_Monomer_Modification, enum_representer)

@dataclass
class NRP_Monomer:
    residue: MonomerResidue
    methylated: bool = False
    chirality: Chirality = Chirality.UNKNOWN
    is_pks_hybrid: bool = False

    def __hash__(self):
        return hash((self.residue, self.methylated, self.chirality, self.is_pks_hybrid))

UNKNOWN_MONOMER = NRP_Monomer(residue=UNKNOWN_RESIDUE,
                              methylated=False,  # TODO: should we have a separate UNKNOWN_METHYLATION?,
                              chirality=Chirality.UNKNOWN)
PKS_MONOMER = NRP_Monomer(residue=PKS,
                          methylated=False,  # TODO: should we have a separate UNKNOWN_METHYLATION?,
                          chirality=Chirality.UNKNOWN)

Prob = NewType('Prob', float)  # TODO: import

class MonomersDefaultFrequencies(NamedTuple):
    residue: Dict[MonomerResidue, Prob]
    methylation: Prob
    d_chirality: Prob


@dataclass
class MonomerNamesHelper:
    names_table: pd.DataFrame
    supported_residues: List[MonomerResidue]
    pks_names: List[NorineMonomerName]
    norine_monomer_cnts: Dict[NorineMonomerName, int]

    _cache: Dict[Tuple[str, str], NRP_Monomer] = None
    default_frequencies: MonomersDefaultFrequencies = None

    mon_to_int: Dict[NRP_Monomer, MonCode] = None
    int_to_mon: Dict[MonCode, NRP_Monomer] = None


    def __post_init__(self):
        self._cache = {}

        # Validation
        assert all(res in self.names_table['core'].values or res == UNKNOWN_RESIDUE
                   for res in self.supported_residues), \
            f"Some supported residues are absent in the names table: " \
            f"{set(self.supported_residues) - set(self.names_table['core'].values)}"
        residues_from_norine = {self.parsed_name(norine_name, 'norine').residue
                                for norine_name in self.norine_monomer_cnts}
        assert set(self.supported_residues) == residues_from_norine, \
            f"Supported residues {set(self.supported_residues) - residues_from_norine} are missing in norine"

        self._set_default_frequencies()

        self.int_to_mon = {}
        self.mon_to_int = {}
        for mon_res_int, mon_res in enumerate(self.supported_residues):
            for meth_int, methylated in enumerate([False, True]):
                for chir_int, chirality in enumerate([Chirality.L, Chirality.D, Chirality.UNKNOWN]):
                    for is_pks_hybrid_int, is_pks_hybrid in enumerate([False, True]):
                        mon = NRP_Monomer(residue=MonomerResidue(mon_res),
                                          methylated=methylated,
                                          chirality=chirality,
                                          is_pks_hybrid=is_pks_hybrid)
                        mon_int = mon_res_int * 12 + meth_int * 6 + chir_int * 2 + is_pks_hybrid_int
                        self.mon_to_int[mon] = MonCode(mon_int)
                        self.int_to_mon[MonCode(mon_int)] = mon

    def _set_default_frequencies(self):
        """
        Set the default frequencies for residues, methylation, and D-chirality.
        """
        parsed_monomers_cnts: Dict[NRP_Monomer, int] = {}
        for monomer_name, cnt in self.norine_monomer_cnts.items():
            parsed_monomer = self.parsed_name(monomer_name, name_format='norine')
            parsed_monomers_cnts[parsed_monomer] = parsed_monomers_cnts.get(parsed_monomer, 0) + cnt

        total_monomers = sum(parsed_monomers_cnts.values())
        total_methylated = sum(cnt for monomer, cnt in parsed_monomers_cnts.items()
                               if monomer.methylated)
        total_d_chirality = sum(cnt for monomer, cnt in parsed_monomers_cnts.items()
                                if monomer.chirality == Chirality.D)
        norine_methylation_freq = Prob(total_methylated / total_monomers)
        norine_d_chirality_freq = Prob(total_d_chirality / total_monomers)

        nerpa_residue_cnts: Dict[MonomerResidue, int] = {}
        for mon, cnt in parsed_monomers_cnts.items():
            nerpa_residue_cnts[mon.residue] = nerpa_residue_cnts.get(mon.residue, 0) + cnt

        residue_freqs = {residue: Prob(cnt / total_monomers)
                         for residue, cnt in nerpa_residue_cnts.items()}

        self.default_frequencies = MonomersDefaultFrequencies(residue=residue_freqs,
                                                              methylation=norine_methylation_freq,
                                                              d_chirality=norine_d_chirality_freq)

    @classmethod
    def parsed_modifications(cls, mods_str: str) -> Tuple[NRP_Monomer_Modification, ...]:
        def parse_mod(mod: str) -> NRP_Monomer_Modification:
            match mod:
                case 'MT': return NRP_Monomer_Modification.METHYLATION
                case other: return NRP_Monomer_Modification.UNKNOWN

        mods_list = eval(mods_str)
        return () if mods_list == [''] else tuple(map(parse_mod, mods_list))

    def parsed_name(self, name: str,
                    name_format: Literal['antismash', 'norine', 'paras']) -> NRP_Monomer:
        if (name, name_format) in self._cache:
            return self._cache[(name, name_format)]

        if name_format == 'paras':
            return NRP_Monomer(residue=paras_residue_to_nerpa_residue(name, self),
                               methylated=False,
                               chirality=Chirality.UNKNOWN)

        if name == 'Ile/aIle':
            name = 'Ile'  # temporary fix for the table
        if name == 'NMe-aIle/NMe-Ile':
            name = 'NMe-Ile'
        if name == 'aThr/Thr':
            name = 'Thr'
        if name == 'NMe-aThr/NMe-Thr':
            name = 'NMe-Thr'
        if name == 'NMe-bMe-Leu/diMe-aIle':
            name = 'NMe-Leu'

        chirality = Chirality.UNKNOWN
        if name_format == 'norine':
            if name.startswith('X'):
                return UNKNOWN_MONOMER
            if name in self.pks_names:
                return PKS_MONOMER
            if name.startswith('D-'):
                chirality = Chirality.D
                name = name[2:]
            else:
                chirality = Chirality.L

        column_name = 'as_short' if name_format == 'antismash' else 'norine_rban_code'

        try:
            row = self.names_table[self.names_table[column_name] == name].iloc[0]
            residue = MonomerResidue(row['core']) if row['core'] in self.supported_residues else UNKNOWN_RESIDUE
        except IndexError:
            #print(f'WARNING: Name {name} not found in the names table')
            return UNKNOWN_MONOMER  # TODO: raise ValueError when the table is ready
            # raise ValueError(f'Name {name} not found in the names table')

        self._cache[(name, name_format)] = NRP_Monomer(residue=residue,
                                                       methylated=NRP_Monomer_Modification.METHYLATION in self.parsed_modifications(row['modifications']),
                                                       chirality=chirality)
        return self._cache[(name, name_format)]


PARAS_RESIDUE = str

def paras_residue_to_nerpa_residue(paras_name: PARAS_RESIDUE,
                                   monomer_names_helper: MonomerNamesHelper) -> MonomerResidue:
    core_map = {
        '2,3-dihydroxybenzoic acid': 'Bza',
        'anthranilic acid': 'Bza',
        'salicylic acid': 'Bza',
        '2,4-diaminobutyric acid': 'Dab',
        '2-aminoadipic acid': 'Aad',
        '2-aminoisobutyric acid': 'Aib',
        '3,5-dihydroxyphenylglycine': 'dHpg',
        '4-hydroxyphenylglycine': 'Hpg',
        'D-alanine': 'Ala',
        'alanine': 'Ala',
        'beta-alanine': 'bAla',
        'N5-hydroxyornithine': 'Orn',
        'N5-formyl-N5-hydroxyornithine': 'Orn',
        'R-beta-hydroxytyrosine': 'Tyr',
        'arginine': 'Arg',
        'asparagine': 'Asn',
        'aspartic acid': 'Asp',
        'cysteine': 'Cys',
        'glutamic acid': 'Glu',
        'glutamine': 'Gln',
        'glycine': 'Gly',
        'histidine': 'His',
        'isoleucine': 'Ile',
        'leucine': 'Leu',
        'lysine': 'Lys',
        'ornithine': 'Orn',
        'phenylalanine': 'Phe',
        'pipecolic acid': 'Pip',
        'proline': 'Pro',
        'serine': 'Ser',
        'threonine': 'Thr',
        'tryptophan': 'Trp',
        'tyrosine': 'Tyr',
        'valine': 'Val',
        # anything not in this dict can be mapped to 'unknown'
    }
    return core_map.get(paras_name, UNKNOWN_RESIDUE)

"""
    paras_name_core = paras_name.split('-')[-1]
    try:
        as_substrate = min((substrate
                       for substrate in KNOWN_SUBSTRATES
                       if paras_name_core in substrate.long),
                       key=lambda substrate: len(substrate.long))
        as_short = as_substrate.short
        result = monomer_names_helper.parsed_name(as_short, name_format='antismash').residue
    except StopIteration:
        print(f"WARNING: paras name {paras_name} not recognized")
        result = UNKNOWN_RESIDUE

    return result
"""