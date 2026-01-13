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
    Tuple, Iterable
)
from src.build_output.yaml_representers import enum_representer
from src.general_type_aliases import Prob
import pandas as pd
import polars as pl
from dataclasses import dataclass
from enum import Enum, auto

import yaml

from src.pipeline.logging.logger import NerpaLogger

antiSMASH_MonomerName = NewType('antiSMASH_MonomerName', str)  # antiSMASH Short
NorineMonomerName = NewType('NorineMonomerName', str)
PARAS_MonomerName = NewType('PARAS_MonomerName', str)
NerpaResidue = NewType('NerpaResidue', str)

# Special "residues"
UNKNOWN_RESIDUE = NerpaResidue('_UNKNOWN')
NOT_NRPS_RESIDUE = NerpaResidue('_NOT_NRPS')  # residue that can't be attracted by an A domain, e.g. PKS
PKS_RESIDUE = NerpaResidue('_PKS')  # PKS residue in a PKS-NRPS hybrid
SPECIAL_RESIDUES = [UNKNOWN_RESIDUE, NOT_NRPS_RESIDUE, PKS_RESIDUE]

MonCode = NewType('MonCode', int)

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
    residue: NerpaResidue
    methylated: bool = False
    chirality: Chirality = Chirality.UNKNOWN
    is_pks_hybrid: bool = False

    def __hash__(self):
        return hash((self.residue, self.methylated, self.chirality, self.is_pks_hybrid))

    def __lt__(self, other):
        if other is None:
            return False
        return (self.residue, self.methylated, self.chirality.name, self.is_pks_hybrid) < \
               (other.residue, other.methylated, other.chirality.name, other.is_pks_hybrid)


UNKNOWN_MONOMER = NRP_Monomer(residue=UNKNOWN_RESIDUE,
                              methylated=False,  # TODO: should we have a separate UNKNOWN_METHYLATION?,
                              chirality=Chirality.UNKNOWN)
PKS_MONOMER = NRP_Monomer(residue=PKS_RESIDUE,
                          methylated=False,  # TODO: should we have a separate UNKNOWN_METHYLATION?,
                          chirality=Chirality.UNKNOWN)
NOT_NRPS_MONOMER = NRP_Monomer(residue=NOT_NRPS_RESIDUE,
                               methylated=False,  # TODO: should we have a separate UNKNOWN_METHYLATION?,
                               chirality=Chirality.UNKNOWN)


class MonomersDefaultFrequencies(NamedTuple):
    residue: Dict[NerpaResidue, Prob]
    methylation: Prob
    d_chirality: Prob


@dataclass
class MonomerNamesHelper:
    '''
    names_table: a dataframe with columns
     MonomerName,NameFormat,antiSMASH_short_core/Nerpa_residue,Modifications,Type,Comment
    '''
    names_table: pl.DataFrame
    supported_residues: List[NerpaResidue]
    norine_monomer_cnts: Dict[NorineMonomerName, int]

    _cache: Dict[Tuple[str, str], NRP_Monomer] = None
    default_frequencies: MonomersDefaultFrequencies = None

    mon_to_int: Dict[NRP_Monomer, MonCode] = None
    int_to_mon: Dict[MonCode, NRP_Monomer] = None


    def __post_init__(self):
        self._cache = {}

        # Validation
        assert all(res in self.names_table['NerpaResidue']
                   for res in self.supported_residues), \
            f"Some supported residues are absent in the names table: " \
            f"{set(self.supported_residues) - set(self.names_table['NerpaResidue'])}"
        residues_from_norine = {self.parsed_name(norine_name, 'rBAN/Norine').residue
                                for norine_name in self.norine_monomer_cnts}
        assert set(self.supported_residues) - residues_from_norine == set(), \
            f"Supported residues {set(self.supported_residues) - residues_from_norine} are missing in Norine"

        self._set_default_frequencies(missing_residue_freq=0.01)  # TODO: move to config

        self.int_to_mon = {}
        self.mon_to_int = {}
        for mon_res_int, mon_res in enumerate(self.supported_residues + SPECIAL_RESIDUES):
            for meth_int, methylated in enumerate([False, True]):
                for chir_int, chirality in enumerate([Chirality.L, Chirality.D, Chirality.UNKNOWN]):
                    for is_pks_hybrid_int, is_pks_hybrid in enumerate([False, True]):
                        mon = NRP_Monomer(residue=NerpaResidue(mon_res),
                                          methylated=methylated,
                                          chirality=chirality,
                                          is_pks_hybrid=is_pks_hybrid)
                        mon_int = mon_res_int * 12 + meth_int * 6 + chir_int * 2 + is_pks_hybrid_int
                        self.mon_to_int[mon] = MonCode(mon_int)
                        self.int_to_mon[MonCode(mon_int)] = mon

    def is_proper_monomer(self, mon: NRP_Monomer) -> bool:
        return all([not mon.is_pks_hybrid,
                    mon.chirality != Chirality.UNKNOWN,
                    mon.residue not in [PKS_RESIDUE, NOT_NRPS_RESIDUE]])

    def proper_monomers(self) -> Iterable[NRP_Monomer]:
        """
        Get a list of all proper monomers (excluding UNKNOWN and PKS).
        """
        return filter(self.is_proper_monomer, self.mon_to_int.keys())

    def monomer_default_freq(self,
                             mon: NRP_Monomer,
                             handle_pks_hybrids: bool = False) -> float:  # actually Prob (unannotated to avoid circular import)
        return math.e ** self.monomer_default_log_freq(mon, handle_pks_hybrids)

    def monomer_default_log_freq(self,
                                 mon: NRP_Monomer,
                                 handle_pks_hybrids: bool = False) -> float:  # actually LogProb (unannotated to avoid circular import)
        res = mon.residue if not mon.is_pks_hybrid or handle_pks_hybrids else UNKNOWN_RESIDUE
        residue_log_freq = log(self.default_frequencies.residue[res])
        methylation_log_freq = log(self.default_frequencies.methylation) \
            if mon.methylated else log(1 - self.default_frequencies.methylation)

        d_chr_freq = self.default_frequencies.d_chirality
        match mon.chirality:
            case Chirality.D:
                chirality_log_freq = log(d_chr_freq)
            case Chirality.L:
                chirality_log_freq = log(1 - d_chr_freq)
            case Chirality.UNKNOWN:
                chirality_log_freq = (d_chr_freq * log(d_chr_freq)
                                    + (1 - d_chr_freq) * log(1 - d_chr_freq))

        return residue_log_freq + methylation_log_freq + chirality_log_freq

    def _set_default_frequencies(self, missing_residue_freq: Optional[Prob] = None):
        """
        Set the default frequencies for residues, methylation, and D-chirality.

        Args:
            missing_residue_freq: Optional default frequency for residues that are in
                                 supported_residues but not found in Norine data.
                                 If None, such residues will have frequency 0.
        """
        parsed_monomers_cnts: Dict[NRP_Monomer, int] = {}
        for monomer_name, cnt in self.norine_monomer_cnts.items():
            parsed_monomer = self.parsed_name(monomer_name, name_format='rBAN/Norine')
            parsed_monomers_cnts[parsed_monomer] = parsed_monomers_cnts.get(parsed_monomer, 0) + cnt

        total_monomers = sum(parsed_monomers_cnts.values())
        total_methylated = sum(cnt for monomer, cnt in parsed_monomers_cnts.items()
                               if monomer.methylated)
        total_d_chirality = sum(cnt for monomer, cnt in parsed_monomers_cnts.items()
                                if monomer.chirality == Chirality.D)
        norine_methylation_freq = Prob(total_methylated / total_monomers)
        norine_d_chirality_freq = Prob(total_d_chirality / total_monomers)

        nerpa_residue_cnts: Dict[NerpaResidue, int] = {}
        for mon, cnt in parsed_monomers_cnts.items():
            nerpa_residue_cnts[mon.residue] = nerpa_residue_cnts.get(mon.residue, 0) + cnt

        residue_freqs = {residue: Prob(cnt / total_monomers)
                         for residue, cnt in nerpa_residue_cnts.items()}

        # Handle missing residues if a default frequency is provided
        if missing_residue_freq is not None:
            missing_residues = [res for res in self.supported_residues + SPECIAL_RESIDUES
                                if res not in residue_freqs]

            # Calculate total frequency to redistribute
            total_missing_freq = missing_residue_freq * len(missing_residues)

            scale_factor = (1.0 - total_missing_freq) / sum(residue_freqs.values())
            residue_freqs = {res: Prob(freq * scale_factor)
                             for res, freq in residue_freqs.items()}

            for res in missing_residues:
                residue_freqs[res] = missing_residue_freq

        assert all(res in residue_freqs
                   for res in self.supported_residues + SPECIAL_RESIDUES), \
            f'Some supported residues are missing in the default frequencies: ' \
            f'{set(self.supported_residues + SPECIAL_RESIDUES) - set(residue_freqs.keys())}'

        assert math.isclose(sum(residue_freqs.values()), 1.0), \
            f'Sum of residue frequencies is {sum(residue_freqs.values())}, expected 1.0'

        #print(f'Default frequencies set: residues \n',
        #      '\n'.join(f'{res}: {freq}' for res, freq in dict(residue_freqs).items()))

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
                    name_format: Literal['antiSMASH_short', 'rBAN/Norine', 'PARAS'],
                    log: Optional[NerpaLogger] = None) -> NRP_Monomer:
        if name_format not in ('antiSMASH_short', 'rBAN/Norine', 'PARAS'):
            raise ValueError(f'Unsupported name format: {name_format}')
        if (name, name_format) in self._cache:
            return self._cache[(name, name_format)]

        if name_format == 'rBAN/Norine' and name.startswith('X'):
            return UNKNOWN_MONOMER

        rows = (self.names_table
                .filter((pl.col('MonomerName') == name)
                         & (pl.col('NameFormat') == name_format)))

        if rows.is_empty():  # monomer name not found
            '''
            if log is not None:
                log.error(f'Name {name} not found in the names table. Parsing as UNKNOWN_MONOMER')
            else:
                #with open('missing_norine_names.txt', 'a') as f:
                #    f.write(f'{name}\n')
                print(f'WARNING: Name {name} not found in the names table. Parsing as UNKNOWN_MONOMER')
            '''
            return UNKNOWN_MONOMER

        row = rows.to_dicts()[0]  # it should be exactly one row

        mon_type = row['Type']
        if mon_type == 'PKS':
            return PKS_MONOMER

        if mon_type == 'NOT_NRPS':
            return NOT_NRPS_MONOMER

        assert mon_type in ('NRPS', 'NA', 'PKS_Hybrid'), \
            f'Unexpected monomer type: {name} -- {mon_type}'

        residue = NerpaResidue(row['NerpaResidue'])
        if residue not in self.supported_residues:
            residue = UNKNOWN_RESIDUE
        methylated = 'MT' in eval(row['Modifications'])
        chirality = Chirality.D if 'E' in eval(row['Modifications']) else Chirality.L   # TODO: make a separate column for chirality

        self._cache[(name, name_format)] = NRP_Monomer(residue=residue,
                                                       methylated=methylated,
                                                       chirality=chirality,
                                                       is_pks_hybrid=mon_type == 'PKS_Hybrid')
        return self._cache[(name, name_format)]
