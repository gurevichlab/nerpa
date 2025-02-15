from __future__ import annotations
from typing import (
    Any,
    Dict,
    FrozenSet,
    List,
    NamedTuple,
    Optional,
    Tuple
)
from src.monomer_names_helper import (
    MonomerResidue,
    NorineMonomerName,
    NRP_Monomer_Modification,
    Chirality,
    NRP_Monomer
)
from src.antismash_parsing.antismash_parser_types import GeneId
from src.antismash_parsing.location_features import ModuleLocFeatures, ModuleLocFeature
from src.rban_parsing.rban_monomer import rBAN_Monomer
from src.monomer_names_helper import enum_representer
from enum import auto, Enum
from dataclasses import asdict, dataclass
from itertools import chain, permutations, product
import yaml

Logger = Any  # some magical logger used throughout the pipeline

SMILES = str
LogProb = float
ResidueScores = Dict[MonomerResidue, LogProb]


class BGC_Module_Modification(Enum):
    EPIMERIZATION = auto()
    METHYLATION = auto()

yaml.add_representer(BGC_Module_Modification, enum_representer)


@dataclass
class BGC_Module:
    gene_id: GeneId
    fragment_idx: int
    module_loc: ModuleLocFeatures
    a_domain_idx: int  # not the same as module_idx because modules with no a_domain are skipped
    residue_score: ResidueScores
    modifications: Tuple[BGC_Module_Modification, ...]
    iterative_module: bool
    iterative_gene: bool
    aa10_code: Optional[str] = None  # for debug only
    aa34_code: Optional[str] = None  # for debug only

    @classmethod
    def from_yaml_dict(cls, data: dict) -> BGC_Module:  # TODO: use dacite from_dict or smth
        try:
            res = cls(gene_id=data['gene_id'],
                       fragment_idx=data['fragment_idx'],
                       a_domain_idx=data['a_domain_idx'],
                       module_loc=tuple(ModuleLocFeature[loc_feature]
                                        for loc_feature in data['module_loc']),
                       residue_score=data['residue_score'],
                       modifications=tuple(BGC_Module_Modification[mod]
                                           for mod in data['modifications']),
                       iterative_module=data['iterative_module'],
                       iterative_gene=data['iterative_gene'],
                       aa10_code=data.get('aa10_code') if data.get('aa10_code') != '---'  else None,
                       aa34_code=data.get('aa34_code') if data.get('aa34_code') != '---' else None)
        except:
            raise
        return res


BGC_MODULE_DUMMY = BGC_Module(gene_id=GeneId(''), fragment_idx=-1, module_loc=(), a_domain_idx=-1,
                              residue_score={}, modifications=(), iterative_module=False, iterative_gene=False)


BGC_Fragment = List[BGC_Module]


@dataclass
class BGC_Variant:
    variant_idx: int
    genome_id: str
    contig_idx: int
    bgc_idx: int
    modules: List[BGC_Module]

    @classmethod
    def from_yaml_dict(cls, data: dict) -> BGC_Variant:
        return cls(variant_idx=data['variant_idx'],
                   genome_id=data['genome_id'],
                   bgc_idx=data['bgc_idx'],
                   contig_idx=data.get('contig_idx', 0),  # for compatibility with old yaml files
                   modules=[BGC_Module.from_yaml_dict(module_dict)
                               for module_dict in data['modules']])


    def has_pks_domains(self) -> bool:
        return any('PKS' in feature.name
                   for module in self.modules
                   for feature in module.module_loc)


@dataclass
class NRP_Fragment:
    monomers: List[rBAN_Monomer]
    is_cyclic: bool

    @classmethod
    def from_yaml_dict(cls, data: dict) -> NRP_Fragment:
        return cls(is_cyclic=data['is_cyclic'],
                   monomers=list(map(rBAN_Monomer.from_yaml_dict, data['monomers'])))


@dataclass
class NRP_Variant:
    variant_idx: int
    nrp_id: str
    fragments: List[NRP_Fragment]
    isolated_unknown_monomers: List[rBAN_Monomer] = None

    @classmethod
    def from_yaml_dict(cls, data: dict) -> NRP_Variant:
        return cls(variant_idx=data['variant_idx'],
                   nrp_id=data['nrp_id'],
                   fragments=list(map(NRP_Fragment.from_yaml_dict, data['fragments'])),
                   isolated_unknown_monomers=list(map(rBAN_Monomer.from_yaml_dict, data['isolated_monomers']))
                   if data.get('isolated_monomers') else None)

    def get_linearizations(self, max_num_fragments_to_permute: int = 3) -> List[List[rBAN_Monomer]]:
        def fragment_linearizations(fragment: NRP_Fragment) -> List[List[rBAN_Monomer]]:
            if fragment.is_cyclic:
                return [fragment.monomers[i:] + fragment.monomers[:i]
                        for i in range(len(fragment.monomers))]
            else:
                return [fragment.monomers]

        fragments_permutations = permutations(self.fragments) \
            if len(self.fragments) <= max_num_fragments_to_permute else [self.fragments]

        linearizations = []
        for fragments in fragments_permutations:
            for fragments_linearization in product(*(fragment_linearizations(fragment)
                                                     for fragment in fragments)):
                linearizations.append(list(chain(*fragments_linearization)))
        return linearizations


class Match_BGC_Variant_Info(NamedTuple):
    genome_id: str
    contig_idx: int
    bgc_idx: int
    variant_idx: int

    @classmethod
    def from_bgc_variant(cls, bgc_variant: BGC_Variant) -> Match_BGC_Variant_Info:
        return cls(genome_id=bgc_variant.genome_id,
                   contig_idx=bgc_variant.contig_idx,
                   bgc_idx=bgc_variant.bgc_idx,
                   variant_idx=bgc_variant.variant_idx)

    def get_antismash_id(self) -> str:
        return f'r{self.contig_idx}c{self.bgc_idx}'


class Match_NRP_Variant_Info(NamedTuple):
    nrp_id: str
    variant_idx: int

    @classmethod
    def from_nrp_variant(cls, nrp_variant: NRP_Variant) -> Match_NRP_Variant_Info:
        return cls(nrp_id=nrp_variant.nrp_id,
                   variant_idx=nrp_variant.variant_idx)
