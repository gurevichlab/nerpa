from __future__ import annotations
from typing import (
    Any,
    Dict,
    FrozenSet,
    List,
    Optional,
    Tuple
)
from src.monomer_names_helper import (
    MonomerResidue,
    NorineMonomerName,
    NRP_Monomer_Modification,
)
from src.antismash_parsing.antismash_parser_types import GeneId
from src.antismash_parsing.location_features import ModuleLocFeatures, ModuleLocFeature
from src.rban_parsing.monomer_features_types import MonomerFeature, MonomerFeatures
from enum import auto, Enum
import yaml
import os
from dataclasses import asdict, dataclass

Logger = Any  # some magical logger used throughout the pipeline

SMILES = str
LogProb = float
ResidueScores = Dict[MonomerResidue, LogProb]


def enum_representer(dumper, e: Enum):
    return dumper.represent_scalar(f'!{e.__class__.__name__}', e.name)


class Chirality(Enum):
    L = auto()
    D = auto()
    UNKNOWN = auto()

yaml.add_representer(Chirality, enum_representer)
yaml.add_representer(NRP_Monomer_Modification, enum_representer)
yaml.add_representer(ModuleLocFeature, enum_representer)
yaml.add_representer(MonomerFeature, enum_representer)


@dataclass
class NRP_Monomer:
    residue: MonomerResidue
    modifications: Tuple[NRP_Monomer_Modification, ...]
    chirality: Chirality
    rban_name: NorineMonomerName
    rban_idx: int
    monomer_features: MonomerFeatures

    @classmethod
    def from_yaml_dict(cls, data: dict) -> NRP_Monomer:
        return cls(residue=MonomerResidue(data['residue']),
                   modifications=tuple(NRP_Monomer_Modification[mod]
                                       for mod in data['modifications']),
                   chirality=Chirality[data['chirality']],
                   rban_name=NorineMonomerName(data['rban_name']),
                   rban_idx=data['rban_idx'],
                   monomer_features=tuple(MonomerFeature[loc_feature]
                                          for loc_feature in data['monomer_features']))


NRP_MONOMER_DUMMY = NRP_Monomer(residue=MonomerResidue(''), modifications=(), chirality=Chirality.UNKNOWN,
                                rban_name=NorineMonomerName(''), rban_idx=-1, monomer_features=())


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
        return cls(gene_id=data['gene_id'],
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


BGC_MODULE_DUMMY = BGC_Module(gene_id=GeneId(''), fragment_idx=-1, module_loc=(), a_domain_idx=-1,
                              residue_score={}, modifications=(), iterative_module=False, iterative_gene=False)


BGC_Fragment = List[BGC_Module]


@dataclass
class BGC_Variant:
    variant_idx: int
    genome_id: str
    bgc_idx: int
    modules: List[BGC_Module]

    @classmethod
    def from_yaml_dict(cls, data: dict) -> BGC_Variant:
        return cls(variant_idx=data['variant_idx'],
                   genome_id=data['genome_id'],
                   bgc_idx=data['bgc_idx'],
                   modules=[BGC_Module.from_yaml_dict(module_dict)
                               for module_dict in data['modules']])

    def has_pks_domains(self) -> bool:
        return any('PKS' in feature.name
                   for module in self.modules
                   for feature in module.module_loc)


@dataclass
class NRP_Fragment:
    monomers: List[NRP_Monomer]
    is_cyclic: bool

    @classmethod
    def from_yaml_dict(cls, data: dict) -> NRP_Fragment:
        return cls(is_cyclic=data['is_cyclic'],
                   monomers=list(map(NRP_Monomer.from_yaml_dict, data['monomers'])))


@dataclass
class NRP_Variant:
    variant_idx: int
    nrp_id: str
    fragments: List[NRP_Fragment]
    isolated_unknown_monomers: List[NRP_Monomer] = None

    @classmethod
    def from_yaml_dict(cls, data: dict) -> NRP_Variant:
        return cls(variant_idx=data['variant_idx'],
                   nrp_id=data['nrp_id'],
                   fragments=list(map(NRP_Fragment.from_yaml_dict, data['fragments'])),
                   isolated_unknown_monomers=list(map(NRP_Monomer.from_yaml_dict, data['isolated_monomers']))
                   if data.get('isolated_monomers') else None)
