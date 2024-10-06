from __future__ import annotations
from typing import (
    Any,
    Dict,
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


@dataclass
class NRP_Monomer:
    residue: MonomerResidue
    modifications: Tuple[NRP_Monomer_Modification, ...]
    chirality: Chirality
    rban_name: NorineMonomerName
    rban_idx: int
    is_hybrid: bool = False

    @classmethod
    def from_yaml_dict(cls, data: dict) -> NRP_Monomer:
        return cls(residue=MonomerResidue(data['residue']),
                   modifications=tuple(NRP_Monomer_Modification[mod]
                                       for mod in data['modifications']),
                   chirality=Chirality[data['chirality']],
                   rban_name=NorineMonomerName(data['rban_name']),
                   rban_idx=data['rban_idx'],
                   is_hybrid=data.get('hybrid', False))


class BGC_Module_Modification(Enum):
    EPIMERIZATION = auto()
    METHYLATION = auto()

yaml.add_representer(BGC_Module_Modification, enum_representer)


@dataclass
class ModuleLocation:
    start_of_gene: bool
    end_of_gene: bool
    start_of_fragment: bool
    end_of_fragment: bool
    pks_upstream: bool
    pks_downstream: bool


@dataclass
class BGC_Module:
    gene_id: GeneId
    module_loc: ModuleLocation
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
                   a_domain_idx=data['a_domain_idx'],
                   residue_score=data['residue_score'],
                   modifications=tuple(BGC_Module_Modification[mod]
                                       for mod in data['modifications']),
                   iterative_module=data['iterative_module'],
                   iterative_gene=data['iterative_gene'],
                   aa10_code=data.get('aa10_code') if data.get('aa10_code') != '---'  else None,
                   aa34_code=data.get('aa34_code') if data.get('aa34_code') != '---' else None)


BGC_Fragment = List[BGC_Module]


@dataclass
class BGC_Variant:
    variant_idx: int
    genome_id: str
    bgc_idx: int
    fragments: List[BGC_Fragment]
    has_pks_domains: bool = False

    @classmethod
    def from_yaml_dict(cls, data: dict) -> BGC_Variant:
        return cls(variant_idx=data['variant_idx'],
                   genome_id=data['genome_id'],
                   bgc_idx=data['bgc_idx'],
                   fragments=[[BGC_Module.from_yaml_dict(module_dict)
                               for module_dict in fragment_list]
                               for fragment_list in data['fragments']],
                   has_pks_domains=data.get('has_pks_domains', False))


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

    @classmethod
    def from_yaml_dict(cls, data: dict) -> NRP_Variant:
        return cls(variant_idx=data['variant_idx'],
                   nrp_id=data['nrp_id'],
                   fragments=list(map(NRP_Fragment.from_yaml_dict, data['fragments'])))
