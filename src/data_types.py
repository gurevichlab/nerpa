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
from src.antismash_parsing.antismash_parser_types import GeneId, BGC_ID, antiSMASH_metadata
from src.antismash_parsing.genomic_context import ModuleGenomicContext, ModuleGenomicContextFeature
from src.rban_parsing.rban_monomer import rBAN_Monomer
from src.monomer_names_helper import enum_representer, AA10, AA34
from enum import auto, Enum
from dataclasses import asdict, dataclass
from itertools import chain, permutations, product
import yaml


SMILES = str
Prob = float
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
    a_domain_idx: int  # not the same as module_idx because modules with no a_domain are skipped
    genomic_context: ModuleGenomicContext
    residue_score: ResidueScores
    modifications: Tuple[BGC_Module_Modification, ...]
    iterative_module: bool
    iterative_gene: bool
    aa10_code: Optional[str] = None  # for debug only
    aa34_code: Optional[str] = None  # for debug only

    @classmethod
    def from_dict(cls, data: dict) -> BGC_Module:  # TODO: use dacite from_dict or smth
        try:
            res = cls(gene_id=data['gene_id'],
                      fragment_idx=data['fragment_idx'],
                      a_domain_idx=data['a_domain_idx'],
                      genomic_context=tuple(ModuleGenomicContextFeature[context_feature]
                                            for context_feature in data['genomic_context']),
                      residue_score=data['residue_score'],
                      modifications=tuple(BGC_Module_Modification[mod]
                                           for mod in data['modifications']),
                      iterative_module=data['iterative_module'],
                      iterative_gene=data['iterative_gene'],
                      aa10_code=data.get('aa10_code') if data.get('aa10_code') != '---'  else None,
                      aa34_code=data.get('aa34_code') if data.get('aa34_code') != '---' else None)
        except:
            raise  # for debugging
        return res


class A_Domain_ID(NamedTuple):
    gene_id: GeneId
    a_domain_idx: int  # not the same as module_idx because modules with no a_domain are skipped

    @classmethod
    def from_module(cls, module: BGC_Module) -> A_Domain_ID:
        return cls(gene_id=module.gene_id,
                   a_domain_idx=module.a_domain_idx)


BGC_MODULE_DUMMY = BGC_Module(gene_id=GeneId(''), fragment_idx=-1, genomic_context=(), a_domain_idx=-1,
                              residue_score={}, modifications=(), iterative_module=False, iterative_gene=False)


BGC_Fragment = List[BGC_Module]



class BGC_Variant_ID(NamedTuple):
    bgc_id: BGC_ID
    variant_idx: int

    @classmethod
    def from_dict(cls, data: dict) -> BGC_Variant_ID:
        return cls(
            bgc_id=BGC_ID.from_dict(data["bgc_id"]),
            variant_idx=data["variant_idx"]
        )

    def to_dict(self) -> dict:
        return {'bgc_id': self.bgc_id.to_dict(),
                'variant_idx': self.variant_idx}



@dataclass
class BGC_Variant:
    bgc_variant_id: BGC_Variant_ID
    modules: List[BGC_Module]
    metadata: Optional[antiSMASH_metadata]

    @classmethod
    def from_yaml_dict(cls, data: dict) -> BGC_Variant:
        try:
            metadata = antiSMASH_metadata(**data['metadata'])
        except:
            metadata = None
        return cls(
            bgc_variant_id=BGC_Variant_ID.from_dict(data["bgc_variant_id"]),
            modules=[BGC_Module.from_dict(module) for module in data["modules"]],
            metadata=metadata
        )

    def to_dict(self) -> Dict[str, Any]:
        data = asdict(self)
        data['bgc_variant_id'] = self.bgc_variant_id.to_dict()  # Convert NamedTuple manually
        return data


    def has_pks_domains(self) -> bool:
        return any('PKS' in feature.name
                   for module in self.modules
                   for feature in module.genomic_context)

    def get_genome_id(self) -> str:
        if self.metadata:
           return self.metadata.sequence_file.rsplit('.', 1)[0]
        else:
            return self.bgc_variant_id.bgc_id.antiSMASH_file.stem


@dataclass
class NRP_Fragment:
    monomers: List[rBAN_Monomer]
    is_cyclic: bool

    @classmethod
    def from_yaml_dict(cls, data: dict) -> NRP_Fragment:
        return cls(is_cyclic=data['is_cyclic'],
                   monomers=[rBAN_Monomer.from_list(mon_data_lst)
                             for mon_data_lst in data['monomers']])


class NRP_Variant_ID(NamedTuple):
    nrp_id: str
    variant_idx: int


@dataclass
class NRP_metadata:
    name: Optional[str]
    smiles: Optional[SMILES]
    origin: Optional[str]
    inchikey: Optional[str]
    source: Optional[str]

    @classmethod
    def from_dict(cls, data: dict) -> NRP_metadata:
        fields = cls.__dataclass_fields__.keys()
        class_data = {field: None for field in fields}
        for key, value in data.items():
            if key.lower() in fields:
                class_data[key.lower()] = data[key]
        return cls(**class_data)


@dataclass
class NRP_Variant:
    nrp_variant_id: NRP_Variant_ID
    fragments: List[NRP_Fragment]
    isolated_unknown_monomers: List[rBAN_Monomer]
    metadata: NRP_metadata

    @classmethod
    def from_yaml_dict(cls, data: dict) -> NRP_Variant:
        if 'metadata' in data:
            metadata = NRP_metadata.from_dict(data['metadata'])
        else:
            metadata = NRP_metadata(name=None, smiles=None, origin=None, inchikey=None, source=None)
        return cls(nrp_variant_id=NRP_Variant_ID(*data['nrp_variant_id']),
                   fragments=list(map(NRP_Fragment.from_yaml_dict, data['fragments'])),
                   isolated_unknown_monomers=list(map(rBAN_Monomer.from_yaml_dict, data['isolated_monomers']))
                   if data.get('isolated_monomers') else [],
                   metadata=metadata)

