from __future__ import annotations
from enum import Enum, auto

import math
import yaml
from typing import (
    Any,
    Dict,
    List,
    NamedTuple,
    Optional,
    Tuple,
)
from dataclasses import asdict, dataclass

from src.antismash_parsing.antismash_parser_types import GeneId, BGC_ID, antiSMASH_metadata
from src.antismash_parsing.genomic_context import ModuleGenomicContext, ModuleGenomicContextFeature
from src.build_output.yaml_representers import enum_representer
from src.general_type_aliases import LogProb
from src.monomer_names_helper import MonomerResidue

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

    def __eq__(self, other: BGC_Module) -> bool:
        features_coincide = all([
            self.genomic_context == other.genomic_context,
            self.modifications == other.modifications,
            self.iterative_module == other.iterative_module,
            self.iterative_gene == other.iterative_gene])

        residue_scores_coincide = all([
            math.isclose(score, other.residue_score[res],
                         abs_tol=1e-3)
            for res, score in self.residue_score.items()])

        return features_coincide and residue_scores_coincide


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
        data['metadata'] = self.metadata.to_dict()
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

    def __eq__(self, other: BGC_Variant) -> bool:
        if not isinstance(other, BGC_Variant):
            return NotImplemented
        return self.modules == other.modules


class BGC_Variants_Info(NamedTuple):
    bgc_variants: List[BGC_Variant]
    bgc_id_to_repr_id: Dict[BGC_Variant_ID, BGC_Variant_ID]

    def get_representative_bgc_variants(self) -> List[BGC_Variant]:
        repr_ids = set(self.bgc_id_to_repr_id.values())
        return [bgc for bgc in self.bgc_variants if bgc.bgc_variant_id in repr_ids]

