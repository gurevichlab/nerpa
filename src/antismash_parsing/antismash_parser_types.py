from __future__ import annotations

from pathlib import Path
from typing import (
    Dict,
    List,
    Union,
    NamedTuple,
    NewType,
    Optional,
)
from dataclasses import dataclass
from enum import Enum, auto
from src.monomer_names_helper import (
    antiSMASH_MonomerName,
    AA10,
    AA34
)


antiSMASH_record = NewType('antiSMASH_record', dict)
GeneId = NewType('GeneId', str)

class STRAND(Enum):
    FORWARD = auto()
    REVERSE = auto()


class Coords(NamedTuple):
    start: int
    end: int
    strand: STRAND

    # TODO: I am not sure if this is correct
    @classmethod
    def from_hmm_hit(cls, query_start: int, query_end: int) -> Coords:
        if query_start < query_end:
            return cls(query_start, query_end, STRAND.FORWARD)
        else:
            return cls(query_end, query_start, STRAND.REVERSE)

    def is_in(self, other: Coords) -> bool:
        return self.start >= other.start and self.end <= other.end

    def left_from(self, other: Coords) -> bool:
        return self.end < other.start

    def right_from(self, other: Coords) -> bool:
        return self.start > other.end


class SVM_LEVEL(Enum):
    SINGLE_AMINO = auto()
    SMALL_CLUSTER = auto()
    LARGE_CLUSTER = auto()
    PHYSIOCHEMICAL_CLASS = auto()


class SVM_Prediction(NamedTuple):
    score: float
    substrates: List[antiSMASH_MonomerName]


SVM_Predictions = Dict[SVM_LEVEL, SVM_Prediction]


@dataclass
class A_Domain:
    aa10: AA10
    aa34: AA34
    svm: SVM_Predictions


class DomainType(Enum):
    A = auto()
    PKS = auto()

    PCP = auto()

    C = auto()
    C_STARTER = auto()
    C_LCL = auto()
    C_DCL = auto()
    C_DUAL = auto()


    E = auto()
    MT = auto()

    TE_TD = auto()

    CTERM = auto()
    NTERM = auto()

    def in_c_domain_group(self) -> bool:
        return self in {DomainType.C, DomainType.C_STARTER, DomainType.C_LCL, DomainType.C_DCL, DomainType.C_DUAL}


class BGC_Module_ID(NamedTuple):
    gene_id: GeneId
    module_idx: int


@dataclass
class Module:
    a_domain: Optional[A_Domain]
    domains_sequence: List[DomainType]


@dataclass
class Gene:
    gene_id: GeneId
    coords: Coords
    modules: List[Module]  # modules are in the order of appearance in the gene
    orphan_c_at_start: bool = False
    orphan_c_at_end: bool = False


class BGC_ID(NamedTuple):
    antiSMASH_file: Path
    contig_idx: int
    bgc_idx: int

    def to_dict(self):
        return {'antiSMASH_file': str(self.antiSMASH_file),
                'contig_idx': self.contig_idx,
                'bgc_idx': self.bgc_idx}

    @classmethod
    def from_dict(cls, data: dict) -> BGC_ID:
        return cls(
            antiSMASH_file=Path(data['antiSMASH_file']),
            contig_idx=data['contig_idx'],
            bgc_idx=data['bgc_idx']
        )

    def __str__(self):
        return f'{self.antiSMASH_file}_{self.contig_idx}_{self.bgc_idx}'

    def to_str_short(self):
        return f'{self.antiSMASH_file.stem}_{self.contig_idx}_{self.bgc_idx}'


@dataclass
class antiSMASH_metadata:
    antismash_json: Path
    contig_idx: int
    sequence_file: Optional[str]
    id: Optional[str]
    organism: Optional[str]
    taxonomy: Optional[List[str]]  # e.g. ['Bacteria', 'Proteobacteria', 'Alphaproteobacteria', 'Hyphomicrobiales', 'Rhizobiaceae', 'Rhizobium/Agrobacterium group', 'Agrobacterium', 'Agrobacterium tumefaciens complex']

    @classmethod
    def from_record(cls,
                    record: antiSMASH_record,
                    ctg_idx: int,
                    antismash_json_file: Path) -> antiSMASH_metadata:  # contig_data is an element of antismash_json['records']
        sequence_file = record.get('input_file')
        try:
            contig_data = record['contigs'][ctg_idx]
            annotations = contig_data['annotations']
        except KeyError:
            return cls(antismash_json=antismash_json_file,
                       contig_idx=ctg_idx,
                       sequence_file=sequence_file,
                       id=None, organism=None, taxonomy=None)

        id = annotations.get('id')
        organism = annotations.get('organism')
        taxonomy = annotations.get('taxonomy')
        return cls(antismash_json=antismash_json_file,
                   contig_idx=ctg_idx,
                   sequence_file=sequence_file,
                   id=id,
                   organism=organism,
                   taxonomy=taxonomy)


@dataclass
class BGC_Cluster:
    bgc_id: BGC_ID
    genes: List[Gene]
    metadata: antiSMASH_metadata

    def has_pks_domains(self) -> bool:
        return any(DomainType.PKS in module.domains_sequence
                   for gene in self.genes
                   for module in gene.modules)

    def has_a_domains(self) -> bool:
        return any(module.a_domain is not None
                   for gene in self.genes
                   for module in gene.modules)


@dataclass
class Fragmented_BGC_Cluster:
    bgc_id: BGC_ID
    fragments: List[List[Gene]]

