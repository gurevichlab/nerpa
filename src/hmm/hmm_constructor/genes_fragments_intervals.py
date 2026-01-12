from typing import List, Dict, Tuple
from src.antismash_parsing.bgc_variant_types import (
    BGC_ID,
    BGC_Variant_ID,
    BGC_Variant,
    BGC_Module,
    GeneId,
)
from src.rban_parsing.nrp_variant_types import NRP_Variant, NRP_Variant_ID

from src.general_type_aliases import LogProb
from collections import defaultdict


def get_genes_intervals(assembly_line: List[BGC_Module]) -> Dict[GeneId, Tuple[int, int]]:  # gene_id -> (start, end)
    return defaultdict(lambda: (-1, -1),
                       {gene_id: (next(i
                                       for i, module in enumerate(assembly_line)
                                       if module.gene_id == gene_id),
                                  next(len(assembly_line) - 1 - i
                                       for i, module in enumerate(reversed(assembly_line))
                                       if module.gene_id == gene_id))
                        for gene_id in {module.gene_id for module in assembly_line}})


def get_fragments_intervals(assembly_line: List[BGC_Module]) -> Dict[int, Tuple[int, int]]:  # fragment_idx -> (start, end)
    return defaultdict(lambda: (-1, -1),
                       {fragment_idx: (next(i
                                            for i, module in enumerate(assembly_line)
                                            if module.fragment_idx == fragment_idx),
                                       next(len(assembly_line) - 1 - i
                                            for i, module in enumerate(reversed(assembly_line))
                                            if module.fragment_idx == fragment_idx))
                        for fragment_idx in {module.fragment_idx for module in assembly_line}})
