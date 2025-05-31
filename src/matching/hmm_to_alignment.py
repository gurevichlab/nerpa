from __future__ import annotations
from typing import TYPE_CHECKING

from src.testing.check_matches import bgc_modules_coincide

if TYPE_CHECKING:
    from src.matching.detailed_hmm import DetailedHMM
from typing import Dict, List, Tuple
from src.matching.genes_fragments_intervals import get_genes_intervals, get_fragments_intervals
from src.matching.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    DetailedHMMStateType
)
from itertools import pairwise
from src.matching.alignment_type import Alignment
from src.matching.alignment_step_type import (
    AlignmentStep,
    AlignmentStep_BGC_Module_Info,
)
from src.rban_parsing.rban_monomer import rBAN_Monomer
from pathlib import Path
from itertools import count


def hmm_path_to_alignment(hmm: DetailedHMM,
                          path: List[int],
                          nrp_monomers: List[rBAN_Monomer]) -> Alignment:
    #free_idx = next(i for i in count()
    #                if not Path(f'{hmm.bgc_variant.bgc_variant_id.bgc_id.genome_id}_path_{i}.png').exists())
    #hmm.draw(Path(f'{hmm.bgc_variant.bgc_variant_id.bgc_id.genome_id}_path_{free_idx}.png'), path)  # for debugging
    #hmm.draw(Path(f'hmm_path.png'), path)  # for debugging
    #print('bgc: ', hmm.bgc_variant.bgc_variant_id.bgc_id.genome_id)
    #print('path: ', path)
    #print('nrp monomers: ', [mon.residue for mon in nrp_monomers])
    ET = DetailedHMMEdgeType  # for convenience

    genes_intervals = get_genes_intervals(hmm.bgc_variant.modules)
    alignment = []
    mon_idx = 0  # next nrp_monomer to be matched
    module_idx = 0  # next bgc_module to be matched
    for i, (edge_from, edge_to) in enumerate(pairwise(path)):
        edge = hmm.transitions[edge_from][edge_to]
        bgc_module = hmm.bgc_variant.modules[module_idx] \
            if module_idx < len(hmm.bgc_variant.modules) else None
        bgc_module_info = AlignmentStep_BGC_Module_Info.from_bgc_module(bgc_module) \
            if bgc_module is not None else None
        nrp_mon = nrp_monomers[mon_idx] \
            if mon_idx < len(nrp_monomers) else None

        match edge.edge_type:
            case (ET.SKIP_MODULE_AT_START | ET.SKIP_MODULE | ET.SKIP_MODULE_AT_END | ET.SKIP_MODULE_END_MATCHING):
                alignment.append(AlignmentStep(
                    bgc_module=bgc_module_info,
                    nrp_monomer=None,
                    score=edge.weight,
                    match_detailed_score=None,
                    step_type=edge.edge_type))
                module_idx += 1

            case ET.SKIP_GENE:
                    gene_id = bgc_module.gene_id
                    next_module_idx = genes_intervals[gene_id][1] + 1
                    for skipped_module_idx in range(module_idx, next_module_idx):
                        skipped_module = hmm.bgc_variant.modules[skipped_module_idx]
                        alignment.append(AlignmentStep(
                            bgc_module=AlignmentStep_BGC_Module_Info.from_bgc_module(skipped_module),
                            nrp_monomer=None,
                            score=edge.weight,
                            match_detailed_score=None,
                            step_type=edge.edge_type))
                    module_idx = next_module_idx

            case (ET.START_INSERTING_AT_START
                  | ET.INSERT_AT_START
                  | ET.START_INSERTING
                  | ET.INSERT
                  | ET.START_INSERTING_AT_END
                  | ET.INSERT_AT_END):
                alignment.append(AlignmentStep(
                    bgc_module=None,
                    nrp_monomer=nrp_mon,
                    score=edge.weight,
                    match_detailed_score=None,
                    step_type=edge.edge_type))
                mon_idx += 1

            case ET.MATCH:
                match_detailed_score = hmm.hmm_helper.normalized_match_detailed_score(bgc_module,
                                                                                      nrp_mon.to_base_mon(),
                                                                                      hmm.bgc_variant.has_pks_domains())
                emission_score = sum(match_detailed_score)
                alignment.append(AlignmentStep(bgc_module=bgc_module_info,
                                               nrp_monomer=nrp_mon,
                                               score=emission_score + edge.weight,
                                               match_detailed_score=match_detailed_score,
                                               step_type=edge.edge_type))
                module_idx += 1
                mon_idx += 1

            case ET.ITERATE_MODULE:
                # iterate module
                alignment.append(AlignmentStep(bgc_module=None,
                                               nrp_monomer=None,
                                               score=edge.weight,
                                               match_detailed_score=None,
                                               step_type=edge.edge_type))
                module_idx -= 1
            case ET.ITERATE_GENE:
                prev_module = hmm.bgc_variant.modules[module_idx - 1]
                gene_start = genes_intervals[prev_module.gene_id][0]
                alignment.append(AlignmentStep(bgc_module=None,
                                               nrp_monomer=None,
                                               score=edge.weight,
                                               match_detailed_score=None,
                                               step_type=edge.edge_type))
                module_idx = gene_start

        # auxiliary edges
            case (ET.START_SKIP_MODULES_AT_START
                  | ET.START_MATCHING | ET.END_MATCHING
                  | ET.NO_INSERTIONS | ET.END_INSERTING
                  | ET.NO_ITERATION):
                alignment.append(AlignmentStep(bgc_module=None,
                                               nrp_monomer=None,
                                               score=edge.weight,
                                               match_detailed_score=None,
                                               step_type=edge.edge_type))

            case _:
                raise ValueError(f'Unexpected edge type: {edge.edge_type}')

    return alignment