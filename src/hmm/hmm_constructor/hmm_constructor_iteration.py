from __future__ import annotations

from src.hmm.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    DetailedHMMStateType,
    DetailedHMMState,
    StateIdx,
)

ET = DetailedHMMEdgeType
ST = DetailedHMMStateType


def add_iteration_subgraph(constructor,  # type: HMM_Constructor,
                           module_idx: int,
                           node_after_subgraph: StateIdx) -> StateIdx:
    '''
    Builds a subgraph that chooses whether to iterate the module and/or gene
    and then transitions to node_after_subgraph.
    Returns the index of the root of the subgraph.
    '''
    module = constructor.bgc_variant.modules[module_idx]
    iterative_module = module.iterative_module  # for symmetry with iterative_gene
    iterative_gene = (
        module.iterative_gene and
        constructor.gene_borders[module.gene_id][1] == module_idx
    )  # should be the last module in the gene to allow gene iteration

    match (iterative_module, iterative_gene):
        case (False, False):  # no iteration
            return node_after_subgraph
        case (True, False) | (False, True):  # only one type of iteration
            if iterative_module:
                iteration_return_point = constructor.module_idx_to_subgraph_root[module_idx]
                iteration_state_type = ST.CHOOSE_IF_ITERATE_MODULE
                iteration_edge_type = ET.ITERATE_MODULE
            else:
                fst_module_in_gene_idx = constructor.gene_borders[module.gene_id][0]
                iteration_return_point = constructor.module_idx_to_subgraph_root[fst_module_in_gene_idx]
                iteration_state_type = ST.CHOOSE_IF_ITERATE_GENE
                iteration_edge_type = ET.ITERATE_GENE

            constructor.states.append(DetailedHMMState(state_type=iteration_state_type,
                                                       related_module_idx=module_idx))
            choose_if_iterate_root = len(constructor.states) - 1

            constructor.add_edge(choose_if_iterate_root, iteration_return_point, iteration_edge_type)
            constructor.add_edge(choose_if_iterate_root, node_after_subgraph, ET.NO_ITERATION)
            return choose_if_iterate_root
        case (True, True):  # both module and gene can be iterated
            constructor.states.extend([
                DetailedHMMState(state_type=ST.CHOOSE_IF_ITERATE_MODULE,
                                 related_module_idx=module_idx),
                DetailedHMMState(state_type=ST.CHOOSE_IF_ITERATE_GENE,
                                 related_module_idx=module_idx),
            ])

            choose_if_iterate_module = len(constructor.states) - 2
            choose_if_iterate_gene = len(constructor.states) - 1
            module_iteration_return_point = constructor.module_idx_to_subgraph_root[module_idx]
            fst_module_in_gene_idx = constructor.gene_borders[module.gene_id][0]
            gene_iteration_return_point = constructor.module_idx_to_subgraph_root[fst_module_in_gene_idx]

            constructor.add_edge(choose_if_iterate_module, module_iteration_return_point, ET.ITERATE_MODULE)
            constructor.add_edge(choose_if_iterate_module, choose_if_iterate_gene, DetailedHMMEdgeType.NO_ITERATION,)
            constructor.add_edge(choose_if_iterate_gene, gene_iteration_return_point, ET.ITERATE_GENE)
            constructor.add_edge(choose_if_iterate_gene, node_after_subgraph, DetailedHMMEdgeType.NO_ITERATION,)

            return choose_if_iterate_module
