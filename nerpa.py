import sys
from typing import List, Optional, NamedTuple, Tuple

import networkx as nx

from src.pipeline.pipeline_helper import PipelineHelper
from src.pipeline.logger import NerpaLogger
from src.rban_parsing.rban_parser import Parsed_rBAN_Record
from src.monomer_names_helper import MonomerNamesHelper, UNKNOWN_RESIDUE
from src.data_types import BGC_Variant
from src.write_results import write_yaml
from pathlib import Path
import yaml
from collections import defaultdict


def node_match(n1, n2):
    name1, name2 = n1.get('name'), n2.get('name')
    return any([name1 == name2,  # same rban name
                name1.startswith('X') and name2.startswith('X'),  # both unknown
                ':' in name1 and ':'])  # both lipid tails


def edge_match(e1, e2):
    return e1.get('bond_type') == e2.get('bond_type')


class EditInfo(NamedTuple):
    edge_changes: int
    num_indels: int
    subs: List[Tuple[str, str]]

    def to_dict(self):
        return {'edge_changes': self.edge_changes,
                'num_indels': self.num_indels,
                'subs': self.subs}


def make_edits(G1: nx.DiGraph,
               G2: nx.DiGraph,
               edit_info: Tuple[list, list, float]) -> EditInfo:
    node_maps, edge_maps, cost = edit_info

    # 1: count num edge changes
    num_edge_changes = 0
    for edge_1, edge_2 in edge_maps:
        if edge_1 is None or edge_2 is None:
            num_edge_changes += 1
        elif not edge_match(G1.edges[edge_1], G2.edges[edge_2]):
            num_edge_changes += 1

    # 2: list node changes
    node_subs = []
    node_indels = 0
    for node_1, node_2 in node_maps:
        if node_1 is None or node_2 is None:
            node_indels += 1
        elif not node_match(G1.nodes[node_1], G2.nodes[node_2]):
            node_subs.append((G1.nodes[node_1]['name'], G2.nodes[node_2]['name']))

    return EditInfo(num_edge_changes, node_indels, node_subs)

def graph_edits(G1: nx.DiGraph,
                G2: nx.DiGraph,
                max_dist: int = None) -> Optional[EditInfo]:

    edit_dist = nx.graph_edit_distance(G1, G2,
                                       node_match=node_match,
                                       edge_match=edge_match,
                                       upper_bound=max_dist)
    if edit_dist is None:
        return None
    edit_info = next(nx.optimize_edit_paths(G1, G2,
                                            node_match=node_match,
                                            edge_match=edge_match,
                                            upper_bound=edit_dist))
    return make_edits(G1, G2, edit_info)


def compute_compound_similarity(rban_records: List[Parsed_rBAN_Record],
                                max_dist: int = 3):
    print('Computing compound similarity')
    rban_graphs = [(rban_record.compound_id, rban_record.to_nx_monomer_graph())
                   for rban_record in rban_records]
    similarity_table = {}
    for i, (compound_id1, graph1) in enumerate(rban_graphs):
        print(f'{i}/{len(rban_graphs)}')
        similarity_table[compound_id1] = []
        for compound_id2, graph2 in rban_graphs:
            if compound_id1 == compound_id2:
                continue
            if (edits := graph_edits(graph1, graph2, max_dist=3)) is not None:
                similarity_table[compound_id1].append((compound_id2, edits.to_dict()))

        write_yaml(similarity_table, Path('similarity_table.yaml'))
    return similarity_table


def compute_compound_info_table(rban_records: List[Parsed_rBAN_Record],
                                monomer_helper: MonomerNamesHelper):
    # q: create pandas table with 3 columns: id, num of recognized monomers, num of nerpa2-supported monomers
    # q: how to count recognized monomers? recognized monomers are those which rban name doesn't start with 'X' or have ':' inside
    # q: how to count nerpa2-supported monomers? nerpa2-supported which residues are not UNKNOWN_RESIDUE
    print('Computing compound info table')
    table = []
    for rban_record in rban_records:
        recognized_monomers = 0
        nerpa2_supported_monomers = 0
        for monomer in rban_record.monomers.values():
            rban_name = monomer.name
            if not rban_name.startswith('X') and ':' not in rban_name:
                recognized_monomers += 1
            if monomer_helper.parsed_name(rban_name, name_format='norine').residue != UNKNOWN_RESIDUE:
                nerpa2_supported_monomers += 1
        table.append({'id': rban_record.compound_id,
                      'recognized_monomers': recognized_monomers,
                      'nerpa2_supported_monomers': nerpa2_supported_monomers})
    with open('compound_info_table.tsv', 'w') as f:
        f.write('id\trecognized_monomers\tnerpa2_supported_monomers\n')
        for row in table:
            f.write(f"{row['id']}\t{row['recognized_monomers']}\t{row['nerpa2_supported_monomers']}\n")


def num_recognized_monomers(rban_record: Parsed_rBAN_Record):
    recognized_monomers = 0
    for monomer in rban_record.monomers.values():
        rban_name = monomer.name
        if not rban_name.startswith('X') and ':' not in rban_name:
            recognized_monomers += 1
    return recognized_monomers


def filter_mibig(bgc_variants: List[BGC_Variant],
                 rban_records: List[Parsed_rBAN_Record]):
    rban_records_by_bgc_id = defaultdict(list)
    for rban_record in rban_records:
        bgc_id = rban_record.compound_id.split('.')[0]
        rban_records_by_bgc_id[bgc_id].append(rban_record)

    filtered_bgcs = set()
    for bgc_variant in bgc_variants:
        if len(bgc_variant.modules) < 3:
            continue
        if all(num_recognized_monomers(rban_record) < 3
               for rban_record in rban_records_by_bgc_id[bgc_variant.genome_id]):
            continue
        filtered_bgcs.add(bgc_variant.genome_id)

    with open('filtered_bgcs.txt', 'w') as f:
        for bgc_id in filtered_bgcs:
            f.write(f'{bgc_id}\n')


def compute_compound_isomorphism(rban_records: List[Parsed_rBAN_Record]):
    print('Computing compound isomorphism')
    rban_graphs = [(rban_record.compound_id, rban_record.to_nx_monomer_graph())
                     for rban_record in rban_records]
    isomorphism_table = {}
    for i, (compound_id1, graph1) in enumerate(rban_graphs):
        print(f'{i}/{len(rban_graphs)}')
        isomorphism_table[compound_id1] = []
        for compound_id2, graph2 in rban_graphs:
            if compound_id1 == compound_id2:
                continue
            if Parsed_rBAN_Record.graphs_isomorphic(graph1, graph2):
                isomorphism_table[compound_id1].append(compound_id2)

        write_yaml(isomorphism_table, Path('isomorphism_table_unknown_chr.yaml'))


def main(log: NerpaLogger):  # log is passed as an argument to make it easier to write log in case of exception
    pipeline_helper = PipelineHelper(log)

    nrp_variants, rban_records = pipeline_helper.get_nrp_variants_and_rban_records()
    #compute_compound_isomorphism(rban_records)
    #compute_compound_similarity(rban_records)
    #exit(0)

    bgc_variants = pipeline_helper.pipeline_helper_antismash.get_bgc_variants()
    #filter_mibig(bgc_variants, rban_records)

    hmms = pipeline_helper.construct_hmms(bgc_variants)
    #for hmm in hmms:
    #    hmm.draw(Path(f'{hmm.bgc_variant.genome_id}.png'))

    # compute_compound_info_table(rban_records, pipeline_helper.monomer_names_helper)
    # exit(0)
    nrp_linearizations = pipeline_helper.get_nrp_linearizations(nrp_variants)

    if pipeline_helper.args.only_preprocessing:
        matches = []
    else:
        matches = pipeline_helper.get_matches(hmms, nrp_linearizations)
    pipeline_helper.write_results(matches, bgc_variants, nrp_variants, rban_records)


if __name__ == "__main__":
    log = NerpaLogger()
    try:
        main(log)
    except Exception as e:
        _, exc_value, _ = sys.exc_info()
        log.exception(exc_value)
    finally:
        # TODO: clean up: remove all intermediate files
        pass
