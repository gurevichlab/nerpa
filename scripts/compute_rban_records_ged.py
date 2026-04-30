from __future__ import annotations
from typing import List, Optional, Tuple
from src.rban_parsing.rban_parser import Parsed_rBAN_Record

import argparse
from collections import Counter
from pathlib import Path
from typing import Any

import networkx as nx
import yaml

def rban_record_to_graph(record: Parsed_rBAN_Record) -> nx.Graph:
    g: nx.Graph = nx.Graph()

    # Keep monomer indices as node IDs (unique), but "label" them with the monomer name.
    for mon_id, mon_info in record.monomers.items():
        g.add_node(mon_id, label=mon_info.name)

    # Undirected monomer connectivity (multiple bonds collapse in a plain Graph).
    for (u, v) in record.monomer_bonds.keys():
        g.add_edge(u, v)

    return g


def load_graphs_by_id(yaml_file: Path) -> List[Tuple[str, nx.Graph]]:
    data = yaml.safe_load(yaml_file.read_text())
    if not isinstance(data, list) or not all(isinstance(x, dict) for x in data):
        raise SystemExit('YAML must be a list of dicts')

    rban_records = [
        Parsed_rBAN_Record.from_dict(data_item)
        for data_item in data
    ]
    return [
        (record.compound_id, rban_record_to_graph(record))
        for record in rban_records
    ]
    

LABEL_KEY = 'label'

def _node_label_hist(g: nx.Graph, label_key: str) -> Counter[Any]:
    return Counter(
        data.get(label_key)
        for _, data in g.nodes(data=True)
    )


def ged(
    g1: nx.Graph,
    g2: nx.Graph,
    max_allowed_distance: Optional[int],
    label_key: str,
) -> Optional[int]:
    """
    Exact GED with unit costs, but return None if GED > max_allowed_distance.

    Model:
      - node substitution (relabel) cost = 0 iff labels match else 1
      - node insert/delete cost = 1
      - edge substitution cost = 0 (unlabeled)
      - edge insert/delete cost = 1

    Cheap pre-pruning (valid lower bounds for this model):
      - node lower bound: required node ins/del + best-case relabels
      - edge lower bound: |m1 - m2| because substitutions don't change edge count
    """
    if max_allowed_distance is None:
        max_allowed_distance = (
              g1.number_of_nodes()
            + g1.number_of_edges()
            + g2.number_of_nodes()
            + g2.number_of_edges()
        )

    # Lower bound for node edits:
    # In the best case, we match as many node labels as possible.
    # Everything else becomes either a relabel (substitution) or an ins/del if sizes differ.
    n1 = _node_label_hist(g1, label_key=label_key)
    n2 = _node_label_hist(g2, label_key=label_key)
    max_label_matches = sum(min(n1[l], n2[l])
                            for l in (n1.keys() | n2.keys()))
    min_node_changes = max(g1.number_of_nodes(), g2.number_of_nodes()) - max_label_matches

    # Lower bound for edge edits:
    # With unlabeled edges and substitution cost 0, the only way to change edge count is ins/del.
    min_edges_changes = abs(g1.number_of_edges() - g2.number_of_edges())

    if min_node_changes + min_edges_changes > max_allowed_distance:
        return None

    def node_subst_cost(a: dict[str, Any], b: dict[str, Any]) -> float:
        return 0.0 if a.get(label_key) == b.get(label_key) else 1.0

    def node_del_cost(_: dict[str, Any]) -> float:
        return 1.0

    def node_ins_cost(_: dict[str, Any]) -> float:
        return 1.0

    def edge_subst_cost(_: dict[str, Any], __: dict[str, Any]) -> float:
        return 0.0

    def edge_del_cost(_: dict[str, Any]) -> float:
        return 1.0

    def edge_ins_cost(_: dict[str, Any]) -> float:
        return 1.0

    d = nx.graph_edit_distance(
        g1,
        g2,
        node_subst_cost=node_subst_cost,
        node_del_cost=node_del_cost,
        node_ins_cost=node_ins_cost,
        edge_subst_cost=edge_subst_cost,
        edge_del_cost=edge_del_cost,
        edge_ins_cost=edge_ins_cost,
        upper_bound=float(max_allowed_distance),  # prunes anything that can't end up <= k
        timeout=None,
    )

    return (
        round(d)
        if d is not None
        else None
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description='Compute pairwise thresholded graph edit distances (GED) for Parsed_rBAN_Record objects.')
    p.add_argument('--input', type=Path, help='YAML file containing a list of Parsed_rBAN_Record objects (as dicts)')
    p.add_argument('--max-allowed-distance',
                   type=int,
                   default=3,
                   help='All distances exceeding that are reported as "-1" (for speed up)')
    p.add_argument('-o',
                   type=Path,
                   help='output file')
                   
    return p.parse_args()


def main() -> None:
    args = parse_args()

    graphs_by_id = load_graphs_by_id(args.input)

    distances = [[] * len(graphs_by_id) for _ in range(len(graphs_by_id))]
    for i, (id1, g1) in enumerate(graphs_by_id):
        for j in range(i):
            id2, g2 = graphs_by_id[j]
            d = ged(g1, g2, max_allowed_distance=args.max_allowed_distance, label_key=LABEL_KEY)
            if d is None:
                d = -1
            distances[i].append(d)

        distances[i].append(0)  # distance to self

    for i in range(len(graphs_by_id)):
        for j in range(i + 1, len(graphs_by_id)):
            distances[i].append(distances[j][i])  # symmetric
            
    with open(args.o, 'w') as out:
        out.write('\t'.join(compound_id for compound_id, _ in graphs_by_id) + '\n')
        for i in range(len(graphs_by_id)):
            out.write(graphs_by_id[i][0] + '\t')
            out.write('\t'.join(map(str, distances[i])) + '\n')


if __name__ == '__main__':
    main()
