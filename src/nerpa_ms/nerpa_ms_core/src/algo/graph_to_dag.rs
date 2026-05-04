use std::collections::HashMap;

use crate::data_types::common_types::MonomerIdx;
use crate::data_types::dag::{DAG, Edge, VertexLabel, VertexId};
use crate::data_types::graph_modifications::{GraphModification, InsertionSite};
use crate::data_types::monomer_graph::{MonomerFeatures, MonomerGraph};
use crate::data_types::{monomers_db::MonomersDB, parsed_rban_record::Parsed_rBAN_Record};

pub fn get_insert_mods_leaf<'a>(
    monomer_graph: &MonomerGraph,
    monomer_idx: MonomerIdx,
    monomers_db: &'a MonomersDB,
    max_inserts: usize
) -> Vec<GraphModification<'a>> {
    monomer_graph.possible_insertions_at_leaf(monomer_idx, monomers_db)
	.into_iter()
	.take(max_inserts)
	.map(|entry| GraphModification::Insert {
	    site: InsertionSite::Leaf(monomer_idx),
	    mon_db_entry: entry,
	})
	.collect()
}

pub fn get_insert_mods_edge<'a>(
	monomer_graph: &MonomerGraph,
	monomer_idx1: MonomerIdx,
	monomer_idx2: MonomerIdx,
	monomers_db: &'a MonomersDB,
	max_inserts: usize
) -> Vec<GraphModification<'a>> {
    monomer_graph.possible_insertions_between(monomer_idx1, monomer_idx2, monomers_db)
	.into_iter()
	.take(max_inserts)
	.map(|entry| GraphModification::Insert {
	    site: InsertionSite::Edge(monomer_idx1, monomer_idx2),
	    mon_db_entry: entry,
	})
	.collect()
}

pub fn add_inserts_to_subgraph_root<'a>(
    subgraph_root: usize,
    insert_mods: &[GraphModification<'a>],
    labels: &mut HashMap<usize, VertexLabel>,
    out_edges: &mut HashMap<usize, Vec<Edge<'a>>>
) {
    let free_node_idx = labels.len();
    for (i, gm) in insert_mods.iter().enumerate() {
	let insert_node_idx = free_node_idx + i;
	match &gm {
	    GraphModification::Insert { site: _, mon_db_entry } => {
		labels.insert(insert_node_idx, VertexLabel {
		    monomer_code: Some(mon_db_entry.monomer.features.mon_code.clone()),
		    name: mon_db_entry.monomer.features.name.0.clone()
		});
	    },
	    _ => unreachable!("expected Insert"),
	}
	out_edges.get_mut(&subgraph_root)
	    .unwrap()
	    .push(Edge {
		to: insert_node_idx,
		weight: 1,
		modification: Some(gm.clone()),
	    });
	out_edges.insert(insert_node_idx, vec![Edge {
	    to: subgraph_root,
	    weight: 0,
	    modification: None,
	}]);
    }
	
}

pub fn create_dag<'mon_db>(monomer_graph: &MonomerGraph,
		      linearization: &Vec<MonomerIdx>,
		      monomers_db: &'mon_db MonomersDB) -> DAG<'mon_db> {
    // for debugging purposes, limit the number of modifications
    let max_subs = 2; 
    let max_inserts = 2;

    // HashMap instead of Vec for less headache
    // the keys are actually continuous 0,1,...
    let mut labels: HashMap<usize, VertexLabel> = HashMap::new();
    let mut out_edges: HashMap<usize, Vec<Edge>> = HashMap::new();
    let mut subgraph_root = 0usize;

    for (lin_idx, monomer_idx) in linearization.iter().cloned().enumerate() {
	let monomer_info: &MonomerFeatures = {
	    &monomer_graph.monomers.get(&monomer_idx)
		.expect("linearization contains monomer index that is not present in rban_record.monomers")
		.features
	};
	labels.insert(subgraph_root,
		      VertexLabel{
			  monomer_code: None,
			  name: format!("mon{}", monomer_idx.0),
		      });
	out_edges.insert(subgraph_root, Vec::new());


        let subs: Vec<GraphModification> = {
            let sub_entries = monomer_graph.possible_substitutions(monomer_idx, monomers_db);
            // filter out substitutions that are identical to the original monomer
            sub_entries.into_iter()
                .take(max_subs)
                .map(|entry| GraphModification::Substitute {
		    monomer_idx,
		    mon_db_entry: entry,
		})
                .collect()
        };
	let inserts: Vec<GraphModification> = {
	    if lin_idx == 0 {
		if monomer_graph.degree(monomer_idx) != 1 {
		    Vec::new()
		} else {
		    get_insert_mods_leaf(monomer_graph, monomer_idx, monomers_db, max_inserts)
		}
	    } else {
		let prev_monomer_idx = linearization[lin_idx - 1];
		if monomer_graph.get_bond(monomer_idx, prev_monomer_idx).is_none() {
		    println!("Warning: no bond between monomers {} and {} in the linearization. Skipping insertions at edge.", monomer_idx.0, prev_monomer_idx.0);
		    Vec::new()
		} else {
		    get_insert_mods_edge(monomer_graph, prev_monomer_idx, monomer_idx, monomers_db, max_inserts)
		}
	    }
	};

	// +2 instead of +1 because of the "no modifications" node
        let next_subgraph_root = subgraph_root + inserts.len() + subs.len() + 2;

	// ===== No modifications (weight 0)
        labels.insert(subgraph_root + 1, VertexLabel {
	    monomer_code: Some(monomer_info.mon_code.clone()),
	    name: monomer_info.name.0.clone() + "*",
	});
        out_edges.get_mut(&subgraph_root)
	    .unwrap()
	    .push(Edge {
		to: subgraph_root + 1,
		weight: 0,
		modification: Some(GraphModification::KeepAsIs { monomer_idx }),
            });
        out_edges.insert(subgraph_root + 1, vec![Edge {
            to: next_subgraph_root,
            weight: 0,
            modification: None,
        }]);
			
	// ===== Insertions
	add_inserts_to_subgraph_root(subgraph_root, &inserts, &mut labels, &mut out_edges);

	// ===== Substitutions
        for (i, gm) in subs.into_iter().enumerate() {
	    let sub_node_idx = subgraph_root + inserts.len() + 2 + i;
            match &gm {
                GraphModification::Substitute { mon_db_entry, .. } => {
                    labels.insert(sub_node_idx, VertexLabel {
			monomer_code: Some(mon_db_entry.monomer.features.mon_code.clone()),
			name: mon_db_entry.monomer.features.name.0.clone()
		    });
                },
                _ => unreachable!("expected Substitute"),
            }

            out_edges.get_mut(&subgraph_root)
                .unwrap()
		.push(Edge {
                to: sub_node_idx,
                weight: 1,
                modification: Some(gm),
            });
            out_edges.insert(sub_node_idx, vec![Edge {
                to: next_subgraph_root,
                weight: 0,
                modification: None,
            }]);
        }

	// ===== Deletion
	if monomer_graph.can_remove(monomer_idx, monomers_db) {
	    let del_mod = GraphModification::Remove { monomer_idx: monomer_idx };
	    out_edges.get_mut(&subgraph_root)
		.unwrap()
		.push(Edge {
		to: next_subgraph_root,
		weight: 1,
		modification: Some(del_mod),
	    });
	}

        subgraph_root = next_subgraph_root;
    }

    let final_node_idx = subgraph_root;
    labels.insert(final_node_idx, VertexLabel{
	monomer_code: None,
	name: "FINAL".to_string(),
    });
    out_edges.insert(final_node_idx, Vec::new());

    // Add insertions at the end of the linearization
    let last_monomer_idx = linearization.last().unwrap().clone();
    if monomer_graph.degree(last_monomer_idx) == 1 {
	let inserts_at_end = get_insert_mods_leaf(monomer_graph, last_monomer_idx, monomers_db, max_inserts);
	add_inserts_to_subgraph_root(final_node_idx, &inserts_at_end, &mut labels, &mut out_edges);
    }
    
    // Convert HashMaps to Vecs. The keys of the HashMaps are actually continuous 0,1,...,labels.len()-1, so we can just iterate over the keys in order.
    let num_nodes = labels.len();
    debug_assert_eq!(out_edges.len(), num_nodes, "out_edges should have the same number of entries as labels");
    let labels_vec: Vec<VertexLabel> = {
	(0..num_nodes)
	    .map(|i| labels[&i].clone())
	    .collect()
    };
    let out_edges_vec: Vec<Vec<Edge>> = {
	(0..num_nodes)
	    .map(|i| out_edges[&i].clone())
	    .collect()
    };

    DAG { 
	nrp_variant_id: monomer_graph.compound_id.clone(),
	labels: labels_vec,
	out_edges: out_edges_vec,
	start: 0,
	finish: final_node_idx,
    }
}
