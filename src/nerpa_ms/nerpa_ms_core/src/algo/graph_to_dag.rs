use std::collections::HashMap;

use crate::data_types::bond_consts::{AMINO_BINDING_SITE_C, AMINO_BINDING_SITE_N};
use crate::data_types::bonds::{BindingSiteType, Bond};
use crate::data_types::common_types::MonomerIdx;
use crate::data_types::dag::{DAG, Edge, VertexLabel, VertexId};
use crate::data_types::graph_modifications::{GraphModification, InsertionSite};
use crate::data_types::monomer_graph::{MonomerFeatures, MonomerGraph};
use crate::data_types::{monomers_db::MonomersDB, parsed_rban_record::Parsed_rBAN_Record};

pub fn get_insert_mods_leaf<'a>(
    monomer_graph: &MonomerGraph,
    monomer_idx: MonomerIdx,
    bs_type: &BindingSiteType,
    monomers_db: &'a MonomersDB,
    max_inserts: usize
) -> Vec<GraphModification<'a>> {
    monomer_graph.possible_attachments(
	monomer_idx,
	bs_type,
	monomers_db
    )
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
    bond: &Bond,
    monomers_db: &'a MonomersDB,
    max_inserts: usize
) -> Vec<GraphModification<'a>> {
    monomer_graph.possible_insertions_inside_bond(bond, monomers_db)
	.into_iter()
	.take(max_inserts)
	.map(|entry| GraphModification::Insert {
	    site: InsertionSite::Edge(bond.monomers.0, bond.monomers.1),
	    mon_db_entry: entry,
	})
	.collect()
}

pub fn get_inserts_for_linearization<'a>(
	monomer_graph: &MonomerGraph,
	linearization: &[MonomerIdx],
	monomers_db: &'a MonomersDB,
	max_inserts: usize
) -> Vec<Vec<GraphModification<'a>>> {
    // for linearization of length n returns exactly n+1 vectors of insert modifications:
    // before the first monomer, between each pair of consecutive monomers, and after the last monomer (some may be empty if no insertions are possible)
    let mut inserts_per_position: Vec<Vec<GraphModification>> = Vec::new();

    for (lin_idx, monomer_idx) in linearization.iter().cloned().enumerate() {
	// inserting BEFORE the current monomer
	let mut inserts: Vec<GraphModification> = Vec::new();
	let maybe_prev_idx = if lin_idx > 0 {
	    Some(linearization[lin_idx - 1])
	} else {
	    None
	};

	// 1. Try inserting in the edge between linearization[lin_idx-1] and linearization[lin_idx] (if lin_idx > 0)
	if let Some(prev_idx) = maybe_prev_idx
	    && let Some(bond) = monomer_graph.get_bond(prev_idx, monomer_idx) {
		inserts.extend(
		    get_insert_mods_edge(
			monomer_graph,
			bond,
			monomers_db,
			max_inserts
		    )
		);
		continue;
	    }

	// 2. No C->N bond (prev_idx, mon_idx) but a C->N bond (parent, mon_idx) -- try inserting there
	let cn_bonds_to_mon: Vec<Bond> = {
	    monomer_graph.bonds_by_bs_type(monomer_idx)
		.iter()
		.filter_map(|(bst, bond)|
			    if *bst == *AMINO_BINDING_SITE_N {
				Some(bond)
			    } else {
				None
			    })
		.cloned()
		.collect()
	};

	for bond_to_mon in cn_bonds_to_mon.iter() {
	    inserts.extend(
		get_insert_mods_edge(
		    monomer_graph,
		    bond_to_mon,
		    monomers_db,
		    max_inserts
		)
	    );
	}

	// 3. No C->N bond (prev_idx, mon_idx) but a C->N bond (prev_idx, child) -- try inserting there
	if let Some(prev_idx) = maybe_prev_idx {
	    let cn_bonds_from_prev: Vec<Bond> = {
		monomer_graph.bonds_by_bs_type(prev_idx)
		    .iter()
		    .filter(|(bst, _bond)| *bst == *AMINO_BINDING_SITE_C)
		    .map(|(_bst, bond)| bond)
		    .cloned()
		    .collect()
	    };
	    for bond_from_prev in cn_bonds_from_prev.iter() {
		inserts.extend(
		    get_insert_mods_edge(
			monomer_graph,
			bond_from_prev,
			monomers_db,
			max_inserts,
		    ));
	    }
	}

	// 4. If monomer is a leaf, try attaching a new monomer with an AMINO bond
	// Applied to leaves only not to produce unrealistic graphs
	if monomer_graph.degree(monomer_idx) == 1 {
	    inserts.extend(
		get_insert_mods_leaf(
		    monomer_graph,
		    monomer_idx,
		    &*AMINO_BINDING_SITE_N,
		    monomers_db,
		    max_inserts
		)
	    );
	}


	// 5. If prev_monomer_idx is a leaf, try attaching a new monomer to it with an AMINO bond
	// Applied to leaves only not to produce unrealistic graphs
	if let Some(prev_idx) = maybe_prev_idx
	    && monomer_graph.degree(prev_idx) == 1 {
		inserts.extend(
		    get_insert_mods_leaf(
			monomer_graph,
			prev_idx,
			&*AMINO_BINDING_SITE_C,
			monomers_db,
			max_inserts
		    )
		);
	    }

	inserts_per_position.push(inserts);
    }

    let mut inserts_after_last: Vec<GraphModification> = Vec::new();
    let last = *linearization.last().unwrap();
    let bonds_from_last: Vec<Bond> = {
	monomer_graph.bonds_by_bs_type(last)
	    .iter()
	    .filter_map(|(bst, bond)|
			if *bst == *AMINO_BINDING_SITE_C {
			    Some(bond)
			} else {
			    None
			})
	    .cloned()
	    .collect()
    };
    for bond_from_last in bonds_from_last.iter() {
	inserts_after_last.extend(
	    get_insert_mods_edge(
		monomer_graph,
		bond_from_last,
		monomers_db,
		max_inserts
	    )
	);
    }
    if monomer_graph.degree(last) == 1 {
	inserts_after_last.extend(
	    get_insert_mods_leaf(
		monomer_graph,
		last,
		&*AMINO_BINDING_SITE_C,
		monomers_db,
		max_inserts
	    )
	);
    }
	

    inserts_per_position
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
		weight: 0,
		modification: None,
	    });
	// It's important that the BACK edge has weight 1, as DP assumes that edges with weight 0 increase the vertex index
	out_edges.insert(insert_node_idx, vec![Edge {
	    to: subgraph_root,
	    weight: 1,
	    modification: Some(gm.clone()),
	}]);
    }
	
}

pub fn create_dag<'mon_db>(monomer_graph: &MonomerGraph,
		      linearization: &Vec<MonomerIdx>,
		      monomers_db: &'mon_db MonomersDB) -> DAG<'mon_db> {
    // forbid/limit certain nodes/edges for debugging purposes
    let debug_output = false;
    let max_subs = 2;
    let max_inserts = 0;
    let allow_deletions = true;

    // let debug_output = false;
    // let allow_deletions = true;
    // let max_subs = usize::MAX;
    // let max_inserts = usize::MAX;

    // HashMap instead of Vec for less headache
    // the keys are actually continuous 0,1,...
    let mut labels: HashMap<usize, VertexLabel> = HashMap::new();
    let mut out_edges: HashMap<usize, Vec<Edge>> = HashMap::new();
    let mut subgraph_root = 0usize;

    let inserts_per_position =
	get_inserts_for_linearization(
	    monomer_graph,
	    linearization,
	    monomers_db,
	    max_inserts
	);

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
            monomer_graph.possible_substitutions(monomer_idx, monomers_db, false)
		.into_iter()
                .take(max_subs)
                .map(|entry| GraphModification::Substitute {
		    monomer_idx,
		    mon_db_entry: entry,
		})
                .collect()
        };
	let inserts = &inserts_per_position[lin_idx];

	if debug_output {
	    let inserts_str = {
		let mut strs = inserts.iter()
		    .map(|m| match m {
			GraphModification::Insert { site, mon_db_entry } => {
			    mon_db_entry.monomer.features.name.0.clone()
			},
			_ => unreachable!("expected Insert"),
		    })
		    .collect::<Vec<_>>();

		strs.sort();
		strs.join(", ")
	    };
	    println!("Inserts before monomer {}:\n{}\n", monomer_idx.0, inserts_str);
	}

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
	add_inserts_to_subgraph_root(
	    subgraph_root,
	    &inserts,
	    &mut labels,
	    &mut out_edges
	);

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
            out_edges.insert(
		sub_node_idx,
		vec![Edge {
                    to: next_subgraph_root,
                    weight: 0,
                    modification: None,
		}]);
        }

	// ===== Deletion
	if monomer_graph.can_remove(monomer_idx, monomers_db) && allow_deletions {
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
    let inserts_at_end = inserts_per_position.last().unwrap();
    add_inserts_to_subgraph_root(
	final_node_idx,
	inserts_at_end,
	&mut labels,
	&mut out_edges
    );
    
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
