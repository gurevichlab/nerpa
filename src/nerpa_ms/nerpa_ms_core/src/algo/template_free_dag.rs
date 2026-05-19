use std::collections::HashMap;

use crate::data_types::common_types::MonomerIdx;
use crate::data_types::dag::{DAG, Edge, VertexLabel, VertexId};
use crate::data_types::graph_modifications::{GraphModification, InsertionSite};
use crate::data_types::monomer_graph::{MonomerFeatures, MonomerGraph, Monomer};
use crate::data_types::monomers_db::{MonomersDB, MonomersDB_Entry};
use crate::data_types::bonds::{Bond, BindingSitesProfile};
use crate::data_types::parsed_rban_record::NRP_Metadata;
use crate::data_types::bond_consts::{AMINO_BINDING_SITE_C, AMINO_BINDING_SITE_N, AMINO_C_END_PROFILE, AMINO_CN_ESTER_O_PROFILE, AMINO_MIDDLE_PROFILE, AMINO_N_END_PROFILE, AMINO_N_ESTER_C_PROFILE, ESTER_BINDING_SITE_C, ESTER_BINDING_SITE_O};
use crate::algo::monomer_graph_operations::splice_bonds;

pub enum GraphType {
	Linear,
	Cyclic,
	Rho(usize),
}

fn get_any_db_entry(
    profile: &BindingSitesProfile,
    monomers_db: &MonomersDB,
) -> &'a MonomersDB_Entry {
	// just get any entry from the monomers_db to use as a template for building the monomer graph
	// we will overwrite the monomer_idx and shift the atom ids, so it doesn't matter which entry we take
    monomers_db.get(profile)
	.expect("monomers_db should have at least one entry for the given binding sites profile")
	.first()
	.expect("monomers_db should have at least one entry for the given binding sites profile")
}

fn fix_entries(entries: &[&MonomersDB_Entry]) -> Vec<MonomersDB_Entry> {
    // assign monomer_idx to each entry corresponding to its index in the entries slice, and shift the atom ids so that they are unique across entries
    let mut fixed_entries: Vec<MonomersDB_Entry> = Vec::new();
    for (i, entry) in entries.iter().enumerate() {
	let mut fixed_entry = entry.clone();
	fixed_entry.set_monomer_idx(MonomerIdx(i as u32));
	let max_id = if i == 0 {
	    AtomId(0)
	} else {
	    fixed_entries
		.last().unwrap()
		.monomer
		.atoms.iter()
		.map(|atom| atom.id)
		.max().unwrap()
	};
	fixed_entry.shift_atom_ids(max_id.0 + 1);
	fixed_entries.push(fixed_entry);
    }
    fixed_entries
}

impl MonomerGraph {
    fn build_linear(
	entries: &[MonomersDB_Entry],
	compound_id: String,
	metadata: NRP_Metadata,
    ) -> MonomerGraph {
	// connects entries via C->N chain
	debug_assert!(
	    entries.iter()
		.enumerate()
		.all(|(i, entry)| entry.get_monomer_idx() == MonomerIdx(i as u32)),
	    "All entries should have monomer_idx corresponding to their index in the entries slice"
	);

	let monomers: HashMap<MonomerIdx, Monomer> = {
	    entries
		.iter()
		.enumerate()
		.map(|(i, entry)| (MonomerIdx(i as u32), entry.monomer.clone()))
		.collect()
	};

	let monomer_bonds: Vec<Bond> = {
	    let mut bonds: Vec<Bond> = Vec::new();
	    for i in 0..(monomers.len() - 1) {
		let bond_left = {
		    entries[i]
			.bonds_by_bs
			.get_by_bs_type(&*AMINO_BINDING_SITE_C)
			.expect("monomer is missing C binding site")
		};
		let bond_right = {
		    entries[i + 1]
			.bonds_by_bs
			.get_by_bs_type(&*AMINO_BINDING_SITE_N)
			.expect("monomer is missing N binding site")
		};
		let bond = splice_bonds((MonomerIdx(i as u32), bond_left),
					(MonomerIdx((i + 1) as u32), bond_right));
		bonds.push(bond);
	    }
	    bonds
	};

	MonomerGraph {
	    compound_id,
	    monomers,
	    monomer_bonds,
	    metadata,
	}
    }
	
    pub fn build_cyclic(
	entries: &[MonomersDB_Entry],
	compound_id: String,
	metadata: NRP_Metadata,
    ) -> MonomerGraph {
	// connects entries via C->N chain, and also connects the last entry back to the first to form a cycle
	debug_assert!(
	    entries.iter()
		.enumerate()
		.all(|(i, entry)| entry.get_monomer_idx() == MonomerIdx(i as u32)),
	    "All entries should have monomer_idx corresponding to their index in the entries slice"
	);
	let mut monomer_graph = MonomerGraph::build_linear(entries, compound_id, metadata);

	let bond_left = {
	    entries.last().unwrap()
		.bonds_by_bs
		.get_by_bs_type(&*AMINO_BINDING_SITE_C)
		.expect("monomer is missing C binding site")
	};
	let bond_right = {
	    entries.first().unwrap()
		.bonds_by_bs
		.get_by_bs_type(&*AMINO_BINDING_SITE_N)
		.expect("monomer is missing N binding site")
	};
	let cycle_bond = splice_bonds(
	    (MonomerIdx((entries.len() - 1) as u32), bond_left),
	    (MonomerIdx(0), bond_right)
	);

	monomer_graph.monomer_bonds.push(cycle_bond);
	monomer_graph
    }

    fn build_rho(
	entries: &[MonomersDB_Entry],
	rho_idx: usize,
	compound_id: String,
	metadata: NRP_Metadata,
    ) -> MonomerGraph {
	// 0 -> (C->N) -> 1 -> ... -> (C->N) -> n-1 -> (C->O) -> rho_idx
	debug_assert!(
	    rho_idx < entries.len(),
	    "rho_idx should be a valid index in entries"
	);
	debug_assert!(
	    entries.iter()
		.enumerate()
		.all(|(i, entry)| entry.get_monomer_idx() == MonomerIdx(i as u32)),
	    "All entries should have monomer_idx corresponding to their index in the entries slice"
	);
	let mut monomer_graph = MonomerGraph::build_linear(entries, compound_id, metadata);

	let bond_left = {
	    entries.last().unwrap()
		.bonds_by_bs
		.get_by_bs_type(&*ESTER_BINDING_SITE_C)
		.expect("monomer is missing C binding site")
	};
	let bond_right = {
	    entries[rho_idx]
		.bonds_by_bs
		.get_by_bs_type(&*ESTER_BINDING_SITE_O)
		.expect("monomer is missing O binding site")
	};
	let rho_bond = splice_bonds(
	    (MonomerIdx((entries.len() - 1) as u32), bond_left),
	    (MonomerIdx(rho_idx as u32), bond_right)
	);

	monomer_graph.monomer_bonds.push(rho_bond);
	monomer_graph
    }

    pub fn build_some_graph(
	graph_type: GraphType,
	n: usize,
	monomers_db: &MonomersDB,
	compound_id: String,
	metadata: NRP_Metadata,
	) -> MonomerGraph {
	match graph_type {
	    GraphType::Linear => {
		let mut entries_raw: Vec<&MonomersDB_Entry> = Vec::new();
		entries_raw.push(get_any_db_entry(&*AMINO_C_END_PROFILE, monomers_db));
		for _ in 1..(n-1) {
		    entries_raw.push(get_any_db_entry(&*AMINO_MIDDLE_PROFILE, monomers_db));
		}
		entries_raw.push(get_any_db_entry(&*AMINO_N_END_PROFILE, monomers_db));
		let entries = fix_entries(&entries_raw);
		MonomerGraph::build_linear(&entries, compound_id, metadata)
	    },
	    GraphType::Cyclic => {
		let entries_raw: Vec<&MonomersDB_Entry> = (0..n)
		    .map(|_| get_any_db_entry(&*AMINO_MIDDLE_PROFILE, monomers_db))
		    .collect();
		let entries = fix_entries(&entries_raw);
		MonomerGraph::build_cyclic(&entries, compound_id, metadata)
	    },
	    GraphType::Rho(rho_idx) => {
		let mut entries_raw: Vec<&MonomersDB_Entry> = Vec::new();
		entries_raw.push(get_any_db_entry(&*AMINO_C_END_PROFILE, monomers_db));
		for _ in 1..rho_idx {
		    entries_raw.push(get_any_db_entry(&*AMINO_MIDDLE_PROFILE, monomers_db));
		}
		entries_raw.push(get_any_db_entry(&*AMINO_CN_ESTER_O_PROFILE, monomers_db));
		for _ in (rho_idx+1)..(n-1) {
		    entries_raw.push(get_any_db_entry(&*AMINO_MIDDLE_PROFILE, monomers_db));
		}
		entries_raw.push(get_any_db_entry(&*AMINO_N_ESTER_C_PROFILE, monomers_db));

		let entries = fix_entries(&entries_raw);
		MonomerGraph::build_rho(&entries, rho_idx, compound_id, metadata)
	    },
	}
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
	let inserts: Vec<GraphModification> = {
	    let mut inserts: Vec<GraphModification> = Vec::new();
	    if monomer_graph.degree(monomer_idx) == 1 {
		inserts.extend(
		    get_insert_mods_leaf(monomer_graph, monomer_idx, monomers_db, max_inserts)
		);
	    }
	    if let Some(parent) = monomer_graph.get_amino_chain_parent(monomer_idx) {
		inserts.extend(
		    get_insert_mods_edge(monomer_graph, parent, monomer_idx, monomers_db, max_inserts)
		);
	    }
	    if lin_idx > 0 {
		let prev_monomer_idx = linearization[lin_idx - 1];
		if monomer_graph.degree(prev_monomer_idx) == 1 {
		    inserts.extend(
			get_insert_mods_leaf(monomer_graph, prev_monomer_idx, monomers_db, max_inserts)
		    );
		}
		if let Some(prev_child) = monomer_graph.get_amino_chain_child(prev_monomer_idx) {
		    if prev_child != monomer_idx {
			inserts.extend(
			    get_insert_mods_edge(monomer_graph, prev_monomer_idx, prev_child, monomers_db, max_inserts)
			);
		    }
		}
	    }

	    inserts

	};

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
    let last_monomer_idx = linearization.last().unwrap().clone();
    let inserts_at_end: Vec<GraphModification> = {
	if monomer_graph.degree(last_monomer_idx) == 1 {
	    get_insert_mods_leaf(monomer_graph, last_monomer_idx, monomers_db, max_inserts)
	} else {
	    if let Some(child) = monomer_graph.get_amino_chain_child(last_monomer_idx) {
		get_insert_mods_edge(monomer_graph, last_monomer_idx, child, monomers_db, max_inserts)
	    }
	    else {
		Vec::new()
	    }
	}
    };
    add_inserts_to_subgraph_root(final_node_idx, &inserts_at_end, &mut labels, &mut out_edges);
    
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
