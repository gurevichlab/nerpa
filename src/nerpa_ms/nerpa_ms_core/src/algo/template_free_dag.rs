use std::collections::HashMap;

use crate::data_types::common_types::MonomerIdx;
use crate::data_types::dag::{DAG, Edge, VertexLabel, VertexId};
use crate::data_types::graph_modifications::{GraphModification, InsertionSite};
use crate::data_types::monomer_graph::{MonomerFeatures, MonomerGraph, Monomer};
use crate::data_types::monomers_db::{MonomersDB, MonomersDB_Entry};
use crate::data_types::bonds::{Bond, BindingSitesProfile};
use crate::data_types::parsed_rban_record::{AtomId, NRP_Metadata};
use crate::data_types::bond_consts::{AMINO_BINDING_SITE_C, AMINO_BINDING_SITE_N, AMINO_C_END_PROFILE, AMINO_CN_ESTER_O_PROFILE, AMINO_MIDDLE_PROFILE, AMINO_N_END_PROFILE, AMINO_N_ESTER_C_PROFILE, ESTER_BINDING_SITE_C, ESTER_BINDING_SITE_O};
use crate::algo::monomer_graph_operations::splice_bonds;

#[derive(Debug, Clone, Copy)]
pub enum GraphType {
	Linear,
	Cyclic,
	Rho(usize),
}

impl std::fmt::Display for GraphType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GraphType::Linear => write!(f, "Linear"),
            GraphType::Cyclic => write!(f, "Cyclic"),
            GraphType::Rho(rho_idx) => write!(f, "Rho({})", rho_idx),
        }
    }
}

fn get_any_db_entry<'a>(
    profile: &BindingSitesProfile,
    monomers_db: &'a MonomersDB,
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
	let mut fixed_entry = (*entry).clone();
	fixed_entry.set_monomer_idx(MonomerIdx(i as u32));
	let max_id = if i == 0 {
	    AtomId(0)
	} else {
	    fixed_entries
		.last().unwrap()
		.monomer
		.atoms.iter()
		.map(|atom| atom.id.clone())
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

    pub fn build_some_graph_with_linearization(
	graph_type: &GraphType,
	n: usize,
	monomers_db: &MonomersDB,
	compound_id: String,
	metadata: NRP_Metadata,
	) -> (MonomerGraph, Vec<MonomerIdx>) {
	let monomer_graph = {
	    match *graph_type {
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
	};

	let linearization: Vec<MonomerIdx> = (0..n).map(|i| MonomerIdx(i as u32)).collect();
	(monomer_graph, linearization)
    }
}

pub fn create_template_free_dag<'mon_db>(
    graph_type: &GraphType,
    n: usize,
    monomers_db: &'mon_db MonomersDB
) -> DAG<'mon_db> {
    let (monomer_graph, linearization) = MonomerGraph::build_some_graph_with_linearization(
	&graph_type,
	n,
	monomers_db,
	format!("{}_{}", &graph_type, n),
	NRP_Metadata::default(),
    );
    
	
    // forbid/limit certain nodes/edges for debugging purposes
    // let max_subs = 2;
    let max_subs = usize::MAX;

    // HashMap instead of Vec for less headache
    // the keys are actually continuous 0,1,...
    let mut labels: HashMap<usize, VertexLabel> = HashMap::new();
    let mut out_edges: HashMap<usize, Vec<Edge>> = HashMap::new();
    let mut subgraph_root = 0usize;

    for monomer_idx in linearization.iter().cloned() {
	labels.insert(subgraph_root,
		      VertexLabel{
			  monomer_code: None,
			  name: format!("mon{}", monomer_idx.0),
		      });
	out_edges.insert(subgraph_root, Vec::new());


        let subs: Vec<GraphModification> = {
            monomer_graph.possible_substitutions(monomer_idx, monomers_db, true)
		.into_iter()
                .take(max_subs)
                .map(|entry| GraphModification::Substitute {
		    monomer_idx,
		    mon_db_entry: entry,
		})
                .collect()
        };

	// +2 instead of +1 because of the "no modifications" node
        let next_subgraph_root = subgraph_root + subs.len() + 2;

	// ===== Substitutions
        for (i, gm) in subs.into_iter().enumerate() {
	    let sub_node_idx = subgraph_root + 1 + i;
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
                weight: 0,
                modification: Some(gm),
            });
            out_edges.insert(sub_node_idx, vec![Edge {
                to: next_subgraph_root,
                weight: 0,
                modification: None,
            }]);
        }

        subgraph_root = next_subgraph_root;
    }

    let final_node_idx = subgraph_root;
    labels.insert(final_node_idx, VertexLabel{
	monomer_code: None,
	name: "FINAL".to_string(),
    });
    out_edges.insert(final_node_idx, Vec::new());

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
