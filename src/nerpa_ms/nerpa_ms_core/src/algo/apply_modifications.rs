use crate::data_types::{common_types::MonomerIdx, config::DebugConfig, mod_graph::{ModGraph, Edge, VertexId}, graph_modifications::{GraphModification, InsertionSite}, monomer_graph::MonomerGraph, monomers_db::{MonomersDB, MonomersDB_Entry}, parsed_rban_record::Parsed_rBAN_Record};

use std::collections::{HashMap, HashSet};
use chrono::Local;


#[derive(Debug, Clone)]
pub struct AlteredMonomerGraph {
    pub new_monomer_graph: MonomerGraph,
    pub linearization: Vec<MonomerIdx>,

    // (old_idx, new_idx) pairs
    // old_idx == None --> inserted monomer
    // new_idx == None --> removed monomer
    // neither is None --> substituted monomer
    // Monomers kept as is are not included!
    pub old_to_new_mon_map: Vec<(Option<MonomerIdx>, Option<MonomerIdx>)>,
}

pub fn modifications_consistent(modifications: &[GraphModification]) -> bool {
    let mods_by_index = {
	let mut index_to_mods: HashMap<MonomerIdx, Vec<&GraphModification>> = HashMap::new();
	for m in modifications {
	    match m {
		GraphModification::KeepAsIs { monomer_idx } | GraphModification::Substitute { monomer_idx, .. } | GraphModification::Remove { monomer_idx } => {
		    index_to_mods.entry(*monomer_idx).or_default().push(m);
		},
		GraphModification::Insert { site, .. } => {
		    match site {
			InsertionSite::Edge(mon_idx1, mon_idx2) => {
			    index_to_mods.entry(*mon_idx1).or_default().push(m);
			    index_to_mods.entry(*mon_idx2).or_default().push(m);
			},
			InsertionSite::Leaf(mon_idx) => {
			    index_to_mods.entry(*mon_idx).or_default().push(m);
			}
		    }
		}
	    }
	}
	index_to_mods
    };

    for (_mon_idx, mods) in mods_by_index {
	// if has remove, shouldn't have anything else
	// only one modification per monomer_idx except for inserts on edges
	let has_remove = mods.iter().any(|m| matches!(m, GraphModification::Remove { .. }));
	if has_remove && mods.len() > 1 {
	    return false;
	}
	let num_not_insert = {
	    mods.iter()
		.filter(|m| !matches!(m, GraphModification::Insert { .. }))
		.count()
	};
	if num_not_insert > 1 {
	    return false;
	}
    }

    true
}

enum ModificationsBatch<'batch: 'mon_db, 'mon_db> {
    EdgeInserts(Vec<&'batch GraphModification<'mon_db>>),
    Other(&'batch GraphModification<'mon_db>),
}

pub fn apply_modifications(
    monomer_graph: &MonomerGraph,
    modifications: &[GraphModification],
    monomers_db: &MonomersDB,
    debug_stdout: bool,
) -> Option<AlteredMonomerGraph> {
    let mut new_monomer_graph = (*monomer_graph).clone();
    let mut old_to_new_mon_map: Vec<(Option<MonomerIdx>, Option<MonomerIdx>)> = Vec::new();
    let mut linearization: Vec<MonomerIdx> = Vec::new();

    let modifications_str = {
	modifications.iter()
        .filter(|m| {!matches!(m, GraphModification::KeepAsIs { .. })})
	.map(|m| m.to_str_short())
	.collect::<Vec<_>>()
	.join("\n")
    };

    if debug_stdout {
	let now = Local::now().format("%Y-%m-%d %H:%M:%S%.3f");
	println!(
	    "\n[{}] Creating a new monomer graph by applying the following modifications:\n{}",
	    now,
	    modifications_str
	);
    }

    if !modifications_consistent(modifications) {
	if debug_stdout {
	    println!("Incompatible modifications detected. No monomer graph can be produced.");
	}
	return None;
    }

    // We have to group insertions to the the same edge and deal with them together not to mess up the indices
    let non_insert_mods: Vec<&GraphModification> = {
	modifications
	    .iter()
	    .filter(|m| {
		!matches!(m, GraphModification::Insert { site: InsertionSite::Edge { .. }, .. })
	    })
	    .collect()
    };
    let insert_by_edge: HashMap<(MonomerIdx, MonomerIdx), Vec<&GraphModification>> = {
	let mut edge_to_inserts: HashMap<(MonomerIdx, MonomerIdx), Vec<&GraphModification>> = HashMap::new();
	for m in modifications {
	    if let GraphModification::Insert { site: InsertionSite::Edge(mon_idx1,
									 mon_idx2), .. } = m {
		edge_to_inserts.entry((*mon_idx1, *mon_idx2)).or_default().push(m);
	    }
	}
	edge_to_inserts
    };
	
    // Apply edge insertions
    for (edge, edge_inserts) in insert_by_edge {
	let db_entries: Vec<&MonomersDB_Entry> = {
	    edge_inserts
		.iter()
		.map(|m| {
		    match m {
			GraphModification::Insert { site: _, mon_db_entry } => *mon_db_entry,
			_ => unreachable!("expected only edge inserts in this batch"),
		    }
		})
		.collect()
	};
	
	let mut last_inserted_idx = edge.0;
	for (i, db_entry) in db_entries.iter().enumerate() {
	    let insert_res = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
		new_monomer_graph.insert_between(last_inserted_idx, edge.1, db_entry)
	    }));
	    let new_idx = match insert_res {
		Ok(idx) => idx,
		Err(_panic) => {
		    eprintln!(
			"insert_between panicked for edge ({}, {}) originally modification {:?}, {}-th insert in the batch", 
			last_inserted_idx, edge.1,
			edge_inserts[i],
			i + 1,
		    );
		    panic!("insert_between panicked");
		}
	    };

	    old_to_new_mon_map.push((None, Some(new_idx)));
	    linearization.push(new_idx);
	    last_inserted_idx = new_idx;
	};
    }

    // Apply other modifications
    for modification in non_insert_mods {
	match modification {
	    GraphModification::KeepAsIs { monomer_idx } => {
		linearization.push(*monomer_idx);
	    },
	    GraphModification::Substitute { monomer_idx, mon_db_entry } => {
		new_monomer_graph.substitute(*monomer_idx, mon_db_entry);
		old_to_new_mon_map.push((Some(*monomer_idx), Some(*monomer_idx)));
		linearization.push(*monomer_idx);
	    },
	    GraphModification::Remove { monomer_idx } => {
		new_monomer_graph.remove(*monomer_idx, monomers_db);
		old_to_new_mon_map.push((Some(*monomer_idx), None));
	    },
	    GraphModification::Insert { site, mon_db_entry } => {
		match site {
		    InsertionSite::Edge{ .. } => unreachable!("edge insertions should have been handled in a separate batch"),
		    InsertionSite::Leaf(mon_idx) => {
			// let new_idx = new_monomer_graph.attach_leaf(
			//     *mon_idx,
			//     mon_db_entry,
			//     monomers_db
			// );
			// old_to_new_mon_map.push((None, Some(new_idx)));
			// linearization.push(new_idx);
			unimplemented!("leaf insertions are not supported for now")
		    }
		}
	    }
	}
    }

    if debug_stdout {
	let now = Local::now().format("%Y-%m-%d %H:%M:%S%.3f");
	println!(
	    "\n[{}] New monomer graph created.\n",
	    now,
	);
    }

    Some(AlteredMonomerGraph {
	new_monomer_graph,
	linearization,
	old_to_new_mon_map,
    })
}
