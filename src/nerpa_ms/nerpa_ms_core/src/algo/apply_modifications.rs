use std::collections::HashSet;

use crate::data_types::{common_types::MonomerIdx, dag::{DAG, Edge, VertexId}, graph_modifications::{GraphModification, InsertionSite}, monomer_graph::MonomerGraph, monomers_db::{MonomersDB, MonomersDB_Entry}, parsed_rban_record::Parsed_rBAN_Record};


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
    let indices_removed: HashSet<MonomerIdx> = {
	modifications.iter()
	    .filter_map(|m| {
		if let GraphModification::Remove { monomer_idx } = m {
		    Some(*monomer_idx)
		} else {
		    None
		}
	    })
	    .collect()
    };

    let indices_that_should_remain: HashSet<MonomerIdx> = {
	modifications.iter()
	    .map(|m| {
		match m {
		    GraphModification::KeepAsIs { monomer_idx } => vec![*monomer_idx],
		    GraphModification::Substitute { monomer_idx, .. } => vec![*monomer_idx],
		    GraphModification::Insert { site, mon_db_entry } => {
			match site {
			    InsertionSite::Edge(mon_idx1, mon_idx2) => vec![*mon_idx1, *mon_idx2],
			    InsertionSite::Leaf(mon_idx) => vec![*mon_idx],
			}
		    },
		    _ => vec![],
		}
	    })
	    .flatten()
	    .collect()
	};

    indices_removed.is_disjoint(&indices_that_should_remain)
		    
}

enum ModificationsBatch<'batch: 'mon_db, 'mon_db> {
    EdgeInserts(Vec<&'batch GraphModification<'mon_db>>),
    Other(&'batch GraphModification<'mon_db>),
}

pub fn apply_modifications(
    monomer_graph: &MonomerGraph,
    modifications: &[GraphModification],
    monomers_db: &MonomersDB,
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
    println!("\nCreating a new monomer graph by applying the following modifications:\n{}", modifications_str);

    if !modifications_consistent(modifications) {
	println!("Incompatible modifications detected. No monomer graph can be produced.");
	return None;
    }

    // We have to group insertions to the the same edge and deal with them together not to mess up the indices
    let batches = {
	let mut batches: Vec<ModificationsBatch> = Vec::new();
	let mut edge_inserts: Vec<&GraphModification> = Vec::new();
	for m in modifications {
	    match m {
		GraphModification::Insert {
		    site: InsertionSite::Edge { .. },
		    ..
		} => edge_inserts.push(m),
		_ => {
		    if !edge_inserts.is_empty() {
			batches.push(ModificationsBatch::EdgeInserts(edge_inserts.clone()));
			edge_inserts.clear();
		    }
		    batches.push(ModificationsBatch::Other(m));
		}
	    }
	}

	if !edge_inserts.is_empty() {
	    batches.push(ModificationsBatch::EdgeInserts(edge_inserts.clone()));
	}
	 batches
    };
	
    for batch in batches {
	match batch {
	    ModificationsBatch::EdgeInserts(edge_inserts) => {
		// all inserts in the batch should affect the same edge
		let mut edges: HashSet<(MonomerIdx, MonomerIdx)> = HashSet::new();
		let mut db_entries: Vec<&MonomersDB_Entry> = Vec::new();
		for m in edge_inserts {
		    match m {
			GraphModification::Insert {
			    site: InsertionSite::Edge(mon_idx1, mon_idx2),
			    mon_db_entry: db_entry,
			} => {
			    edges.insert((*mon_idx1, *mon_idx2));
			    db_entries.push(*db_entry);
			},
			_ => unreachable!("expected only edge inserts in this batch"),
		    }
		}
		if edges.len() != 1 {
		    panic!("Expected all edge inserts in the batch to affect the same edge, but found inserts for edges: {:?}", edges);
		}
		let edge = edges.into_iter().next().unwrap();
		let mut last_inserted_idx = edge.0;
		for db_entry in db_entries {
		    let new_idx = new_monomer_graph.insert_between(last_inserted_idx, edge.1, db_entry);
		    old_to_new_mon_map.push((None, Some(new_idx)));
		    linearization.push(new_idx);
		    last_inserted_idx = new_idx;
		}
	    },
	    ModificationsBatch::Other(modification) => {
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
			    InsertionSite::Edge{ .. } => panic!("..."),
			    InsertionSite::Leaf(mon_idx) => {
				let new_idx = new_monomer_graph.attach_to_leaf(
				    *mon_idx,
				    mon_db_entry,
				    monomers_db
				);
				old_to_new_mon_map.push((None, Some(new_idx)));
				linearization.push(new_idx);
			    }
			}
		    }
		}
	    }
	}		    
    }

    Some(AlteredMonomerGraph {
	new_monomer_graph,
	linearization,
	old_to_new_mon_map,
    })
}
