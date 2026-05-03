use crate::data_types::{common_types::MonomerIdx, graph_modifications::{GraphModification, InsertionSite}, monomer_graph::MonomerGraph, monomers_db::{MonomersDB, MonomersDB_Entry}, parsed_rban_record::Parsed_rBAN_Record};


#[derive(Debug, Clone)]
pub struct AlteredMonomerGraph {
    pub new_monomer_graph: MonomerGraph,
    pub old_to_new_mon_map: Vec<(Option<MonomerIdx>, Option<MonomerIdx>)>,
}

pub fn apply_modifications(
    monomer_graph: &MonomerGraph,
    modifications: &[GraphModification],
    monomers_db: &MonomersDB,
) -> AlteredMonomerGraph {
    let mut new_monomer_graph = (*monomer_graph).clone();
    let mut old_to_new_mon_map: Vec<(Option<MonomerIdx>, Option<MonomerIdx>)> = Vec::new();
    for modification in modifications {
	match modification {
	    GraphModification::Substitute { monomer_idx, mon_db_entry } => {
		new_monomer_graph.substitute(*monomer_idx, mon_db_entry);
		old_to_new_mon_map.push((Some(*monomer_idx), Some(*monomer_idx)));
	    },
	    GraphModification::Remove { monomer_idx } => {
		new_monomer_graph.remove(*monomer_idx, monomers_db);
		old_to_new_mon_map.push((Some(*monomer_idx), None));
	    },
	    GraphModification::Insert { site, mon_db_entry } => {
		let new_idx;
		match site {
		    InsertionSite::Edge(mon_idx1, mon_idx2) => {
			new_idx = new_monomer_graph.insert_between(
			    *mon_idx1,
			    *mon_idx2,
			    mon_db_entry
			);
		    },
		    InsertionSite::Leaf(mon_idx) => {
			new_idx = new_monomer_graph.attach_to_leaf(
			    *mon_idx,
			    mon_db_entry,
			    monomers_db
			);
		    },
		}
		old_to_new_mon_map.push((None, Some(new_idx)));

	    }
	}
    }

    AlteredMonomerGraph {
	new_monomer_graph,
	old_to_new_mon_map,
    }
}
