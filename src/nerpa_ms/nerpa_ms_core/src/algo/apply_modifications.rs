use crate::data_types::{common_types::MonomerIdx, graph_modifications::GraphModification, monomer_graph::MonomerGraph, monomers_db::MonomersDB_Entry, parsed_rban_record::Parsed_rBAN_Record};


#[derive(Debug, Clone)]
pub struct AlteredMonomerGraph {
    pub new_monomer_graph: MonomerGraph,
    pub old_to_new_mon_map: Vec<(Option<MonomerIdx>, Option<MonomerIdx>)>,
}

pub fn apply_modifications(monomer_graph: &MonomerGraph,
			   modifications: &[GraphModification]) -> AlteredMonomerGraph {
    let mut new_monomer_graph = (*monomer_graph).clone();
    let mut old_to_new_mon_map: Vec<(Option<MonomerIdx>, Option<MonomerIdx>)> = Vec::new();
    for modification in modifications {
	match modification {
	    GraphModification::Substitute { monomer_idx, mon_db_entry } => {
		new_monomer_graph.substitute(*monomer_idx, mon_db_entry);
		old_to_new_mon_map.push((Some(*monomer_idx), Some(*monomer_idx)));
	    },
	    _ => unimplemented!("Only substitution modifications are supported for now"),
	}
    }

    AlteredMonomerGraph {
	new_monomer_graph,
	old_to_new_mon_map,
    }
}
