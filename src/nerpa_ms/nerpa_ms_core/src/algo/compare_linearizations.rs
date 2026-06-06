use std::collections::HashSet;

use crate::algo::graph_to_dag::get_inserts_for_linearization;
use crate::data_types::common_types::MonomerIdx;
use crate::data_types::graph_modifications::GraphModification;
use crate::data_types::monomer_graph::{MonomerGraph};
use crate::data_types::monomers_db::MonomersDB_Entry;
use crate::data_types::parsed_rban_record::NorineMonomerName;
use crate::data_types::{monomers_db::MonomersDB};

#[derive(Debug, Clone)]
pub enum LinearizationModification {
    Substitute {
        pos: usize,
        new_name: NorineMonomerName,
    },
    Remove {
        pos: usize,
    },
    Insert {
	pos: usize,
	new_name: NorineMonomerName,
    },
}

impl <'a> LinearizationModification {
    pub fn try_to_graph_modification(
	&self,
	monomer_graph: &MonomerGraph,
	linearization: &[MonomerIdx],
	monomers_db: &'a MonomersDB,
    ) -> Option<GraphModification<'a>> {
	match self {
	    LinearizationModification::Remove { pos } => {
		let mon_idx = *linearization.get(*pos)
		    .expect("Invalid position in linearization modification");
		if monomer_graph.can_remove(mon_idx, monomers_db) {
		    Some(GraphModification::Remove { monomer_idx: mon_idx })
		}
		else { None }
	    },
	    LinearizationModification::Insert { pos, new_name } => {
		// check that there exists a possible insertion at this position with this new monomer name
		let inserts_per_position = get_inserts_for_linearization(
		    monomer_graph,
		    linearization,
		    monomers_db,
		    usize::MAX
		);
		inserts_per_position
		    .get(*pos)
		    .expect(&format!("Invalid position in linearization modification: {} in linearization of length {}", pos, linearization.len()))
		    .iter()
		    .find(|m| {
			if let GraphModification::Insert { site: _, mon_db_entry } = m {
			    mon_db_entry.monomer.features.name == *new_name
			} else {
			    false
			}
		    })
		    .cloned()
	    },
	    LinearizationModification::Substitute { pos, new_name } => {
		let mon_idx = *linearization.get(*pos)
		    .expect("Invalid position in linearization modification");	
		let maybe_entry: Option<&MonomersDB_Entry> = {
		    monomer_graph
			.possible_substitutions(mon_idx, monomers_db, true)
			.iter().cloned()
			.find(|entry| entry.monomer.features.name == *new_name)
		};
		maybe_entry.map(|entry| GraphModification::Substitute {
		    monomer_idx: mon_idx,
		    mon_db_entry: entry,
		})
	    }
	}
    }
}

#[derive(Debug, Clone, Copy)]
pub enum LinearizationComparisonMode {
    Free, // any modification allowed
    NewMonomerNamesExist, // new monomers' names should be supported (i.e. exist in in the monomers_db)
    EditsCanBePerformed, // suggested edits should be applicable to the monomer graph
}

pub fn linearization_edits<'a>(
    original_graph: &MonomerGraph,
    original_linearization: &[MonomerIdx],
    target_graph: &MonomerGraph,
    target_linearization: &[MonomerIdx],
    monomers_db: &'a MonomersDB,
    mode: LinearizationComparisonMode,
) -> Option<Vec<LinearizationModification>> {
    let mon_db_names: HashSet<NorineMonomerName> = {
	monomers_db
	    .values()
	    .flat_map(|entries| {
		entries.iter().map(|entry| entry.monomer.features.name.clone())
	    })
	    .collect()
    };

    let original_names: Vec<NorineMonomerName> = {
	original_linearization
	    .iter()
	    .map(|mon_idx| {
		original_graph.monomers
		    .get(mon_idx).unwrap()
		    .features.name
		    .clone()
	    })
	    .collect()
    };

    let target_names: Vec<NorineMonomerName> = {
	target_linearization
	    .iter()
	    .map(|mon_idx| {
		target_graph.monomers
		    .get(mon_idx).unwrap()
		    .features.name
		    .clone()
	    })
	    .collect()
    };
    
    let lin_mod_is_allowed = |lin_mod: &LinearizationModification| -> bool {
        match mode {
            LinearizationComparisonMode::Free => true,

            LinearizationComparisonMode::NewMonomerNamesExist => {
		match lin_mod {
                    LinearizationModification::Remove { .. } => true,
                    (LinearizationModification::Insert { new_name, .. }
                    | LinearizationModification::Substitute { new_name, .. }) => mon_db_names.contains(new_name)
		}
	    },

            LinearizationComparisonMode::EditsCanBePerformed => {
                LinearizationModification::try_to_graph_modification(
                    lin_mod,
                    original_graph,
                    original_linearization,
                    monomers_db,
                )
                .is_some()
            }
        }
    };

    let linearization_modifications: Vec<LinearizationModification> = {
        shortest_linearization_modifications(&original_names, &target_names, lin_mod_is_allowed)?
    };

    Some(linearization_modifications)
}



fn shortest_linearization_modifications<F>(
	original_names: &[NorineMonomerName],
	target_names: &[NorineMonomerName],
	lin_mod_is_allowed: F,
) -> Option<Vec<LinearizationModification>>
where
	F: Fn(&LinearizationModification) -> bool,
{
    #[derive(Debug, Clone)]
    struct Transition {
	from: (usize, usize),
	to: (usize, usize),
	weight: usize,
	modification: Option<LinearizationModification>,
    }

    fn apply_transition(
	dp: &mut [Vec<usize>],
	best_transition: &mut [Vec<Option<Transition>>],
	transition: Transition,
    ) {
	let (from_i, from_j) = transition.from;
	let (to_i, to_j) = transition.to;

	let candidate = dp[from_i][from_j].saturating_add(transition.weight);
	if candidate < dp[to_i][to_j] {
	    dp[to_i][to_j] = candidate;
	    best_transition[to_i][to_j] = Some(transition);
	}
    }

    let (m, n) = (original_names.len(), target_names.len());
    let inf: usize = usize::MAX / 4;

    let mut dp: Vec<Vec<usize>> = vec![vec![inf; n + 1]; m + 1];
    let mut best_transition: Vec<Vec<Option<Transition>>> = vec![vec![None; n + 1]; m + 1];
    dp[0][0] = 0;

    for i in 0..=m {
	for j in 0..=n {
	    let from = (i, j);

	    let transitions: Vec<Transition> = {
		let mut transitions = Vec::with_capacity(3);

		// Remove: (i, j) -> (i+1, j)
		if i < m {
		    let lin_mod = LinearizationModification::Remove { pos: i };
		    transitions.push(Transition {
			from,
			to: (i + 1, j),
			weight: 1,
			modification: Some(lin_mod),
		    });
		}

		// Insert: (i, j) -> (i, j+1)
		if j < n {
		    let lin_mod = LinearizationModification::Insert {
			pos: i,
			new_name: target_names[j].clone(),
		    };
		    transitions.push(Transition {
			from,
			to: (i, j + 1),
			weight: 1,
			modification: Some(lin_mod),
		    });
		}

		// Substitute: (i, j) -> (i+1, j+1)
		if i < m && j < n {
		    let lin_mod = LinearizationModification::Substitute {
			pos: i,
			new_name: target_names[j].clone(),
		    };
		    transitions.push(Transition {
			from,
			to: (i + 1, j + 1),
			weight: 1,
			modification: Some(lin_mod),
		    });
		}

		// Match: (i, j) -> (i+1, j+1) 
		if i < m && j < n && original_names[i] == target_names[j] {
		    transitions.push(Transition {
			from,
			to: (i + 1, j + 1),
			weight: 0,
			modification: None,
		    });
		}

		// Filter out invalid transitions
		transitions
		    .into_iter()
		    .filter(|t| {
			if let Some(m) = t.modification.as_ref() {
			    lin_mod_is_allowed(m)
			}
			else { true }
		    })
		    .collect()
	    };

	    for transition in transitions {
		apply_transition(&mut dp, &mut best_transition, transition);
	    }
	}
    }

    if dp[m][n] >= inf {
	return None;
    }

    let mut modifications_rev: Vec<LinearizationModification> = Vec::new();
    let mut at = (m, n);

    while at != (0, 0) {
	let (i, j) = at;
	let transition = best_transition[i][j].as_ref()?;

	if let Some(lin_mod) = &transition.modification {
	    modifications_rev.push(lin_mod.clone());
	}

	at = transition.from;
    }

    modifications_rev.reverse();
    Some(modifications_rev)
}
