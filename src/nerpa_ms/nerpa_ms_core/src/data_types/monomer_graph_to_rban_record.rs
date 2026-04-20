use std::collections::HashMap;

use crate::data_types::{bonds::Bond, common_types::MonomerIdx, monomer_graph::{AtomData, AtomicEdgeData, Monomer, MonomerGraph}, parsed_rban_record::{AtomId, AtomInfo, AtomicEdgeInfo, MonomerEdgeInfo, MonomerEdgeInfoSingle, MonomerInfo, Parsed_rBAN_Record}};

impl From<&Monomer> for MonomerInfo {
    fn from(monomer: &Monomer) -> Self {
		let atoms: Vec<AtomId> = monomer
			.atoms
			.iter()
			.map(|atom_data| atom_data.id.clone())
			.collect();
	MonomerInfo {
	    name: monomer.features.name.clone(),
	    nerpa_core: monomer.features.nerpa_core.clone(),
	    chirality: monomer.features.chirality,
	    methylated: monomer.features.methylated,
	    is_pks_hybrid: monomer.features.is_pks_hybrid,
	    mon_code: monomer.features.mon_code,
	    atoms,
	}
    }
}

impl From<&AtomData> for AtomInfo {
    fn from(atom_data: &AtomData) -> Self {
	AtomInfo {
	    name: atom_data.name.clone(),
	    hydrogens: atom_data.hydrogens,
	}
    }
}

impl From<&AtomicEdgeData> for AtomicEdgeInfo {
    fn from(edge_data: &AtomicEdgeData) -> Self {
	AtomicEdgeInfo {
	    arity: edge_data.arity.parse::<f64>().unwrap(),
	    bond_type: edge_data.bond_type.clone(),
	}
    }
}

impl From<&Bond> for MonomerEdgeInfo{
    fn from(bond: &Bond) -> Self {
        let (mon_idx1, mon_idx2) = bond.monomers;

        // Build the per-atomic-edge info from the bond template.
        let atomic_edges = bond.get_atomic_edges();
	let monomer_edges_single: Vec<MonomerEdgeInfoSingle> = {
	    atomic_edges
		.iter()
		.map(|atomic_edge| {
		    let (atom_id1, atom_id2) = &atomic_edge.atom_ids;
		    let mut monomer_to_atom: HashMap<MonomerIdx, AtomId> = HashMap::new();
		    monomer_to_atom.insert(mon_idx1, atom_id1.clone());
		    monomer_to_atom.insert(mon_idx2, atom_id2.clone());
		    MonomerEdgeInfoSingle {
			monomer_to_atom,
			arity: atomic_edge.arity.clone(),
			bond_type: atomic_edge.bond_type.clone(),
		    }
		})
		.collect()
	};
	let fst_edge = monomer_edges_single.first().unwrap().clone();
	MonomerEdgeInfo {
	    monomer_to_atom: fst_edge.monomer_to_atom,
	    arity: fst_edge.arity,
	    bond_type: fst_edge.bond_type,
	    all_edges: monomer_edges_single,
	}
    }
}

impl From<&MonomerGraph> for Parsed_rBAN_Record {
    fn from(monomer_graph: &MonomerGraph) -> Self {
	let monomers: HashMap<MonomerIdx, MonomerInfo> = {
	    monomer_graph
		.monomers
		.iter()
		.map(|(idx, monomer)| (*idx, MonomerInfo::from(monomer)))
		.collect()
	};

	let monomer_bonds: HashMap<(MonomerIdx, MonomerIdx), MonomerEdgeInfo> = {
	    monomer_graph
	    .monomer_bonds
	    .iter()
	    .map(|bond| (bond.monomers.clone(), MonomerEdgeInfo::from(bond)))
	    .collect()
	};

	let atoms: HashMap<AtomId, AtomInfo> = {
	    monomer_graph
	    .monomers
	    .values()
	    .flat_map(|monomer| &monomer.atoms)
	    .map(|atom_data| (atom_data.id.clone(), AtomInfo::from(atom_data)))
	    .collect()
	};

	let atomic_bonds: HashMap<(AtomId, AtomId), AtomicEdgeInfo> = {
	    let mut atomic_bonds_vec: Vec<((AtomId, AtomId), AtomicEdgeInfo)> = Vec::new();
	    // Add intra-monomer atomic bonds
	    for monomer in monomer_graph.monomers.values() {
		for atomic_bond in &monomer.atomic_bonds {
		    let edge_info = AtomicEdgeInfo::from(atomic_bond);
		    atomic_bonds_vec.push((atomic_bond.atom_ids.clone(), edge_info));
		}
	    }
	    // Add inter-monomer atomic bonds from monomer_bonds
	    for bond in &monomer_graph.monomer_bonds {
		let atomic_edges = bond.get_atomic_edges();
		for atomic_edge in atomic_edges {
		    let edge_info = AtomicEdgeInfo::from(&atomic_edge);
		    atomic_bonds_vec.push((atomic_edge.atom_ids.clone(), edge_info));
		}
	    }
	    atomic_bonds_vec.into_iter().collect()
	};

	Parsed_rBAN_Record {
	    monomers,
	    monomer_bonds,
	    atoms,
	    atomic_bonds,
	    compound_id: monomer_graph.compound_id.clone(),
	    metadata: monomer_graph.metadata.clone(),
	}
    }
}
