use std::collections::HashMap;

impl From<&MonomerInfo> for MonomerFeatures {
	fn from(monomer_info: &MonomerInfo) -> Self {
		MonomerFeatures {
			name: monomer_info.name.clone(),
			nerpa_core: monomer_info.nerpa_core.clone(),
			chirality: monomer_info.chirality,
			methylated: monomer_info.methylated,
			is_pks_hybrid: monomer_info.is_pks_hybrid,
			mon_code: monomer_info.mon_code,
		}
	}
}

impl Monomer {
    pub fn from_rban_record(
	rban_record: &Parsed_rBAN_Record,
	idx: MonomerIdx,
    ) -> Self {
	let monomer_info = rban_record.monomers.get(&idx)
	    .expect("Missing MonomerInfo for monomer index {:?}", idx);
	let atoms_with_data = monomer_info
	    .atoms
	    .iter()
	    .map(|atom_id| {
		let atom_info = rban_record.atoms.get(atom_id).unwrap();
		AtomData {
		    id: atom_id.clone(),
		    name: atom_info.name.clone(),
		    hydrogens: atom_info.hydrogens,
		}
	    })
	    .collect();
	let atoms_set: HashSet<AtomId> = monomer_info
	    .atoms
	    .iter()
	    .cloned()
	    .collect();
	let atomic_bonds = rban_record.atomic_bonds
	    .iter()
	    .filter_map(|((a1, a2), edge_info)| {
		if atoms_set.contains(a1) && atoms_set.contains(a2) {
		    Some(AtomicEdgeData {
			atom_ids: (a1.clone(), a2.clone()),
			arity: format!("{}", edge_info.arity),
			bond_type: edge_info.bond_type.clone(),
		    })
		}
		else {
		    None
		}
	    })
	    .collect();

	Monomer {
	    atoms: atoms_with_data,
	    atomic_bonds,
	    features: MonomerFeatures::from(monomer_info),
	}
    }
}

impl From<&Parsed_rBAN_Record> for MonomerGraph {
    fn from(rban_record: &Parsed_rBAN_Record) -> Self {
        let monomers: HashMap<MonomerIdx, Monomer> = {
	    rban_record
	    .monomers
	    .iter()
	    .map(|(&idx, monomer_info)| {
		let monomer = Monomer::from_rban_record(rban_record, idx);
		(idx, monomer)
	    })
	    .collect()
	};

        // Deterministic order is handy for tests/debugging.
        let edges: Vec<(MonomerIdx, MonomerIdx)> = {
	    let mut edges: Vec<(MonomerIdx, MonomerIdx)> =
		rban_record.monomer_bonds.keys().cloned().collect();
            edges.sort();
	    edges
	};

        let monomer_bonds = edges
            .iter()
            .map(|edge| rban_record.edge_to_bond(edge))
            .collect();

        MonomerGraph {
            monomers,
            monomer_bonds,
            compound_id: rban_record.compound_id.clone(),
            metadata: rban_record.metadata.clone(),
        }
    }
}
