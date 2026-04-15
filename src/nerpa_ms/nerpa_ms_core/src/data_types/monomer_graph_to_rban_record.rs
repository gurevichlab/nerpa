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
        let (mon_idx1, mon_idx2) = b.monomers;

        // Build the per-atomic-edge info from the bond template.
        let atomic_edges = b.get_atomic_edges();
    }
}
