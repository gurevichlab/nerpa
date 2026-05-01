use crate::data_types::bond_consts::{AMINO_C_END_PROFILE, AMINO_MIDDLE_PROFILE, AMINO_N_END_PROFILE};
use crate::data_types::monomer_graph::{Monomer, MonomerGraph};
use crate::data_types::monomers_db::{MonomersDB, MonomersDB_Entry, get_entry_by_profile_and_name};
use crate::data_types::bonds::{
    AMINO_BOND_TEMPLATE, BindingSiteType, BindingSitesProfile, Bond, BondSide, BondsByBSType
};
use crate::data_types::common_types::MonomerIdx
use crate::data_types::parsed_rban_record::NorineMonomerName;
use std::collections::HashMap;

pub fn get_entries_for_insertion_between<'a>(
	monomers_db: &'a MonomersDB,
) -> Vec<&'a MonomersDB_Entry> {

    monomers_db
        .get(&*AMINO_MIDDLE_PROFILE)
        .unwrap_or_else(|| panic!("Expected the monomers DB to contain entries with the amino middle profile"))
	.iter()
	.collect()
}

pub fn splice_bonds(mon_with_bond1: (MonomerIdx, &Bond),
		    mon_with_bond2: (MonomerIdx, &Bond)) -> Bond {
    // Takes two monomers with their respective bonds, where the two bonds have the same template and the monomers appear on opposite sides in those bonds, and returns a new bond that "splices" them.
    // E.g. (m1, C1)->(m2, N2) and (m3, C3)->(m4, N4) yields (m1, C1)->(m4, N4)
    let (mon1, bond1) = mon_with_bond1;
    let (mon2, bond2) = mon_with_bond2;
    debug_assert_eq!(bond1.bond_templ, bond2.bond_templ);

    let side1 = if bond1.monomers.0 == mon1 {
	BondSide::Left
    } else if bond1.monomers.1 == mon1 {
	BondSide::Right
    } else {
	panic!("mon_with_bond1 does not contain the monomer that is part of the bond");
    };

    let side2 = if bond2.monomers.0 == mon2 {
	BondSide::Left
    } else if bond2.monomers.1 == mon2 {
	BondSide::Right
    } else {
	panic!("mon_with_bond2 does not contain the monomer that is part of the bond");
    };

    let bond_templ = bond1.bond_templ.clone();
    let monomers;
    let label_to_atom;
    match (side1, side2) {
	(BondSide::Left, BondSide::Right) => {
	    monomers = (
		bond1.monomers.0,
		bond2.monomers.1,
	    );

	    label_to_atom = (
		bond1.label_to_atom.0.clone(),
		bond2.label_to_atom.1.clone(),
	    );
	},
	(BondSide::Right, BondSide::Left) => {
	    monomers = (
		bond2.monomers.0,
		bond1.monomers.1,
	    );

	    label_to_atom = (
		bond2.label_to_atom.0.clone(),
		bond1.label_to_atom.1.clone(),
	    );
	},
	_ => panic!("The two bonds should have opposite sides for the monomers being connected"),
    };

    Bond {
	bond_templ,
	monomers,
	label_to_atom,
    }
}

fn _profile_without_bs_type(
    profile: BindingSitesProfile,
    bs_type_to_remove: &BindingSiteType,
) -> BindingSitesProfile {
    let mut v = profile.into_vec();
    let idx_to_remove = v
        .iter()
        .position(|bs| bs == bs_type_to_remove)
        .expect("Expected the binding site type to exist in the profile");
    v.remove(idx_to_remove);
    BindingSitesProfile::new(v)
}

struct _DataForRemovalDeg1<'a> {
    parent_idx: MonomerIdx,
    parent_sub: &'a MonomersDB_Entry,
}

impl MonomerGraph {
    pub fn substitute(&mut self, monomer_idx: MonomerIdx, mon_db_entry: &MonomersDB_Entry) {
	let mon_bonds_by_bs = self.bonds_by_bs_type(monomer_idx);
	if !mon_bonds_by_bs.compatible_with(&mon_db_entry.bonds_by_bs) {
	    panic!("Attempting to substitute monomer {} with a monomer that has an incompatible binding sites profile", monomer_idx);
	}

	// Shifting the atom IDs in the db entry to avoid collisions with existing atom IDs in the graph
	let corrected_db_entry = {
	    let max_id: u32 = self.monomers.values()
		.flat_map(|mon| mon.atoms.iter().map(|atom| atom.id.0))
		.max()
		.unwrap_or(0);
	    let mut db_entry = mon_db_entry.clone();
	    db_entry.shift_atom_ids(max_id + 1);
	    db_entry.set_monomer_idx(monomer_idx);
	    db_entry
	};

	let new_bonds_by_mon_indices = {
	    let mut new_bonds_by_mons = HashMap::new();
	    for i in 0..mon_bonds_by_bs.len() {
		let (mon_bs, mon_bond) =
		    mon_bonds_by_bs.get(i).unwrap();
		let (db_entry_bs, db_entry_bond) =
		    corrected_db_entry.bonds_by_bs.get(i).unwrap();
		debug_assert_eq!(*mon_bs, *db_entry_bs);

		let other_mon_idx = if mon_bond.monomers.0 == monomer_idx {
		    mon_bond.monomers.1
		} else {
		    mon_bond.monomers.0
		};
		let new_bond = splice_bonds((monomer_idx, db_entry_bond),
					    (other_mon_idx, mon_bond));

		new_bonds_by_mons.insert(new_bond.monomers, new_bond);
	    }
	    new_bonds_by_mons
	};

	for old_bond in self.monomer_bonds.iter_mut() {
	    if let Some(new_bond) = new_bonds_by_mon_indices.get(&old_bond.monomers) {
		*old_bond = new_bond.clone();
	    }
	}

	self.monomers.insert(monomer_idx, corrected_db_entry.monomer.clone());

    }

    fn _leaf_parent_with_bs_type(
        &self,
        monomer_idx: MonomerIdx,
    ) -> (MonomerIdx, BindingSiteType) {
	let bonds_by_bs = self.bonds_by_bs_type(monomer_idx);
	debug_assert_eq!(bonds_by_bs.len(), 1);
	let bs_and_bond = bonds_by_bs.get(0).unwrap();
	
        let bond = &bs_and_bond.1;

        let other_mon = if bond.monomers.0 == monomer_idx {
            bond.monomers.1
        } else {
            bond.monomers.0
        };

        let other_bs_type = BindingSiteType {
            bond_templ: bond.bond_templ.clone(),
            side: match bs_and_bond.0.side {
                BondSide::Left => BondSide::Right,
                BondSide::Right => BondSide::Left,
            },
        };

        (other_mon, other_bs_type)
    }

    fn _data_to_remove_deg1<'a>(&self,
				monomer_idx: MonomerIdx,
				monomers_db: &'a MonomersDB) -> Option<_DataForRemovalDeg1<'a>> {
	debug_assert_eq!(self.bonds_by_bs_type(monomer_idx).len(), 1);
	
        let (other_mon, other_bs_type) =
            self._leaf_parent_with_bs_type(monomer_idx);

        let other_bs_profile = self.bonds_by_bs_type(other_mon).get_profile();
        let other_bs_profile_wo_this_bs =
            _profile_without_bs_type(other_bs_profile, &other_bs_type);

        let other_name = &self.monomers.get(&other_mon).unwrap().features.name;

        let maybe_entry = get_entry_by_profile_and_name(monomers_db,
							&other_bs_profile_wo_this_bs,
							other_name);
	match maybe_entry {
	    Some(entry) => Some(_DataForRemovalDeg1 {
		parent_idx: other_mon,
		parent_sub: entry,
	    }),
	    None => None,
	}
    }

    fn _can_remove_deg1(&self,
			monomer_idx: MonomerIdx,
			monomers_db: &MonomersDB) -> bool {
	self._data_to_remove_deg1(monomer_idx, monomers_db).is_some()
    }

    fn _remove_deg1(&mut self,
		    monomer_idx: MonomerIdx,
		    monomers_db: &MonomersDB) {
	debug_assert!(self._can_remove_deg1(monomer_idx, monomers_db));
	let data_for_removal = self._data_to_remove_deg1(monomer_idx, monomers_db)
	    .expect("Expected data for removal to be present since _can_remove_deg1 returned true");

        // 1) Remove the monomer and its only incident bond.
        self.monomer_bonds
            .retain(|b| b.monomers.0 != monomer_idx && b.monomers.1 != monomer_idx);
        self.monomers.remove(&monomer_idx);

        // 2) Now the other monomer's binding-sites profile has changed,
        // so we can substitute it with the corresponding degree-reduced DB entry.
        self.substitute(data_for_removal.parent_idx, data_for_removal.parent_sub);
    }

    fn _can_remove_deg2(&self,
			monomer_idx: MonomerIdx) -> bool {
	let bonds_by_bs = self.bonds_by_bs_type(monomer_idx);
	debug_assert_eq!(bonds_by_bs.len(), 2);
	let fst_bs = &bonds_by_bs.get(0).unwrap().0;
	let snd_bs = &bonds_by_bs.get(1).unwrap().0;
	fst_bs.bond_templ == snd_bs.bond_templ
	    && fst_bs.side != snd_bs.side
    }

    fn _remove_deg2(&mut self,
		    monomer_idx: MonomerIdx) {
	debug_assert!(self._can_remove_deg2(monomer_idx));

        // We are removing a degree-2 monomer which has two bonds with the same template,
        // and the monomer appears on opposite sides (Left/Right) in those two bonds.
        //
        // So we "splice" it out:
        //   (mon_l) --[templ]-- (monomer_idx) --[templ]-- (mon_r)
        // becomes
        //   (mon_l) --[templ]-- (mon_r)

	let bonds_by_bs = self.bonds_by_bs_type(monomer_idx);
	debug_assert_eq!(bonds_by_bs.len(), 2);
        let (_bs_a, bond_a) = bonds_by_bs.get(0).unwrap();
        let (_bs_b, bond_b) = bonds_by_bs.get(1).unwrap();
	let mon_a = if bond_a.monomers.0 == monomer_idx {
	    bond_a.monomers.1
	} else {
	    bond_a.monomers.0
	};
	let mon_b = if bond_b.monomers.0 == monomer_idx {
	    bond_b.monomers.1
	} else {
	    bond_b.monomers.0
	};

	let new_bond = splice_bonds((mon_a, bond_a), (mon_b, bond_b));
	
	// Remove the two old bonds and add the new spliced bond.
	self.monomer_bonds
	    .retain(|b| b.monomers.0 != monomer_idx && b.monomers.1 != monomer_idx);
	self.monomer_bonds.push(new_bond);

	// Finally, remove the monomer itself.
	self.monomers.remove(&monomer_idx);
    }


    pub fn can_remove(&mut self,
		      monomer_idx: MonomerIdx,
		      monomers_db: &MonomersDB) -> bool {
	let bonds_by_bs = self.bonds_by_bs_type(monomer_idx);
	match bonds_by_bs.len() {
	    1 => self._can_remove_deg1(monomer_idx, monomers_db),
	    2 => self._can_remove_deg2(monomer_idx),
	    _ => false, // monomers with degree > 2 cannot be removed
	}
    }

    pub fn remove(&mut self,
		  monomer_idx: MonomerIdx,
    		  monomers_db: &MonomersDB) {
	if !self.can_remove(monomer_idx, monomers_db) {
	    panic!("Attempting to remove monomer {} that cannot be removed", monomer_idx);
	}
	let bonds_by_bs = self.bonds_by_bs_type(monomer_idx);
	match bonds_by_bs.len() {
	    1 => {
		self._remove_deg1(monomer_idx, monomers_db);
	    },
	    2 => {
		self._remove_deg2(monomer_idx);
	    },
	    _ => unreachable!(), // monomers with degree > 2 cannot be removed, so this case should have been ruled out by can_remove
	}
    }


   // -------------------------
   // Insertion helpers
   // -------------------------

   pub fn possible_insertions_between<'a>(
       &self,
       mon1: MonomerIdx,
       mon2: MonomerIdx,
       monomers_db: &'a MonomersDB,
   ) -> Vec<&'a MonomersDB_Entry> {
       let bond = self.get_bond(mon1, mon2)
           .unwrap_or_else(|| panic!("bond between monomers {:?} not found",
				     (mon1, mon2)));

       if &bond.bond_templ != &*AMINO_BOND_TEMPLATE {
           return Vec::new();
       }

       get_entries_for_insertion_between(monomers_db)
   }

    fn prepare_entry_for_insertion(&self,
				   mon_db_entry: &MonomersDB_Entry) -> (MonomerIdx, MonomersDB_Entry) {
	let new_idx = {
	    let max_idx = self.monomers.keys()
		.max()
		.unwrap();
	    MonomerIdx(max_idx.0 + 1)
	};
	let max_atom_id = self.monomers.values()
	    .flat_map(|mon| mon.atoms.iter().map(|atom| atom.id.0))
	    .max()
	    .unwrap_or(0);

	let corrected_db_entry = {
	    let mut new_entry = mon_db_entry.clone();
	    new_entry.set_monomer_idx(new_idx);
	    new_entry.shift_atom_ids(max_atom_id + 1);
	    new_entry
	};

	(new_idx, corrected_db_entry)
    }
	    


   pub fn insert_between(
       &mut self,
       mon1: MonomerIdx,
       mon2: MonomerIdx,
       mon_db_entry: &MonomersDB_Entry,
   ) {
       let (new_idx, corrected_db_entry) = self.prepare_entry_for_insertion(mon_db_entry);

       let old_bond = self.get_bond(mon1, mon2)
	   .unwrap_or_else(|| panic!("bond between monomers {:?} not found",
				     (mon1, mon2)));


       // Normalize orientation: old bond is (left_mon) -> (right_mon)
       let left_mon = old_bond.monomers.0;
       let right_mon = old_bond.monomers.1;
       
       // Database monomer should have exactly two bonds
       // with the same template as the bond we are inserting into
       // The bond mon_db_bond_left should have mon_db_entry on the Left, and mon_db_bond_right should have it on the Right

       let (mon_db_bond_left, mon_db_bond_right) = {
	   debug_assert_eq!(corrected_db_entry.bonds_by_bs.len(), 2);
	   let (bs_a, bond_a) = &corrected_db_entry.bonds_by_bs.get(0).unwrap();
	   let (bs_b, bond_b) = &corrected_db_entry.bonds_by_bs.get(1).unwrap();
	   debug_assert_eq!(bs_a.bond_templ, bs_b.bond_templ);
	   debug_assert_eq!(bs_a.bond_templ, old_bond.bond_templ);
	   match (bs_a.side, bs_b.side) {
	       (BondSide::Left, BondSide::Right) => (bond_a, bond_b),
	       (BondSide::Right, BondSide::Left) => (bond_b, bond_a),
	       _ => panic!("DB entry bonds should have opposite sides"),
	   }
       };

       let left_mon_new_bond = splice_bonds((left_mon, old_bond),
					    (new_idx, mon_db_bond_right));
       let right_mon_new_bond = splice_bonds((right_mon, old_bond),
					     (new_idx, mon_db_bond_left));

       self.monomer_bonds.retain(|b| {
           !((b.monomers.0 == left_mon && b.monomers.1 == right_mon)
             || (b.monomers.0 == right_mon && b.monomers.1 == left_mon))
       });
       self.monomer_bonds.push(left_mon_new_bond);
       self.monomer_bonds.push(right_mon_new_bond);

       // Insert the new monomer itself.
       self.monomers
           .insert(new_idx, corrected_db_entry.monomer.clone());
   }

    fn leaf_substitute_for_attaching<'a>(&self,
					 monomer_idx: MonomerIdx,
					 monomers_db: &'a MonomersDB) -> Option<&'a MonomersDB_Entry> {
	// Check if the given leaf monomer can be substituted with a DB entry
	// in a way that would allow attaching a new leaf monomer to it.
	let bs_profile = self.bonds_by_bs_type(monomer_idx).get_profile();
	if &bs_profile != &*AMINO_C_END_PROFILE
	    && &bs_profile != &*AMINO_N_END_PROFILE {
	    return None;
	}

	let mon = self.monomers.get(&monomer_idx).unwrap();
	get_entry_by_profile_and_name(monomers_db,
				      &*AMINO_MIDDLE_PROFILE,
				      &mon.features.name)
    }

   pub fn possible_insertions_at_leaf<'a>(
       &self,
       monomer_idx: MonomerIdx,
       monomers_db: &'a MonomersDB,
   ) -> Vec<&'a MonomersDB_Entry> {
       if self.leaf_substitute_for_attaching(monomer_idx, monomers_db).is_none() {
	   return Vec::new();
       }
       let bs_profile = self.bonds_by_bs_type(monomer_idx).get_profile();

       if &bs_profile == &*AMINO_C_END_PROFILE {
	   monomers_db
	       .get(&*AMINO_C_END_PROFILE)
	       .unwrap_or_else(|| panic!("Expected the monomers DB to contain entries with the amino middle profile"))
	       .iter()
	       .collect()
       }
       else if &bs_profile == &*AMINO_N_END_PROFILE {
	   monomers_db
	       .get(&*AMINO_N_END_PROFILE)
	       .unwrap_or_else(|| panic!("Expected the monomers DB to contain entries with the amino middle profile"))
	       .iter()
	       .collect()
	   }
       else {
	   unreachable!("The monomer has either the amino C-end or N-end profile since _data_for_attaching_to_leaf returned Some");
       }

   }

    pub fn attach_to_leaf(&mut self,
			  monomer_idx: MonomerIdx,
			  mon_db_entry: &MonomersDB_Entry,
			  monomers_db: &MonomersDB) {
	// Update the monomer for substituting the current leaf to avoid index collisions
	let old_leaf_sub_entry = {
	    let mut sub_entry = self
		.leaf_substitute_for_attaching(monomer_idx, monomers_db)
		.expect("Expected the leaf monomer to be substitutable for attaching
since can_insert_at_leaf should have been called before and returned true")
		.clone();

	    let max_atom_id = self.monomers.values()
		.flat_map(|mon| mon.atoms.iter().map(|atom| atom.id.0))
		.max()
		.unwrap();
	    sub_entry.set_monomer_idx(monomer_idx);
	    sub_entry.shift_atom_ids(max_atom_id + 1);
	    sub_entry
	};
	debug_assert_eq!(&old_leaf_sub_entry.bonds_by_bs.get_profile(), &*AMINO_MIDDLE_PROFILE);

	let old_leaf_bonds_by_bs = self.bonds_by_bs_type(monomer_idx);
	debug_assert_eq!(old_leaf_bonds_by_bs.len(), 1);
	let (old_leaf_bs, old_leaf_bond) = old_leaf_bonds_by_bs.get(0).unwrap();

	// identify which of the two bonds in the DB entry is the one that will be attached to the old leaf's parent, and which one will be the dangling bond to which the new leaf will be attached
	let old_leaf_sub_attached_bond;
	let old_leaf_sub_dangling_bond;
	{
	    let old_leaf_sub_bonds_by_bs = &old_leaf_sub_entry.bonds_by_bs.as_slice();
	    let (old_leaf_sub_bs_1, old_leaf_sub_bond_1) = &old_leaf_sub_bonds_by_bs[0];
	    let (old_leaf_sub_bs_2, old_leaf_sub_bond_2) = &old_leaf_sub_bonds_by_bs[1];
	    if old_leaf_sub_bs_1 == old_leaf_bs {
		old_leaf_sub_attached_bond = old_leaf_sub_bond_1;
		old_leaf_sub_dangling_bond = old_leaf_sub_bond_2;
	    }
	    else {
		debug_assert_eq!(old_leaf_sub_bs_2, old_leaf_bs);
		old_leaf_sub_attached_bond = old_leaf_sub_bond_2;
		old_leaf_sub_dangling_bond = old_leaf_sub_bond_1;
	    }
	}

	// substitute the old leaf monomer with the DB entry
	{
	    let parent_idx = if old_leaf_bond.monomers.0 == monomer_idx {
		old_leaf_bond.monomers.1
	    } else {
		old_leaf_bond.monomers.0
	    };
	    let new_parent_bond = splice_bonds((monomer_idx, old_leaf_sub_attached_bond),
					       (parent_idx, old_leaf_bond));
	    self.monomer_bonds
		.retain(|b| b.monomers.0 != monomer_idx && b.monomers.1 != monomer_idx);
	    self.monomer_bonds.push(new_parent_bond);
	    self.monomers[&monomer_idx] = old_leaf_sub_entry.monomer;
	}
	
	// attach the new leaf
	let (new_idx, corrected_db_entry) = self.prepare_entry_for_insertion(mon_db_entry);
	let corrected_db_entry_bond = &corrected_db_entry.bonds_by_bs.get(0).unwrap().1;
	let new_leaf_bond = splice_bonds((new_idx, corrected_db_entry_bond),
					 (monomer_idx, old_leaf_sub_dangling_bond));

	self.monomer_bonds.push(new_leaf_bond);
	self.monomers.insert(new_idx, corrected_db_entry.monomer);
    }
}
