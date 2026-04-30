use crate::data_types::bond_consts::AMINO_MIDDLE_PROFILE;
use crate::data_types::monomer_graph::{Monomer, MonomerGraph};
use crate::data_types::monomers_db::{MonomersDB, MonomersDB_Entry};
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

impl MonomerGraph {
    pub fn substitute(&mut self, monomer_idx: MonomerIdx, mon_db_entry: &MonomersDB_Entry) {
	let mon_bonds_by_bs = self.bonds_by_bs_type(monomer_idx);
	if !mon_bonds_by_bs.compatible_with(&mon_db_entry.bonds_by_bs) {
	    panic!("Attempting to substitute monomer {} with a monomer that has an incompatible binding sites profile", monomer_idx);
	}

	// Shifting the atom IDs in the db entry to avoid collisions with existing atom IDs in the graph
	let db_entry_shifted = {
	    let max_id: u32 = self.monomers.values()
		.flat_map(|mon| mon.atoms.iter().map(|atom| atom.id.0))
		.max()
		.unwrap_or(0);
	    mon_db_entry.shift_atom_ids(max_id + 1)
	};

	let new_bonds_by_mon_indices = {
	    let mut new_bonds_by_mons = HashMap::new();
	    for i in 0..mon_bonds_by_bs.len() {
		let (mon_bs, mon_bond) =
		    mon_bonds_by_bs.get(i).unwrap();
		let (db_entry_bs, db_entry_bond) =
		    db_entry_shifted.bonds_by_bs.get(i).unwrap();
		debug_assert_eq!(*mon_bs, *db_entry_bs);
		let side = mon_bs.side;

		let mut new_bond = mon_bond.clone();
		match side {
		    BondSide::Left => {
			new_bond.label_to_atom.0 = db_entry_bond.label_to_atom.0.clone();
		    },
		    BondSide::Right => {
			new_bond.label_to_atom.1 = db_entry_bond.label_to_atom.1.clone();
		    }
		}

		new_bonds_by_mons.insert(mon_bond.monomers, new_bond);
	    }
	    new_bonds_by_mons
	};

	for old_bond in self.monomer_bonds.iter_mut() {
	    if let Some(new_bond) = new_bonds_by_mon_indices.get(&old_bond.monomers) {
		*old_bond = new_bond.clone();
	    }
	}

	self.monomers.insert(monomer_idx, db_entry_shifted.monomer.clone());

    }

    // --- helpers used by both _can_remove_deg1 and _remove_deg1 ---

    fn _deg1_other_monomer_and_bs_type(
        &self,
        monomer_idx: MonomerIdx,
        bs_and_bond: &(BindingSiteType, Bond),
    ) -> (MonomerIdx, BindingSiteType) {
        let bond = &bs_and_bond.1;

        let other_mon = if bond.monomers.0 == monomer_idx {
            bond.monomers.1
        } else if bond.monomers.1 == monomer_idx {
            bond.monomers.0
        } else {
            panic!(
                "Bond {:?} is not incident to monomer {}",
                bs_and_bond.1, monomer_idx
            );
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


    fn _can_remove_deg1(&self,
			monomer_idx: MonomerIdx,
			bs_and_bond: &(BindingSiteType, Bond),
			monomers_db: &MonomersDB) -> bool {
        let (other_mon, other_bs_type) =
            self._deg1_other_monomer_and_bs_type(monomer_idx, bs_and_bond);

        let other_bs_profile = self.bonds_by_bs_type(other_mon).get_profile();
        let other_bs_profile_wo_this_bs =
            _profile_without_bs_type(other_bs_profile, &other_bs_type);

        let other_name = &self.monomers.get(&other_mon).unwrap().features.name;

        let maybe_entry = _find_db_entry_by_profile_and_name(monomers_db,
					   &other_bs_profile_wo_this_bs,
					   other_name);
	maybe_entry.is_some()
    }

    fn _remove_deg1(&mut self,
		    monomer_idx: MonomerIdx,
		    bs_and_bond: &(BindingSiteType, Bond),
		    monomers_db: &MonomersDB) {
	debug_assert!(self._can_remove_deg1(monomer_idx, bs_and_bond, monomers_db));
        let (other_mon, other_bs_type) =
            self._deg1_other_monomer_and_bs_type(monomer_idx, bs_and_bond);

        let other_bs_profile = self.bonds_by_bs_type(other_mon).get_profile();
        let other_bs_profile_wo_this_bs =
            _profile_without_bs_type(other_bs_profile, &other_bs_type);

        let other_name = self.monomers.get(&other_mon).unwrap().features.name.clone();
        let db_entry = _find_db_entry_by_profile_and_name(
            monomers_db,
            &other_bs_profile_wo_this_bs,
            &other_name,
        )
        .expect("Expected a matching DB entry (should have been checked by _can_remove_deg1)");

        // 1) Remove the monomer and its only incident bond.
        self.monomer_bonds
            .retain(|b| b.monomers.0 != monomer_idx && b.monomers.1 != monomer_idx);
        self.monomers.remove(&monomer_idx);

        // 2) Now the other monomer's binding-sites profile has changed,
        // so we can substitute it with the corresponding degree-reduced DB entry.
        self.substitute(other_mon, db_entry);
    }

    fn _can_remove_deg2(&self,
			bonds_by_bs: &BondsByBSType) -> bool {
	debug_assert_eq!(bonds_by_bs.len(), 2);
	let fst_bs = &bonds_by_bs.get(0).unwrap().0;
	let snd_bs = &bonds_by_bs.get(1).unwrap().0;
	fst_bs.bond_templ == snd_bs.bond_templ
	    && fst_bs.side != snd_bs.side
    }

    fn _remove_deg2(&mut self,
		    monomer_idx: MonomerIdx,
		    bonds_by_bs: &BondsByBSType) {
	debug_assert!(self._can_remove_deg2(monomer_idx, bonds_by_bs));

        // We are removing a degree-2 monomer which has two bonds with the same template,
        // and the monomer appears on opposite sides (Left/Right) in those two bonds.
        //
        // So we "splice" it out:
        //   (mon_l) --[templ]-- (monomer_idx) --[templ]-- (mon_r)
        // becomes
        //   (mon_l) --[templ]-- (mon_r)

        let pairs = bonds_by_bs.as_slice();
        debug_assert_eq!(pairs.len(), 2);

        let (bs_a, bond_a) = &pairs[0];
        let (bs_b, bond_b) = &pairs[1];

        // Identify which bond has `monomer_idx` on the Left vs Right.
        let (bond_with_monomer_on_left, bond_with_monomer_on_right) = match (bs_a.side, bs_b.side)
        {
            (BondSide::Left, BondSide::Right) => (bond_a, bond_b),
            (BondSide::Right, BondSide::Left) => (bond_b, bond_a),
            _ => panic!(
                "_remove_deg2 called on a monomer that is not (Left, Right) for the two bonds"
            ),
        };

        // Sanity checks: ensure the tuple ordering matches the claimed side.
        if bond_with_monomer_on_left.monomers.0 != monomer_idx {
            panic!("Expected monomer {} to be on the Left of its Left-side bond", monomer_idx);
        }
        if bond_with_monomer_on_right.monomers.1 != monomer_idx {
            panic!("Expected monomer {} to be on the Right of its Right-side bond", monomer_idx);
        }

        let monomer_left_neighbor = bond_with_monomer_on_right.monomers.0;
        let monomer_right_neighbor = bond_with_monomer_on_left.monomers.1;

        if monomer_left_neighbor == monomer_right_neighbor {
            panic!(
                "Refusing to splice monomer {}: would create a self-loop on monomer {}",
                monomer_idx, monomer_left_neighbor
            );
        }

        // Avoid accidentally creating a duplicate/parallel edge unless you explicitly want that.
        let would_duplicate = self.monomer_bonds.iter().any(|b| {
            (b.monomers == (monomer_left_neighbor, monomer_right_neighbor)
                || b.monomers == (monomer_right_neighbor, monomer_left_neighbor))
        });
        if would_duplicate {
            panic!(
                "Refusing to splice monomer {}: a bond already exists between {:?}",
                monomer_idx,
                (monomer_left_neighbor, monomer_right_neighbor)
            );
        }

        // Build the new bond by taking the label->atom mapping from the *neighbors*:
        // - left neighbor mapping comes from the bond where it is on the Left (bond_with_monomer_on_right)
        // - right neighbor mapping comes from the bond where it is on the Right (bond_with_monomer_on_left)
        let new_bond = Bond {
            bond_templ: bond_with_monomer_on_left.bond_templ.clone(),
            monomers: (monomer_left_neighbor, monomer_right_neighbor),
            label_to_atom: (
                bond_with_monomer_on_right.label_to_atom.0.clone(),
                bond_with_monomer_on_left.label_to_atom.1.clone(),
            ),
        };

        // Remove all bonds incident to the removed monomer (degree is 2, so this removes exactly two).
        self.monomer_bonds
            .retain(|b| b.monomers.0 != monomer_idx && b.monomers.1 != monomer_idx);

        // Insert the new spliced bond.
        self.monomer_bonds.push(new_bond);

        // Finally remove the monomer itself from the graph.
        self.monomers.remove(&monomer_idx);

    }


    pub fn can_remove(&mut self,
		      monomer_idx: MonomerIdx,
		      monomers_db: &MonomersDB) -> bool {
	let bonds_by_bs = self.bonds_by_bs_type(monomer_idx);
	match bonds_by_bs.len() {
	    1 => self._can_remove_deg1(monomer_idx, bonds_by_bs.get(0).unwrap(), monomers_db),
	    2 => self._can_remove_deg2(&bonds_by_bs),
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
		let bs_and_bond = bonds_by_bs.get(0).unwrap();
		self._remove_deg1(monomer_idx, &bs_and_bond, monomers_db);
	    },
	    2 => {
		self._remove_deg2(monomer_idx, &bonds_by_bs);
	    },
	    _ => unreachable!(), // monomers with degree > 2 cannot be removed, so this case should have been ruled out by can_remove
	}
    }


   // -------------------------
   // Insertion helpers
   // -------------------------

   fn _find_bond_idx_between(&self, mon1: MonomerIdx, mon2: MonomerIdx) -> Option<usize> {
       self.monomer_bonds.iter().position(|b| {
           (b.monomers.0 == mon1 && b.monomers.1 == mon2)
               || (b.monomers.0 == mon2 && b.monomers.1 == mon1)
       })
   }

   /// Very small "for simplicity" check: only allow insertion into a single (arity=1) bond.
   /// You said "assume C-N amino bond" — this is the concrete check we can do locally without
   /// having a dedicated "Amino" bond kind.
   fn _is_supported_insertion_bond_templ(bond_templ: &BondTemplate) -> bool {
       let slice = bond_templ.as_slice();
       if slice.len() != 1 {
           return false;
       }
       let at = &slice[0];
       at.arity == "1" && at.bond_type == BondType::Single
   }

   fn _required_deg2_insertion_profile(bond_templ: &BondTemplate) -> BindingSitesProfile {
       BindingSitesProfile::new(vec![
           BindingSiteType {
               bond_templ: bond_templ.clone(),
               side: BondSide::Left, // "C-amino" side in your simplified model
           },
           BindingSiteType {
               bond_templ: bond_templ.clone(),
               side: BondSide::Right, // "N-amino" side in your simplified model
           },
       ])
   }

   fn _next_monomer_idx(&self) -> MonomerIdx {
       // Assumes MonomerIdx is a tuple struct like `pub struct MonomerIdx(pub u32);`
       let next = self
           .monomers
           .keys()
           .map(|idx| idx.0)
           .max()
           .unwrap_or(0)
           + 1;
       MonomerIdx(next)
   }

   fn _max_atom_id(&self) -> u32 {
       self.monomers
           .values()
           .flat_map(|mon| mon.atoms.iter().map(|atom| atom.id.0))
           .max()
           .unwrap_or(0)
   }

   fn _entry_label_map_for_side(
       entry: &MonomersDB_Entry,
       bond_templ: &BondTemplate,
       side: BondSide,
   ) -> HashMap<BondAtomLabel, AtomId> {
       let (bs_type, bond) = entry
           .bonds_by_bs
           .iter()
           .find(|(bs, _b)| bs.bond_templ == *bond_templ && bs.side == side)
           .expect("DB entry is expected to have the requested BindingSiteType (caller should ensure profile fits)");

       debug_assert_eq!(bs_type.bond_templ, *bond_templ);
       debug_assert_eq!(bs_type.side, side);

       match side {
           BondSide::Left => bond.label_to_atom.0.clone(),
           BondSide::Right => bond.label_to_atom.1.clone(),
       }
   }

   // -------------------------
   // Public API: insert between
   // -------------------------

   pub fn possible_insertions_between<'a>(
       &self,
       mon1: MonomerIdx,
       mon2: MonomerIdx,
       monomers_db: &'a MonomersDB,
   ) -> Vec<&'a MonomersDB_Entry> {
       let bond = &self.monomer_bonds
           .iter()
           .find(|b| {
	       (b.monomers.0 == mon1 && b.monomers.1 == mon2)
		   || (b.monomers.0 == mon2 && b.monomers.1 == mon1)
	   })
           .unwrap_or_else(|| panic!("bond between monomers {:?} not found",
				     (mon1, mon2)));

       if &bond.bond_templ != &*AMINO_BOND_TEMPLATE {
           return Vec::new();
       }

       get_entries_for_insertion_between(monomers_db)
   }

   pub fn insert_between(
       &mut self,
       mon1: MonomerIdx,
       mon2: MonomerIdx,
       mon_db_entry: &MonomersDB_Entry,
   ) {
       // Find and clone the old bond (we'll remove it from `self.monomer_bonds`).
       let old_bond_idx = self
           ._find_bond_idx_between(mon1, mon2)
           .expect("insert_between: no bond exists between the given monomers");
       let old_bond = self.monomer_bonds[old_bond_idx].clone();

       if !Self::_is_supported_insertion_bond_templ(&old_bond.bond_templ) {
           panic!("insert_between: only supported for C-N amino-like single bonds (arity=1, BondType::Single)");
       }

       // Normalize orientation: old bond is (left_mon) -- (right_mon)
       let left_mon = old_bond.monomers.0;
       let right_mon = old_bond.monomers.1;

       // Make sure the caller’s (mon1, mon2) refers to this bond (either order).
       if !((left_mon == mon1 && right_mon == mon2) || (left_mon == mon2 && right_mon == mon1)) {
           panic!("insert_between: internal error: found bond does not match requested monomer pair");
       }

       // Create the new monomer index.
       let inserted_idx = self._next_monomer_idx();

       // Shift atom ids in the inserted entry to avoid collisions with existing atoms.
       let shift = self._max_atom_id() + 1;
       let inserted_entry = mon_db_entry.shift_atom_ids(shift);

       // (Caller promised mon_db_entry came from `possible_insertions_between`.)
       let required_profile = Self::_required_deg2_insertion_profile(&old_bond.bond_templ);
       debug_assert_eq!(inserted_entry.bonds_by_bs.get_profile(), required_profile);

       // Pull out the label->atom maps for the inserted monomer on both sides.
       // In the new left bond, inserted monomer is on the Right.
       // In the new right bond, inserted monomer is on the Left.
       let inserted_left_map =
           Self::_entry_label_map_for_side(&inserted_entry, &old_bond.bond_templ, BondSide::Left);
       let inserted_right_map =
           Self::_entry_label_map_for_side(&inserted_entry, &old_bond.bond_templ, BondSide::Right);

       // Old bond’s neighbor-side maps (these must be preserved).
       let old_left_map = old_bond.label_to_atom.0.clone();
       let old_right_map = old_bond.label_to_atom.1.clone();

       let new_left_bond = Bond {
           bond_templ: old_bond.bond_templ.clone(),
           monomers: (left_mon, inserted_idx),
           label_to_atom: (old_left_map, inserted_right_map),
       };

       let new_right_bond = Bond {
           bond_templ: old_bond.bond_templ.clone(),
           monomers: (inserted_idx, right_mon),
           label_to_atom: (inserted_left_map, old_right_map),
       };

       // Remove the old bond (unique between a pair, per your assumption) and insert the two new ones.
       self.monomer_bonds.retain(|b| {
           !((b.monomers.0 == left_mon && b.monomers.1 == right_mon)
               || (b.monomers.0 == right_mon && b.monomers.1 == left_mon))
       });
       self.monomer_bonds.push(new_left_bond);
       self.monomer_bonds.push(new_right_bond);

       // Insert the new monomer itself.
       self.monomers
           .insert(inserted_idx, inserted_entry.monomer.clone());
   }

}
