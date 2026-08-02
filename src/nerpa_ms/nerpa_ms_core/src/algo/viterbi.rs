use std::collections::{HashMap, HashSet};

use crate::data_types::bond_consts::{AMINO_BINDING_SITE_C, AMINO_BINDING_SITE_N};
use crate::data_types::bonds::{BindingSiteType, Bond};
use crate::data_types::common_types::{LogOdds, MonomerIdx};
use crate::data_types::mod_graph::{ModGraph, Edge, VertexLabel, VertexId};
use crate::data_types::graph_modifications::{GraphModification, InsertionSite};
use crate::data_types::monomer_graph::{MonomerFeatures, MonomerGraph, Monomer};
use crate::data_types::parsed_rban_record::{MonomerInfo, NorineMonomerName};
use crate::data_types::{monomers_db::MonomersDB, parsed_rban_record::Parsed_rBAN_Record};
use crate::data_types::hmm::HMM;

use super::gen_new_variants::{NewVariantWithOptPaths, generate_new_variants_with_opt_paths};

impl ModGraph<'_> {
    pub fn trivial(nrp_id: String, monomers: &[&Monomer]) -> Self {
	let labels: Vec<VertexLabel> = {
	    (0..monomers.len())
		.zip(monomers.iter())
		.map(|(i, mon)| VertexLabel {
		    monomer_code: Some(mon.features.mon_code.clone()),
		    name: mon.features.name.0.clone(),
		})
		.collect()
	};

	let mut out_edges: Vec<Vec<Edge<'_>>> = vec![vec![]; monomers.len()];
	for i in 0..monomers.len() - 1 {
	    out_edges[i].push(Edge {
		to: i + 1,
		weight: 0,
		modification: None,
	    });
	}  

	ModGraph {
	    nrp_variant_id: nrp_id,
	    labels,
	    out_edges,
	    start: 0,
	    finish: monomers.len() - 1,
	}
    }
}

impl HMM {
    pub fn score(
	&self,
	monomer_graph: &MonomerGraph,
	linearization: &[MonomerIdx],
    ) -> LogOdds {
    let monomers: Vec<&Monomer> = {
	linearization
	    .iter()
	    .map(|mon_idx| monomer_graph.monomers
		 .get(mon_idx)
		 .expect("MonomerIdx in linearization not found in monomer graph")
	    )
	    .collect()
    };
    let mod_graph = ModGraph::trivial(
	monomer_graph.compound_id.clone(),
	&monomers
    );

    let new_variants_with_opt_paths: Vec<NewVariantWithOptPaths> = {
	generate_new_variants_with_opt_paths(
	    &self,
	    monomer_graph,
	    &mod_graph,
	    0, // max_edits
	    1, // max_variants
	    &HashMap::new(), // monomers_db (not needed since we are not applying any modifications)
	    false, // debug_stdout
	)
    };

    new_variants_with_opt_paths
	.into_iter()
	.next()
	.expect("Variant generation failed to produce any variants")
	.new_variant
	.score
    }
}
