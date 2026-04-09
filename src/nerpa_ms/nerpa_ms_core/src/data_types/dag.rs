use crate::data_types::graph_modifications::GraphModification;

use crate::data_types::common_types::MonomerCode;
use serde::de::{Error as DeError, IgnoredAny};
use serde::{Deserialize, Deserializer};

pub type VertexId = usize;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Edge<'a> {
    pub to: VertexId,
    pub weight: u8, // 0 or 1
    pub modification: Option<GraphModification<'a>>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct VertexLabel {
    pub monomer_code: Option<MonomerCode>,
    pub name: String, // used only for graph drawing
}

#[derive(Debug, Clone)]
pub struct DAG<'a> {
    pub nrp_variant_id: String,
    pub labels: Vec<VertexLabel>,      // labels[v] = label of vertex v
    pub out_edges: Vec<Vec<Edge<'a>>>, // out_edges[v] = list of edges from vertex v
    pub start: VertexId,               // 0
    pub finish: VertexId,              // labels.len() - 1
}

impl<'a> DAG<'a> {
    pub fn num_nodes(&self) -> usize {
        self.labels.len()
    }
}

// ---- Custom Deserialize impls ----
// Supports only deserializing with all modifications set to None, which is used for testing and drawing

#[derive(Deserialize)]
struct EdgeDe {
    to: VertexId,
    weight: u8,
}

impl EdgeDe {
    fn into_edge<'a>(self) -> Edge<'a> {
        Edge {
            to: self.to,
            weight: self.weight,
            modification: None,
        }
    }
}

#[derive(Deserialize)]
struct DAGDe {
    nrp_variant_id: String,
    labels: Vec<VertexLabel>,
    out_edges: Vec<Vec<EdgeDe>>,
    start: VertexId,
    finish: VertexId,
}

impl<'de> Deserialize<'de> for DAG<'static> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let helper = DAGDe::deserialize(deserializer)?;

        if helper.out_edges.len() != helper.labels.len() {
            return Err(D::Error::custom("out_edges.len() must equal labels.len()"));
        }

	let out_edges = helper
	    .out_edges
	    .into_iter()
	    .map(|edges_de| {
		edges_de
		    .into_iter()
		    .map(|e_de| e_de.into_edge())
		    .collect::<Vec<Edge>>()
	    })
	    .collect::<Vec<Vec<Edge>>>();

        Ok(DAG {
            nrp_variant_id: helper.nrp_variant_id,
            labels: helper.labels,
            out_edges: out_edges,
            start: helper.start,
            finish: helper.finish,
        })
    }
}
