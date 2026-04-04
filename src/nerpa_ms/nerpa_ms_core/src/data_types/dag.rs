use crate::data_types::graph_modifications::GraphModification;

use crate::data_types::common_types::MonomerCode;

pub type VertexId = usize;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Edge<'a> {
    pub to: VertexId,
    pub weight: u8, // 0 or 1
    pub modification: Option<GraphModification<'a>>,
}

#[derive(Debug, Clone)]
pub struct VertexLabel {
	pub monomer_code: Option<MonomerCode>,
	pub name: String,  // used only for graph drawing
}

#[derive(Debug, Clone)]
pub struct DAG<'a> {
    pub nrp_variant_id: String,
    pub labels: Vec<VertexLabel>, // labels[v] = label of vertex v
    pub out_edges: Vec<Vec<Edge<'a>>>, // out_edges[v] = list of edges from vertex v
    pub start: VertexId,               // 0
    pub finish: VertexId,              // labels.len() - 1
}

impl<'a> DAG<'a> {
    pub fn num_nodes(&self) -> usize {
        self.labels.len()
    }
}
