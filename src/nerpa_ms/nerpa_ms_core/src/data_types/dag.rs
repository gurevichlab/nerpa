use anyhow::{bail, Result};
use serde::Deserialize;

use crate::data_types::common_types::{MonomerCode, MonomerIdx};
use crate::data_types::graph_modifications::{GraphModification};

pub type VertexId = usize;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Edge<'a> {
    pub to: VertexId,
    pub weight: u8, // 0 or 1
    pub modification: Option<GraphModification<'a>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VertexLabel {
    pub mon_code: MonomerCode,
    pub mon_idx: MonomerIdx,
}

#[derive(Debug, Clone)]
pub struct DAG<'a> {
    pub nrp_variant_id: String,
    pub labels: Vec<Option<String>>, // used only for graph drawing
    pub out_edges: Vec<Vec<Edge<'a>>>,        // out_edges[v] = list of edges from vertex v
    pub start: VertexId,                  // 0
    pub finish: VertexId,                 // labels.len() - 1
}
