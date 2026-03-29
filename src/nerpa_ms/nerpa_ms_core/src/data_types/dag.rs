use anyhow::{bail, Result};
use serde::Deserialize;

use crate::data_types::common_types::{MonomerCode, MonomerIdx};

pub type VertexId = usize;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Edge {
    pub to: VertexId,
    pub weight: u8, // 0 or 1
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VertexLabel {
    pub mon_code: MonomerCode,
    pub mon_idx: MonomerIdx,
}

#[derive(Debug, Clone)]
pub struct DAG {
    pub nrp_variant_id: String,
    pub labels: Vec<Option<VertexLabel>>, // start/finish must be None
    pub out_edges: Vec<Vec<Edge>>,        // out_edges[v] = list of edges from vertex v
    pub start: VertexId,                  // 0
    pub finish: VertexId,                 // labels.len() - 1
}

#[derive(Debug, Clone, Deserialize)]
struct DAG_JSON {
    pub nrp_variant_id: String,
    pub labels: Vec<Option<[usize; 2]>>, // start/finish must be None
    /// Each edge is [u, v, w] where w ∈ {0,1}
    pub edges: Vec<[usize; 3]>,
}

impl DAG {
    pub fn vertex_count(&self) -> usize {
        self.labels.len()
    }

    /// Build a `Dag` from a serde_json::Value that matches the prototype schema.
    pub fn from_json_value(val: &serde_json::Value) -> Result<DAG> {
        let raw: DAG_JSON = serde_json::from_value(val.clone())?;

        if raw.labels.len() < 2 {
            bail!("DAG: labels must have length >= 2");
        }

        let start: VertexId = 0;
        let finish: VertexId = raw.labels.len() - 1;

        let labels = raw
            .labels
            .iter()
            .map(|opt_arr| {
                opt_arr.map(|[mon_code, mon_idx]| VertexLabel {
                    mon_code: MonomerCode(mon_code as u32),
                    mon_idx: MonomerIdx(mon_idx as u32),
                })
            })
            .collect::<Vec<Option<VertexLabel>>>();

        if raw.labels[start].is_some() {
            bail!("DAG: labels[0] (start) must be null");
        }
        if raw.labels[finish].is_some() {
            bail!("DAG: labels[last] (finish) must be null");
        }

        let n = raw.labels.len();
        let mut out_edges: Vec<Vec<Edge>> = vec![Vec::new(); n];

        for (i, [u, v, w]) in raw.edges.iter().copied().enumerate() {
            if u >= n || v >= n {
                bail!("DAG: edge #{i} has out-of-bounds vertex (u={u}, v={v}, labels.len()={n})");
            }
            if w != 0 && w != 1 {
                bail!("DAG: edge #{i} has invalid weight {w} (must be 0 or 1)");
            }

            out_edges[u].push(Edge {
                to: v,
                weight: w as u8,
            });
        }

        Ok(DAG {
            nrp_variant_id: raw.nrp_variant_id,
            labels,
            out_edges,
            start,
            finish,
        })
    }
}
