use anyhow::Result;
use std::env;
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use crate::data_types::dag::{VertexId, DAG};

pub struct Draw_DAG_Config {
	pub node_indexes: bool,
}

impl DAG<'_> {
    pub fn to_dot(&self, cfg: &Draw_DAG_Config) -> String {
        // Prefer labels.len(), but fall back safely if you ever forget to fill it.
        let mut out = String::new();
        out.push_str("digraph G {\n");
        out.push_str("  rankdir=LR;\n");
        out.push_str("  node [shape=ellipse, fontname=\"monospace\"];\n");
        out.push_str("  edge [fontname=\"monospace\"];\n\n");

        // Nodes
        for v in 0..self.num_nodes() {
            let label = {
                if let Some(label) = &self.labels[v] {
		    if cfg.node_indexes {v.to_string() + "." + label}
		    else { label.clone() }
                } else {
		    if cfg.node_indexes {v.to_string() }
		    else { String::new() }
                }
            };

            let _ = writeln!(&mut out, "  v{v} [label=\"{label}\"];");
        }

        out.push('\n');

        // Edges
        for (from, edges) in self.out_edges.iter().enumerate() {
            for e in edges {
                let to: VertexId = e.to;
                // Ignore e.modification for drawing, as requested.
                let _ = writeln!(&mut out, "  v{from} -> v{to} [label=\"{}\"];\n", e.weight);
            }
        }

        out.push_str("}\n");
        out
    }

    pub fn draw_svg(&self, out: &Path, cfg: &Draw_DAG_Config) -> Result<()> {
        // Ensure output directory exists
        if let Some(parent) = out.parent() {
            fs::create_dir_all(parent).expect("create output directory");
        }

        // Write DOT to a temp file
        let dot = self.to_dot(cfg);
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time went backwards")
            .as_nanos();

	let dot_path = out.parent().unwrap().join(self.nrp_variant_id.clone() + ".dot");

        fs::write(&dot_path, dot).expect("write .dot file");

        // Render SVG using Graphviz
        let status = Command::new("dot")
            .args([
                "-Tsvg",
                dot_path.to_str().expect("non-utf8 temp path"),
                "-o",
                out.to_str().expect("non-utf8 output path"),
            ])
            .status();

        // Best-effort cleanup of temp file
        let _ = fs::remove_file(&dot_path);
        match status {
            Ok(s) if s.success() => Ok(()),
            Ok(s) => Err(anyhow::anyhow!("dot exited with status {s}")),
            Err(e) => Err(anyhow::anyhow!("Failed to execute dot: {e}")),
        }
    }
}
