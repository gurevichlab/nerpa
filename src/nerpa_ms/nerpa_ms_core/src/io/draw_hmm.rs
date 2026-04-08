use anyhow::Result;
use std::collections::HashSet;
use std::fmt::Write as _;
use std::fs;
use std::path::Path;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

// Adjust this import path to wherever your HMM types live.
use crate::data_types::hmm::{HMM, StateIdx};

pub struct Draw_HMM_Config {
    pub state_indexes: bool,
}

impl HMM {
    pub fn to_dot(
        &self,
        cfg: &Draw_HMM_Config,
        highlight_path: Option<Vec<StateIdx>>,
    ) -> String {
        let highlight_edges: HashSet<(StateIdx, StateIdx)> = match highlight_path {
            Some(path) => path.windows(2).map(|w| (w[0], w[1])).collect(),
            None => HashSet::new(),
        };

        let n_states = self.transitions.len();

        let mut out = String::new();
        out.push_str("digraph HMM {\n");
        out.push_str("  rankdir=LR;\n");
        out.push_str("  node [shape=ellipse, fontname=\"monospace\"];\n");
        out.push_str("  edge [fontname=\"monospace\"];\n\n");

        // Nodes
        for s in 0..n_states {
            let name = state_name(s, n_states);

            let label = if cfg.state_indexes {
                format!("{s}.{name}")
            } else {
                name.to_string()
            };

            let is_emitting = !self.emissions[s].is_empty();

            if is_emitting {
                let _ = writeln!(
                    &mut out,
                    "  s{s} [label=\"{label}\", style=filled, fillcolor=\"lightblue\"];"
                );
            } else {
                let _ = writeln!(&mut out, "  s{s} [label=\"{label}\"];");
            }
        }

        out.push('\n');

        // Edges
        for (from, edges) in self.transitions.iter().enumerate() {
            for &(to, lp) in edges {
                let attrs = if highlight_edges.contains(&(from, to)) {
                    format!(
                        "label=\"{:.3}\", color=\"red\", fontcolor=\"red\", penwidth=2",
                        lp
                    )
                } else {
                    format!("label=\"{:.3}\"", lp)
                };

                let _ = writeln!(&mut out, "  s{from} -> s{to} [{attrs}];");
            }
        }

        out.push_str("}\n");
        out
    }

    pub fn draw_svg(
        &self,
        out: &Path,
        cfg: &Draw_HMM_Config,
        highlight_path: Option<Vec<StateIdx>>,
    ) -> Result<()> {
        // Ensure output directory exists
        if let Some(parent) = out.parent() {
            fs::create_dir_all(parent).expect("create output directory");
        }

        // Write DOT to a temp file
        let dot = self.to_dot(cfg, highlight_path);

        // Give the temp file a stable-ish name; also avoid weird chars.
        let dot_stem = sanitize_for_filename(&format!("{:?}", self.bgc_variant_id));

        // If bgc_variant_id is super long / empty / whatever, still produce *some* name.
        let dot_stem = if dot_stem.is_empty() { "hmm".to_string() } else { dot_stem };

        let dot_path = out
            .parent()
            .unwrap()
            .join(format!("{dot_stem}.dot"));

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

fn state_name(state: usize, n_states: usize) -> &'static str {
    if state == 0 {
        "START"
    } else if state + 1 == n_states {
        "FINISH"
    } else {
        ""
    }
}

fn sanitize_for_filename(s: &str) -> String {
    // Cheap and cheerful.
    s.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '-' || c == '_' || c == '.' {
                c
            } else {
                '_'
            }
        })
        .collect()
}
