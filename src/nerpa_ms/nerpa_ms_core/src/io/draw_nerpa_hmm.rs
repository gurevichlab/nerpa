use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;
use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};

use anyhow::{anyhow, bail, ensure, Context, Result};
use minilp::{ComparisonOp, OptimizationDirection, Problem};

use crate::data_types::hmm::HMM;
use crate::data_types::common_types::LogProb;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StateType {
    SkippingModulesAtEnd,
    InsertAtEnd,
    MatchingFinished,
    Insert,
    ChooseIfFinishMatching,
    Match,
    MatchPossibleAssemblyFinish,
    MatchLastModule,
    MatchingContinued,
    Initial,
    ModuleSubgraphRoot,
    Final,
    ChooseIfIterateModule,
    ChooseIfIterateGene,
    InsertAtStart,
    ChooseIfStartMatching,
    SkippingModulesAtStartFinished,
}

impl StateType {
    /// Parses Rust `hmm.state_labels[idx]`, which matches Python `state_idx_to_label(idx)`:
    /// - MODULE_SUBGRAPH_ROOT: "{idx}:F{fragment}:{gene}:{a_domain}"
    /// - otherwise: "{idx}:{STATE_TYPE_NAME}"
    pub fn from_label(label: &str) -> Result<Self> {
        let (_idx, rest) = label
            .split_once(':')
            .ok_or_else(|| anyhow!("Invalid state label (no ':'): {label:?}"))?;

        if rest.starts_with('F') {
            return Ok(StateType::ModuleSubgraphRoot);
        }

        let st = match rest {
            "SKIPPING_MODULES_AT_END" => StateType::SkippingModulesAtEnd,
            "INSERT_AT_END" => StateType::InsertAtEnd,
            "MATCHING_FINISHED" => StateType::MatchingFinished,
            "INSERT" => StateType::Insert,
            "CHOOSE_IF_FINISH_MATCHING" => StateType::ChooseIfFinishMatching,
            "MATCH" => StateType::Match,
            "MATCH_POSSIBLE_ASSEMBLY_FINISH" => StateType::MatchPossibleAssemblyFinish,
            "MATCH_LAST_MODULE" => StateType::MatchLastModule,
            "MATCHING_CONTINUED" => StateType::MatchingContinued,
            "INITIAL" => StateType::Initial,
            "FINAL" => StateType::Final,
            "CHOOSE_IF_ITERATE_MODULE" => StateType::ChooseIfIterateModule,
            "CHOOSE_IF_ITERATE_GENE" => StateType::ChooseIfIterateGene,
            "INSERT_AT_START" => StateType::InsertAtStart,
            "CHOOSE_IF_START_MATCHING" => StateType::ChooseIfStartMatching,
            "SKIPPING_MODULES_AT_START_FINISHED" => StateType::SkippingModulesAtStartFinished,
            other => bail!("Unknown state type label suffix: {other:?} (full label: {label:?})"),
        };
        Ok(st)
    }
}

fn dot_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 8);
    for ch in s.chars() {
        match ch {
            '\\' => out.push_str("\\\\"),
            '"' => out.push_str("\\\""),
            '\n' => out.push_str("\\n"),
            '\r' => {}
            _ => out.push(ch),
        }
    }
    out
}

fn state_color(st: StateType) -> &'static str {
    match st {
        StateType::ModuleSubgraphRoot => "lightblue",
        StateType::Initial | StateType::Final => "grey",
        _ => "white",
    }
}

/// Replicates the Python LP layout (same constraints/objective, but solved via `minilp`).
fn create_nodes_layout_lp(
    layers: &[Vec<usize>],
    forward_edges: &[(usize, usize)],
    labels: &[String],
    char_width: f64,
    label_padding: f64,
    edge_min_dx: f64,
    layer_gap: f64,
    layer_y_gap: f64,
) -> Result<HashMap<usize, (f64, f64)>> {
    let n_states = labels.len();
    ensure!(n_states > 0, "HMM has no states");

    // Ensure all states are assigned to a layer.
    let mut layer_of: HashMap<usize, usize> = HashMap::new();
    for (li, layer) in layers.iter().enumerate() {
        for &n in layer {
            ensure!(n < n_states, "Layer contains out-of-range node idx {n}");
            ensure!(
                layer_of.insert(n, li).is_none(),
                "State {n} assigned to multiple layers"
            );
        }
    }
    ensure!(
        layer_of.len() == n_states,
        "Not all states were assigned to layers: assigned={}, total={}",
        layer_of.len(),
        n_states
    );

    let node_size = |i: usize| -> f64 { (labels[i].len() as f64) * char_width + label_padding };

    // Objective in Python: minimize x_last - x_first, with x_first fixed to 0.
    // Therefore minimize x_last.
    let mut pb = Problem::new(OptimizationDirection::Minimize);

    let first = 0usize;
    let last = n_states - 1;

    let mut x_vars = Vec::with_capacity(n_states);
    for i in 0..n_states {
        let cost = if i == last { 1.0 } else { 0.0 };
        x_vars.push(pb.add_var(cost, (f64::NEG_INFINITY, f64::INFINITY)));
    }

    // Fix translation: x_first == 0
    pb.add_constraint(vec![(x_vars[first], 1.0)], ComparisonOp::Eq, 0.0);

    // Edge order constraints: x_v - x_u >= edge_min_dx
    for &(u, v) in forward_edges {
        ensure!(u < n_states && v < n_states, "Edge has out-of-range node");
        pb.add_constraint(
            vec![(x_vars[v], 1.0), (x_vars[u], -1.0)],
            ComparisonOp::Ge,
            edge_min_dx,
        );
    }

    // Same-layer non-overlap constraints for adjacent-by-index nodes in that layer
    for layer in layers {
        for w in layer.windows(2) {
            let u = w[0];
            let v = w[1];
            let sep = (node_size(u) / 2.0) + layer_gap + (node_size(v) / 2.0);
            pb.add_constraint(
                vec![(x_vars[v], 1.0), (x_vars[u], -1.0)],
                ComparisonOp::Ge,
                sep,
            );
        }
    }

    let solution = pb.solve().map_err(|e| anyhow!("LP solve failed: {e:?}"))?;

    let mut pos = HashMap::with_capacity(n_states);
    for i in 0..n_states {
        let xi = solution[x_vars[i]];
        let li = *layer_of.get(&i).expect("layer_of missing state");
        let yi = ((layers.len() - li) as f64) * layer_y_gap;
        pos.insert(i, (xi, yi));
    }
    Ok(pos)
}

/// Draws an HMM as DOT and optionally renders to SVG via graphviz.
///
/// Rendering command (assumes graphviz installed):
///   dot -Kneato -n2 -Tsvg -o <output_path>
pub fn draw_nerpa_hmm(
    hmm: &HMM,
    output_path: Option<&Path>,
    highlight_path: Option<&[usize]>,
    edge_weights: bool,
) -> Result<String> {
    let n_states = hmm.state_labels.len();
    ensure!(n_states > 0, "HMM has no states");
    ensure!(
        hmm.transitions.len() == n_states,
        "Expected transitions.len() == state_labels.len()"
    );
    ensure!(
        hmm.emissions.len() == n_states,
        "Expected emissions.len() == state_labels.len()"
    );

    // Parse state types from labels
    let mut state_types = Vec::with_capacity(n_states);
    for s in &hmm.state_labels {
        state_types.push(StateType::from_label(s)?);
    }

    // Layer definitions (same as Python)
    let layer_types: Vec<Vec<StateType>> = vec![
        vec![StateType::SkippingModulesAtEnd],
        vec![StateType::InsertAtEnd],
        vec![StateType::MatchingFinished],
        vec![StateType::Insert, StateType::ChooseIfFinishMatching],
        vec![
            StateType::Match,
            StateType::MatchPossibleAssemblyFinish,
            StateType::MatchLastModule,
            StateType::MatchingContinued,
        ],
        vec![
            StateType::Initial,
            StateType::ModuleSubgraphRoot,
            StateType::Final,
            StateType::ChooseIfIterateModule,
            StateType::ChooseIfIterateGene,
        ],
        vec![StateType::InsertAtStart],
        vec![
            StateType::ChooseIfStartMatching,
            StateType::SkippingModulesAtStartFinished,
        ],
    ];

    let mut layers: Vec<Vec<usize>> = vec![Vec::new(); layer_types.len()];
    for idx in 0..n_states {
        let st = state_types[idx];
        let mut placed = false;
        for (layer_idx, lts) in layer_types.iter().enumerate() {
            if lts.contains(&st) {
                layers[layer_idx].push(idx);
                placed = true;
                break;
            }
        }
        ensure!(
            placed,
            "State {idx} with type {st:?} not placed into any layer"
        );
    }

    // Forward edges for LP constraints + all edges for drawing
    let mut forward_edges: Vec<(usize, usize)> = Vec::new();
    let mut all_edges: Vec<(usize, usize, LogProb)> = Vec::new();
    let mut edge_set: HashSet<(usize, usize)> = HashSet::new();

    for (from_idx, outs) in hmm.transitions.iter().enumerate() {
        for &(to_idx_u32, lp) in outs {
            let to_idx = to_idx_u32 as usize;
            ensure!(
                to_idx < n_states,
                "Transition to_idx out of range: {to_idx}"
            );
            all_edges.push((from_idx, to_idx, lp));
            edge_set.insert((from_idx, to_idx));
            if to_idx > from_idx {
                forward_edges.push((from_idx, to_idx));
            }
        }
    }

    // Validate highlight_path + compute highlighted edges
    let mut highlighted: HashSet<(usize, usize)> = HashSet::new();
    if let Some(path) = highlight_path {
        for &v in path {
            ensure!(v < n_states, "highlight_path vertex out of range: {v}");
        }
        for w in path.windows(2) {
            let (u, v) = (w[0], w[1]);
            ensure!(
                edge_set.contains(&(u, v)),
                "highlight_path edge does not exist: {u} -> {v}"
            );
            highlighted.insert((u, v));
        }
    }

    // LP layout
    let pos = create_nodes_layout_lp(
        &layers,
        &forward_edges,
        &hmm.state_labels,
        0.6, // char_width
        0.8, // label_padding
        0.5, // edge_min_dx
        0.6, // layer_gap
        5.0, // layer_y_gap
    )?;

    // Pivot layer: the layer containing INITIAL/FINAL (for port direction rule)
    let pivot_layer_idx = layer_types
        .iter()
        .position(|lts| lts.contains(&StateType::Initial) || lts.contains(&StateType::Final))
        .ok_or_else(|| anyhow!("No layer contains INITIAL/FINAL"))?;

    let mut layer_of = vec![usize::MAX; n_states];
    for (li, layer) in layers.iter().enumerate() {
        for &v in layer {
            layer_of[v] = li;
        }
    }

    // Build DOT
    let mut dot = String::new();
    writeln!(&mut dot, "digraph HMM {{")?;
    writeln!(&mut dot, "  graph [")?;
    writeln!(&mut dot, r#"    overlap="false","#)?;
    writeln!(&mut dot, r#"    splines="polyline","#)?;
    writeln!(&mut dot, r#"    esep="+1","#)?;
    writeln!(&mut dot, r#"    sep="+10""#)?;
    writeln!(&mut dot, "  ];")?;
    writeln!(&mut dot, r#"  node [shape="ellipse", style="filled"];"#)?;
    writeln!(&mut dot, r#"  edge [arrowhead="vee"];"#)?;

    for idx in 0..n_states {
        let (x, y) = *pos.get(&idx).ok_or_else(|| anyhow!("Missing pos for node {idx}"))?;
        let label = dot_escape(&hmm.state_labels[idx]);
        let fill = state_color(state_types[idx]);
        writeln!(
            &mut dot,
            r#"  "{idx}" [label="{label}", fillcolor="{fill}", pos="{x},{y}!"];"#
        )?;
    }

    for (from, to, lp) in all_edges {
        let is_hl = highlighted.contains(&(from, to));
        let color = if is_hl { "red" } else { "black" };
        let penwidth = if is_hl { "2" } else { "1" };
        let arrowsize = if is_hl { "1.5" } else { "1" };

        let mut attrs = Vec::<String>::new();
        attrs.push(format!(r#"color="{color}""#));
        attrs.push(format!(r#"penwidth="{penwidth}""#));
        attrs.push(format!(r#"arrowsize="{arrowsize}""#));

        if edge_weights {
            attrs.push(format!(r#"label="{:.2}""#, lp)); // raw logprob as requested
        }

        // Simpler port rule (applied to self-loops only):
        // below pivot layer => down ("s"), otherwise up ("n")
        if from == to {
            let from_layer = layer_of[from];
            ensure!(from_layer != usize::MAX, "Missing layer for node {from}");
            let port = if from_layer > pivot_layer_idx { "s" } else { "n" };
            attrs.push(format!(r#"headport="{port}""#));
            attrs.push(format!(r#"tailport="{port}""#));
        }

        writeln!(
            &mut dot,
            r#"  "{from}" -> "{to}" [{}];"#,
            attrs.join(", ")
        )?;
    }

    writeln!(&mut dot, "}}")?;

    // Optionally render to SVG using graphviz
    if let Some(out) = output_path {
        let mut child = Command::new("dot")
            .arg("-Kneato")
            .arg("-n2")
            .arg("-Tsvg")
            .arg("-o")
            .arg(out)
            .stdin(Stdio::piped())
            .stdout(Stdio::null())
            .stderr(Stdio::inherit())
            .spawn()
            .context("Failed to spawn graphviz 'dot' (is graphviz installed and on PATH?)")?;

        {
            let mut stdin = child.stdin.take().context("Failed to open dot stdin")?;
            stdin
                .write_all(dot.as_bytes())
                .context("Failed to write DOT to graphviz stdin")?;
        }

        let status = child.wait().context("Failed to wait on graphviz process")?;
        ensure!(status.success(), "graphviz 'dot' failed with status {status}");
    }

    Ok(dot)
}
