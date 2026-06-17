use super::find_target::TargetSearchResults;
use anyhow::Result;

/// A single distance bucket like "0", "1", "2", "3", "4+",
/// along with all results that fell into that bucket.
pub struct DistanceBucket<'a> {
    pub label: String,
    pub items: Vec<&'a TargetSearchResults>,
}

/// Split results into distance buckets:
/// - one bucket per exact distance in 0..=max_explicit_distance
/// - plus one overflow bucket for distances >= max_explicit_distance + 1 (labeled "{start}+")
///
/// Returns buckets in a stable order:
/// 0, 1, 2, ..., max_explicit_distance, "{max_explicit_distance+1}+"
pub fn split_results_by_distance<'a>(
    results: &'a [TargetSearchResults],
    max_explicit_distance: u32,
) -> Vec<DistanceBucket<'a>> {
    // Pre-create all buckets so empty ones still appear (useful for stable plots).
    let mut buckets: Vec<DistanceBucket<'a>> = {
        let mut out = Vec::new();

        for d in 0..=max_explicit_distance {
            out.push(DistanceBucket {
                label: d.to_string(),
                items: Vec::new(),
            });
        }

        let overflow_start = max_explicit_distance + 1;
        out.push(DistanceBucket {
            label: format!("{overflow_start}+"),
            items: Vec::new(),
        });

        out
    };

    for result in results {
        let distance = result.features.linearization_distance;

        let bucket_index: usize = if distance <= max_explicit_distance {
            distance as usize
        } else {
            // overflow bucket
            (max_explicit_distance + 1) as usize
        };

        buckets[bucket_index].items.push(result);
    }

    buckets
}

/// One bar inside a cluster (within a single bucket).
///
/// `predicate == true` means "this result should be counted in this bar".
pub struct ClusterBar {
    pub label: String,
    pub predicate: Box<dyn Fn(&TargetSearchResults) -> bool + Send + Sync>,
}

/// Defines how to compute bar heights for one bucket.
///
/// If `cumulative == true`, then a result contributes to bar N only if it also
/// contributed to all bars 0..N-1.
///
/// This matches your “each set is a superset of the previous one” requirement
/// and guarantees non-increasing bar heights.
pub struct ClusterBarScheme {
    pub bars: Vec<ClusterBar>,
    pub cumulative: bool,
}

impl ClusterBarScheme {
    /// Labels for legends etc.
    pub fn bar_labels(&self) -> Vec<String> {
        self.bars.iter().map(|b| b.label.clone()).collect()
    }

    /// Compute bar heights for a single bucket.
    pub fn heights_for_bucket(&self, bucket_items: &[&TargetSearchResults]) -> Vec<u32> {
        let mut heights = vec![0u32; self.bars.len()];

        if self.cumulative {
            for item in bucket_items {
                self.add_item_cumulative(&mut heights, item);
            }
	} else {
	    for item in bucket_items {
                self.add_item_independent(&mut heights, item);
            }
        }

        heights
    }

    fn add_item_independent(&self, heights: &mut [u32], item: &TargetSearchResults) {
        for (i, bar) in self.bars.iter().enumerate() {
            if (bar.predicate)(item) {
                heights[i] += 1;
            }
        }
    }

    fn add_item_cumulative(&self, heights: &mut [u32], item: &TargetSearchResults) {
        for (i, bar) in self.bars.iter().enumerate() {
            if (bar.predicate)(item) {
                heights[i] += 1;
            } else {
                // Once a bar fails, all later (stricter) bars fail too.
                break;
            }
        }
    }
}

/// Plot-ready data for one bucket:
/// the plotting code should only care about:
/// - `bucket_label` (x-axis category)
/// - `bar_heights` (cluster bars left-to-right)
pub struct BucketClusterData {
    pub bucket_label: String,
    pub bar_heights: Vec<u32>,
}

/// Convert distance buckets into plot-ready cluster data.
pub fn compute_bucket_cluster_data(
    buckets: &[DistanceBucket<'_>],
    scheme: &ClusterBarScheme,
) -> Vec<BucketClusterData> {
    let mut out = Vec::new();

    for bucket in buckets {
        let heights = scheme.heights_for_bucket(&bucket.items);

        out.push(BucketClusterData {
            bucket_label: bucket.label.clone(),
            bar_heights: heights,
        });
    }

    out
}


use std::{path::Path, sync::LazyLock};

pub static TARGET_SEARCH_BAR_SCHEME: LazyLock<ClusterBarScheme> = LazyLock::new(|| {
    ClusterBarScheme {
        cumulative: true,
        bars: vec![
            ClusterBar {
                label: "total".to_string(),
                predicate: Box::new(|_r: &TargetSearchResults| true),
            },
            ClusterBar {
                label: "new_monomers_supported".to_string(),
                predicate: Box::new(|r: &TargetSearchResults| r.features.new_monomers_supported),
            },
            ClusterBar {
                label: "edits_applicable".to_string(),
                predicate: Box::new(|r: &TargetSearchResults| r.features.edits_applicable),
            },
            ClusterBar {
                label: "found_match".to_string(),
                predicate: Box::new(|r: &TargetSearchResults| r.found_match.is_some()),
            },
            ClusterBar {
                label: "mon_graph_reconstructable".to_string(),
                predicate: Box::new(|r: &TargetSearchResults| {
                    r.features.target_mon_graph_reconstructable
                }),
            },
            ClusterBar {
                label: "mon_masses_reconstructable".to_string(),
                predicate: Box::new(|r: &TargetSearchResults| {
                    r.features.target_mon_masses_reconstructable
                }),
            },
            ClusterBar {
                label: "smiles_reconstructable".to_string(),
                predicate: Box::new(|r: &TargetSearchResults| {
                    r.features.target_smiles_reconstructable
                }),
            },
        ],
    }
});


// Payload structure for JSON output to be consumed by the plotting code in Python.
use serde::Serialize;

#[derive(Debug, Serialize)]
pub struct PlotPayload {
    pub title: String,
    pub bucket_labels: Vec<String>,
    pub bar_labels: Vec<String>,
    pub heights: Vec<Vec<u32>>,
}

pub fn build_plot_payload(
    title: String,
    buckets: &[BucketClusterData],
    bar_labels: Vec<String>,
) -> PlotPayload {
    let bucket_labels = buckets.iter().map(|b| b.bucket_label.clone()).collect();

    let heights = buckets
        .iter()
        .map(|b| b.bar_heights.clone())
        .collect();

    PlotPayload {
        title,
        bucket_labels,
        bar_labels,
        heights,
    }
}

pub fn write_plot_payload_json(
    payload: &PlotPayload,
    output_json_path: &Path,
) -> Result<()> {
    let json = serde_json::to_string_pretty(payload)?;
    std::fs::create_dir_all(output_json_path.parent().unwrap())?;
    std::fs::write(output_json_path, json)?;
    Ok(())
}

pub fn plot(
    results: &[TargetSearchResults],
    nerpa_ms_dir: &Path,
    svg_output_path: &Path,
) -> Result<()> {
    let buckets = split_results_by_distance(results, 4);
    let cluster_data = compute_bucket_cluster_data(&buckets, &TARGET_SEARCH_BAR_SCHEME);
    let bar_labels = TARGET_SEARCH_BAR_SCHEME.bar_labels();
    let payload = build_plot_payload("Find Siblings Benchmark".to_string(), &cluster_data, bar_labels);
    let input_for_drawing_path = nerpa_ms_dir
	.join("tmp")
	.join("find_siblings_plot_data.json");
    write_plot_payload_json(&payload, &input_for_drawing_path)?;

    let plot_script_path = nerpa_ms_dir
	.join("scripts")
	.join("plot_find_siblings.py");

    // Call the plotting script (assumes Python environment with necessary libraries is set up).
    std::process::Command::new("python")
	.arg(plot_script_path)
	.arg(input_for_drawing_path)
	.arg(svg_output_path)
	.status()?;
    
    Ok(())
}
