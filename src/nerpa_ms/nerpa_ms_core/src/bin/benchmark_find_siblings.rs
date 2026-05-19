#![allow(non_snake_case, non_camel_case_types, non_upper_case_globals)]

use std::path::PathBuf;

use clap::Parser;

#[derive(Debug, Parser)]
#[command(name = "benchmark_find_sibglings")]
pub struct Cli {
    /// Path to Nerpa results on MIBiG: each BGC against its respective NRPs in its own folder
    #[arg(long)]
    pub nerpa_results: PathBuf,

    /// Path to a txt file with the pairwise distances between compounds in the MIBiG dataset
    #[arg(long)]
    pub compound_distances: PathBuf,


    /// Path to the output tsv with the benchmarking results
    #[arg(long)]
    pub out: PathBuf,
}

use anyhow::Result;

use std::collections::HashMap;

fn get_compound_distances(path: &PathBuf) -> Result<HashMap<(String, String), u32>> {
    let text = std::fs::read_to_string(&cli.compound_distances)?;

    let mut lines = text.lines().map(str::trim).filter(|l| !l.is_empty());

    let header_line = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("Empty compound distances file"))?;
    let header_ids: Vec<String> = header_line
        .split_whitespace()
        .map(|s| s.to_string())
        .collect();

    let mut map: HashMap<(String, String), u32> = HashMap::new();

    for line in lines {
        let mut it = line.split_whitespace();
        let row_id = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("Malformed line in compound distances: '{line}'"))?
            .to_string();

        let values: Vec<&str> = it.collect();
        if values.len() != header_ids.len() {
            return Err(anyhow::anyhow!(
                "Wrong number of distances for '{row_id}': expected {}, got {}",
                header_ids.len(),
                values.len()
            ));
        }

        for (col_id, dist_str) in header_ids.iter().zip(values.into_iter()) {
            let d: i32 = dist_str.parse().map_err(|e| {
                anyhow::anyhow!(
                    "Failed to parse distance '{dist_str}' for ({row_id}, {col_id}): {e}"
                )
            })?;
	    if d >= 0 {  // negative values indicate "infinite" distance
		map.insert((row_id.clone(), col_id.clone()), d as u32);
		map.insert((col_id.clone(), row_id.clone()), d as u32);
            }
	}
    }

    Ok(map)
}

fn process_single_output_dir(
    dir: &PathBuf,
    monomers_db: &MonomersDB,
    cfg: &NerpaMS_Config,
    compound_distances: &HashMap<(String, String), u32>,
) -> Vec<FindTargetResults> {
    let input_items: Vec<InputItem> = {
	get_input_from_nerpa_results(
	    dir, None, None
	)
	    .expect(format!("Failed to load Nerpa results from {}",
			    dir))
    };

    find_siblings(&input_items, monomers_db, cfg, compound_distances)
}

fn write_results(results: &[FindTargetResults], out: &PathBuf) -> Result<()> {
    let tsv = FindTargetResults::to_tsv(results);
    std::fs::write(out, tsv)?;
    Ok(())
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let nerpa_output_dirs: Vec<PathBuf> = {
	let root = cli.nerpa_results.clone();

	if root.join("nerpa.log").is_file() {
            vec![root]
	} else {
            fs::read_dir(&root)?
		.filter_map(|e| e.ok())
		.map(|e| e.path())
		.filter(|p| p.is_dir())
		.filter(|p| p.join("nerpa.log").is_file())
		.collect()
	}
    };

    let compound_distances: HashMap<(String, String), u32> = {
	get_compound_distances(&cli.compound_distances)
	    .expect(format!("Failed to load compound distances from {}",
			    cli.compound_distances.display()).as_str())
    };

    let cfg = NerpaMS_Config {};
    let monomers_db = load_monomers_db(&cli.monomers_db_json);

    let all_results: Vec<FindTargetResults> = {
	nerpa_output_dirs
	    .iter()
	    .flat_map(|dir| process_single_output_dir(dir, &monomers_db, &cfg, &compound_distances))
	    .collect()
    };
	
    write_results(&all_results, &cli.out)?;

    Ok(())
}
				     
