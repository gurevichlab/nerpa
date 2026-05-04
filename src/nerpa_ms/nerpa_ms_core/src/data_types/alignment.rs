use std::collections::HashMap;

use crate::data_types::{bgc_variant::BGC_Variant, common_types::{LogProb, MonomerIdx}, hmm::{HMM, StateIdx}, parsed_rban_record::Parsed_rBAN_Record};

pub struct AlignmentStep {
    pub module_idx: Option<usize>,
    pub monomer_idx: Option<MonomerIdx>,
}

pub struct Alignment <'a>{
    pub steps: Vec<AlignmentStep>,
    pub bgc_variant: &'a BGC_Variant,
    pub rban_record: &'a Parsed_rBAN_Record,
}

use std::io::Write;
use tabwriter::TabWriter;

impl<'a> Alignment<'a> {
    pub fn new(
	hmm_path: &[StateIdx],
	linearization: &[MonomerIdx],
	bgc_variant: &'a BGC_Variant,
	hmm: &HMM,
	rban_record: &'a Parsed_rBAN_Record,
    ) -> Self {
	let emitting_states_on_path: Vec<StateIdx> = {
	    hmm_path.iter()
		.filter(|&&state_idx| hmm.is_emitting(state_idx))
	    .cloned()
	    .collect()
	};
	assert_eq!(emitting_states_on_path.len(),
		   linearization.len(),
		   "Number of emitting states on HMM path must match length of linearization");

	let state_idx_to_matched_module_idx: HashMap<StateIdx, usize> = {
	    hmm.module_match_states
		.iter()
		.enumerate()
		.map(|(module_idx, &state_idx)| (state_idx, module_idx))
		.collect()
	};

	let mut last_matched_module_idx: isize = -1;
	let mut steps: Vec<AlignmentStep> = Vec::new();
	for (&state_idx, &monomer_idx) in (emitting_states_on_path.iter()
					 .zip(linearization.iter())) {
	    match state_idx_to_matched_module_idx.get(&state_idx) {
		Some(&module_idx) => {
		    for skipped_module_idx in (last_matched_module_idx + 1)..(module_idx as isize) {
			steps.push(AlignmentStep {
			    module_idx: Some(skipped_module_idx as usize),
			    monomer_idx: None,
			});
		    }

		    steps.push(AlignmentStep {
			module_idx: Some(module_idx),
			monomer_idx: Some(monomer_idx),
		    });
		    last_matched_module_idx = module_idx as isize;
		},
		None => {
		    steps.push(AlignmentStep {
			module_idx: None,
			monomer_idx: Some(monomer_idx),
		    });
		},
	    }
	}

	Alignment {
	    steps,
	    bgc_variant,
	    rban_record,
	}
	
    }

    pub fn to_tsv_string(&self) -> String {
        use std::cmp::Ordering;
        use std::fmt::Write;

        fn logprob_to_percent(logp: &LogProb) -> u64 {
            // LogProb is assumed to be a natural-log probability (ln p).
            // Convert to percent in 0..100: exp(logp) * 100, rounded.
            let percent = ((*logp).exp() * 100.0).round();
            percent.clamp(0.0, 100.0) as u64
        }

        let mut out = String::new();

        out.push_str(
            "GeneID\tModuleIdx\tTop3Predictions\tModuleMT\tMonomerResidue\tMonomerName\tMonomerMT\tMonomerIdx\n",
        );

        for step in &self.steps {
            let module_opt = step
                .module_idx
                .and_then(|idx| self.bgc_variant.modules.get(idx));

            let monomer_opt = step
                .monomer_idx
                .and_then(|idx| self.rban_record.monomers.get(&idx));

            let gene_id = module_opt.map(|m| m.gene_id.as_str()).unwrap_or("-");

            let module_idx_str = step
                .module_idx
                .map(|idx| idx.to_string())
                .unwrap_or_else(|| "-".to_string());

            let top3_predictions = if let Some(module) = module_opt {
                let mut pairs: Vec<_> = module.residue_score.iter().collect();
                pairs.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap_or(Ordering::Equal));

                let mut s = String::new();
                for (i, (residue, score)) in pairs.into_iter().take(3).enumerate() {
                    if i > 0 {
                        s.push_str(",");
                    }
                    let pct = logprob_to_percent(score);
		    let residue_str = if residue.is_unknown() {
			"unk".to_string()
		    } else {
			residue.0.clone()
		    };
                    let _ = write!(s, "{}({})", residue_str, pct);
                }

                if s.is_empty() { "-".to_string() } else { s }
            } else {
                "-".to_string()
            };

            let module_mt = if let Some(module) = module_opt {
                if module.modifications.iter().any(|m| m.contains("METHYLATION")) {
                    "MT"
                } else {
                    "-"
                }
            } else {
                "-"
            };

            let (monomer_residue, monomer_name, monomer_mt, monomer_idx_str) =
                if let Some(monomer) = monomer_opt {
                    let residue_str = if monomer.nerpa_core.is_unknown() {
			"unk".to_string()
		    } else { monomer.nerpa_core.0.clone() };

                    let monomer_residue = if let Some(module) = module_opt {
                        module
                            .residue_score
                            .get(&monomer.nerpa_core)
                            .map(|lp| format!("{}({})", residue_str, logprob_to_percent(lp)))
                            .unwrap_or_else(|| "-".to_string())
                    } else {
                        "-".to_string()
                    };

                    let monomer_name = monomer.name.0.clone();
                    let monomer_mt = if monomer.methylated { "MT" } else { "-" };

                    let monomer_idx_str = step
                        .monomer_idx
                        .map(|idx| idx.0.to_string())
                        .unwrap_or_else(|| "-".to_string());

                    (monomer_residue, monomer_name, monomer_mt, monomer_idx_str)
                } else {
                    let monomer_idx_str = step
                        .monomer_idx
                        .map(|idx| format!("{idx:?}"))
                        .unwrap_or_else(|| "-".to_string());

                    ("-".to_string(), "-".to_string(), "-", monomer_idx_str)
                };

            let _ = writeln!(
                out,
                "{gene_id}\t{module_idx_str}\t{top3_predictions}\t{module_mt}\t{monomer_residue}\t{monomer_name}\t{monomer_mt}\t{monomer_idx_str}"
            );
        }

        out
    }

    pub fn to_tsv_string_aligned(&self) -> String {
	let tsv = self.to_tsv_string();
	let mut tw = TabWriter::new(Vec::new());
	let _ = tw.write_all(tsv.as_bytes());
	let _ = tw.flush();

	String::from_utf8(tw.into_inner().unwrap()).unwrap()
    }
	
}
