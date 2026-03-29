#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]

mod command_line_args;
mod data_types;
mod matcher;
mod viterbi_algorithm;

use crate::data_types::{MatchInfo, MatchingConfig, NRP_Linearizations_Info, HMM};
use clap::Parser;

fn main() {
    let args = command_line_args::Args::parse();
    let cfg = MatchingConfig::load_from_args(&args);

    let hmms = HMM::load_hmms(&args.hmms_json);
    let nrps_info = NRP_Linearizations_Info::load_nrps_linearizations(&args.nrps_json);

    let matches = matcher::get_matches(&hmms, &nrps_info, &cfg, args.threads);

    MatchInfo::dump_json_list(&matches, &args.output);
}
