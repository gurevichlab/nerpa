from typing import (
    Dict,
    List,
    NamedTuple
)

from src.pipeline.command_line_args_helper import (
    CommandLineArgs
)
from src.pipeline.logging.logger import NerpaLogger
from src.config import Config
from src.antismash_parsing.bgc_variant_types import (
    BGC_Variant_ID
)
from src.monomer_names_helper import MonomerNamesHelper

from src.hmm.detailed_hmm import DetailedHMM

from src.rban_parsing.get_linearizations import NRP_Linearizations
from src.generic.other import json_round_floats, json_remove_infinities
from src.matching.hmm_match import HMM_Match
import json
from pathlib import Path
import subprocess
from dataclasses import dataclass


class CppOutput(NamedTuple):
    matches: List[HMM_Match]
    p_values_by_bgc_variant: Dict[BGC_Variant_ID, List[float]]

def reformat_json(json_str: str) -> str:
    # Force newlines for lists
    json_str = json_str.replace("[[", "[\n    [").replace("]]", "]\n]")
    json_str = json_str.replace("], [", "],\n    [")
    json_str = json_str.replace("}, {", "},\n    {")
    json_str = json_str.replace(': {"', ': {\n"')
    json_str = json_str.replace(', "', ',\n"')
    return json_str

@dataclass
class PipelineHelperCpp:
    config: Config
    args: CommandLineArgs
    log: NerpaLogger
    monomer_names_helper: MonomerNamesHelper

    def dump_hmms(self, detailed_hmms: List[DetailedHMM]) -> Path:
        out_file = self.config.output_config.cpp_io_config.hmms_json
        out_file.parent.mkdir(parents=True, exist_ok=True)
        data = [{"bgc_variant_info": detailed_hmm.bgc_variant.bgc_variant_id.to_dict(),
                 "hmm_for_matching": detailed_hmm.to_hmm(unknown_chirality_allowed=True,
                                                         emission_weights_type='LogProb').to_json(),
                 "hmm_for_p_values_estimation": detailed_hmm.to_hmm(unknown_chirality_allowed=False,
                                                                    emission_weights_type='LogProb').to_json()}
                for detailed_hmm in detailed_hmms]

        #pretty_json = reformat_json(json.dumps(data))
        refined_data = json_round_floats(data, ndigits=3)
        refined_data = json_remove_infinities(refined_data, infinity=1e30)
        pretty_json = json.dumps(refined_data)
        #pretty_json = json.dumps(data)
        with open(out_file, 'w') as f:
            f.write(pretty_json)
        return out_file


    def dump_nrp_linearizations(self,
                                nrp_linearizations: List[NRP_Linearizations],
                                any_hmm: DetailedHMM) -> Path:
        out_file = self.config.output_config.cpp_io_config.nrp_linearizations_json
        out_file.parent.mkdir(parents=True, exist_ok=True)
        data = [nrp_linearization.to_mon_codes_json(self.monomer_names_helper, any_hmm)
                for nrp_linearization in nrp_linearizations]

        #pretty_json = reformat_json(json.dumps(data))
        pretty_json = json.dumps(data)
        with open(out_file, 'w') as f:
            f.write(pretty_json)
        return out_file

    def run_cpp_matcher(self,
                        hmms_json: Path,
                        nrp_linearizations_json: Path) -> CppOutput:
        # Prepare the command
        cmd = list(map(str,
                          [
                              self.config.cpp_matcher_exec,  # Path to compiled C++ executable
                              '--hmms_json', hmms_json,
                              '--nrps_json', nrp_linearizations_json,
                              '--max_num_matches_per_bgc', str(self.config.matching_config.max_num_matches_per_bgc),
                              '--max_num_matches_per_nrp', str(self.config.matching_config.max_num_matches_per_nrp),
                              '--min_num_matches_per_bgc', str(self.config.matching_config.min_num_matches_per_bgc),
                              '--min_num_matches_per_nrp', str(self.config.matching_config.min_num_matches_per_nrp),
                              '--max_num_matches', str(self.config.matching_config.max_num_matches),
                              '--threads', str(self.args.threads),
                              '--output', self.config.output_config.cpp_io_config.cpp_output_json
                          ]))
        self.log.info(f"Running C++ matcher with command: {' '.join(cmd)}")
        # q: print current time

        # Run the C++ executable
        run_result = subprocess.run(cmd, capture_output=True, text=True)
        if run_result.returncode != 0:
            self.log.error(f"C++ matcher failed with return code {run_result.returncode}")
            if run_result.stdout.strip():
                self.log.error(f"Stdout: {run_result.stdout}")
            if run_result.stderr.strip():
                self.log.error(f"Stderr: {run_result.stderr}")
            raise RuntimeError("C++ matcher execution failed")
        else:
            self.log.info(f"Stdout: {run_result.stdout}")
            if run_result.stderr.strip():
                self.log.error(f"Stderr: {run_result.stderr}")
            self.log.info("C++ matcher completed successfully")

        # Collect the results
        with open(self.config.output_config.cpp_io_config.cpp_output_json, 'r') as f:
            json_data = json.load(f)
            hmm_matches = [HMM_Match.from_json(match_dict)
                           for match_dict in json_data["matches"]]
            p_values_by_bgc_variant = {BGC_Variant_ID.from_dict(p_values_dict["bgc_variant_info"]):
                                           p_values_dict["p_values"]
                                        for p_values_dict in json_data["p_values"]}
        return CppOutput(matches=hmm_matches,
                         p_values_by_bgc_variant=p_values_by_bgc_variant)

    def get_hmm_matches_and_p_values(self,
                                     detailed_hmms: List[DetailedHMM],
                                     nrp_linearizations: List[NRP_Linearizations]) -> CppOutput:
        if not detailed_hmms or not nrp_linearizations:
            return CppOutput(matches=[], p_values_by_bgc_variant={})
        
        hmms_json = self.dump_hmms(detailed_hmms)
        nrp_linearizations_json = self.dump_nrp_linearizations(nrp_linearizations, detailed_hmms[0])
        return self.run_cpp_matcher(hmms_json, nrp_linearizations_json)

