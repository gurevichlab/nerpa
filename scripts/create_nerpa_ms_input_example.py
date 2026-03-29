from pathlib import Path
import json
import yaml

nerpa_results = Path("/home/ilianolhin/work/tools/nerpa/nerpa_results/trash_single")
hmms_json_path = nerpa_results / 'intermediate_files' / 'hmms.json'

with open(hmms_json_path) as f:
    hmms = json.load(f)

hmm = hmms[0]['hmm_for_matching']

rban_records_yaml = nerpa_results / 'preprocessed_input' / 'parsed_rban_records.yaml'

with open(rban_records_yaml) as f:
    rban_records = yaml.safe_load(f)

rban_record = rban_records[0]

linearization = [7, 1, 6, 5, 4, 3, 2]


json_dict = {
    "hmm": hmm,
    "parsed_rban_record": rban_record,
    "linearization": linearization
}

output_json_path = Path('/home/ilianolhin/work/tools/nerpa/src/nerpa_ms/nerpa_ms_core/fixtures/input_example.json')

with open(output_json_path, 'w') as f:
    json.dump([json_dict], f, indent=4)

    
