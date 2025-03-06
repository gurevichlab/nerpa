from typing import List
from pathlib import Path
from src.training.extract_data_for_training import DataForTraining
from src.training.hmm_infer_emission_params import BGC_Module_Info
import yaml

def load_paras_predictions(paras_results_yaml: Path) -> List[BGC_Module_Info]:
    paras_predictions = []
    data = yaml.safe_load(paras_results_yaml.read_text())
    return [BGC_Module_Info.from_dict(bgc_info)
            for bgc_info in data]


def substitute_predictions_with_paras(data_for_training: DataForTraining) -> DataForTraining:
    paras_yaml = Path('/home/ilianolhin/git/nerpa2/paras/paras34_predictions.yaml')
    paras_predictions = load_paras_predictions(paras_yaml)

    new_match_emissions = []
    for i, (bgc_info, bgc_module, nrp_monomer) in enumerate(data_for_training.match_emissions):
        try:
            paras_prediction = next(paras_prediction
                                    for paras_prediction in paras_predictions
                                    if all([paras_prediction.genome_id == bgc_info.genome_id,
                                            paras_prediction.gene_id == bgc_module.gene_id,
                                            paras_prediction.a_domain_idx == bgc_module.a_domain_idx])
                                           or paras_prediction.aa34 == bgc_module.aa34_code)
        except StopIteration:
            print(f'No prediction for {bgc_info.genome_id} {bgc_module.gene_id} {bgc_module.a_domain_idx}')
            continue
        bgc_module.residue_score = paras_prediction.residue_score
        new_match_emissions.append((bgc_info, bgc_module, nrp_monomer))

    return data_for_training

