from typing import Dict, Optional
from pathlib import Path
from src.antismash_parsing.antismash_name_mappings import KNOWN_SUBSTRATES
from src.data_types import AA34, LogProb, Prob
from src.monomer_names_helper import MonomerResidue, MonomerNamesHelper, UNKNOWN_RESIDUE
from src.generic.functional import timing_decorator
import pandas as pd
from math import log
from Bio import SeqIO
import json

from src.pipeline.logger import NerpaLogger

PARAS_RESIDUE = str

# TODO: put in monomer_names_helper.py
def paras_residue_to_nerpa_residue(residue: PARAS_RESIDUE,
                                   monomer_names_helper: MonomerNamesHelper) -> MonomerResidue:
        def get_paras_name_core(paras_name: PARAS_RESIDUE) -> str:
            match paras_name:
                case '3-(2-nitrocyclopropylalanine)':
                    return 'alanine'
                case '3S-methylaspartic acid branched':
                    return '3S-methylaspartic acid'
                case _:
                    return paras_name.split('-')[-1]

        paras_name_core = get_paras_name_core(residue)
        as_short = next(substrate.short
                        for substrate in KNOWN_SUBSTRATES
                        if paras_name_core in substrate.long)
        result = monomer_names_helper.parsed_name(as_short, name_format='antismash').residue
        if result == UNKNOWN_RESIDUE:
            print(f"Unknown residue: {residue}")
        return result


def extract_paras_results_json(paras_results_json: Path) -> Dict[AA34, Dict[PARAS_RESIDUE, Prob]]:
    predictions_per_aa34 = {}
    for domain_predictions in json.load(paras_results_json.open()):
        aa34 = domain_predictions['extended_signature']
        predictions_per_aa34[aa34] = {prediction['substrate_name']: prediction['probability']
                                      for prediction in domain_predictions['predictions']}

    return predictions_per_aa34


def extract_paras_results(paras_predictions: Path,
                          extended_signatures: Path) -> Dict[AA34, Dict[PARAS_RESIDUE, Prob]]:
    # Step 1: Extract sequence ID to extended signature mapping using BioPython
    seq_id_to_signature = {
        record.id: AA34(str(record.seq))
        for record in SeqIO.parse(extended_signatures, "fasta")
    }

    # Step 2: Read the predictions file into a DataFrame
    predictions_df = pd.read_csv(paras_predictions, sep="\t")

    # Convert DataFrame to a dictionary (seq_id -> substrate -> probability)
    predictions = {
        row[0]: {row[i]: float(row[i + 1]) for i in range(1, len(row) - 1, 2)}
        for row in predictions_df.itertuples(index=False)
    }

    # Step 3: Merge the two dictionaries to map extended signatures to predictions
    return {
        seq_id_to_signature[seq_id]: predictions[seq_id]
        for seq_id in predictions if seq_id in seq_id_to_signature
    }


def paras_predictions_to_nerpa_predictions(paras_predictions: Dict[PARAS_RESIDUE, Prob],
                                           monomer_names_helper: MonomerNamesHelper) -> Dict[MonomerResidue, LogProb]:
    #print(paras_predictions)
    nerpa_predictions = {res: 0.0 for res in monomer_names_helper.supported_residues}

    for paras_residue, prob in paras_predictions.items():
        nerpa_res = paras_residue_to_nerpa_residue(paras_residue, monomer_names_helper)
        nerpa_predictions[nerpa_res] += prob

    for res in nerpa_predictions:
        if nerpa_predictions[res] > 1:
            print(f"Invalid probability: {nerpa_predictions[res]}")
            nerpa_predictions[res] = 1.0
        nerpa_predictions[res] = log(nerpa_predictions[res]) if nerpa_predictions[res] > 0 else float('-inf')

    return nerpa_predictions


@timing_decorator
def get_paras_results_all(paras_results: Path,
                          monomer_names_helper: MonomerNamesHelper,
                          log: NerpaLogger) -> Dict[AA34, Dict[MonomerResidue, LogProb]]:
    log.info("Extracting PARAS results...")

    if not paras_results.exists() or not paras_results.is_dir():
        log.error(f"Path {paras_results} does not exist or is not a directory.")
        return {}

    combined_results = {}

    for folder in paras_results.rglob("*"):
        if folder.is_dir():
            for results_json in filter(lambda f: f.suffix == ".json", folder.iterdir()):
                try:
                    paras_results = extract_paras_results_json(results_json)
                except:
                    log.error(f"Error extracting PARAS results from {results_json}")
                    continue
                combined_results.update(paras_results)

            fasta_file = folder / "run_extended_signatures.fasta"
            results_file = folder / "run_paras_results.txt"

            if fasta_file.exists() and results_file.exists():
                try:
                    paras_results = extract_paras_results(results_file, fasta_file)
                except:
                    log.error(f"Error extracting PARAS results from {results_file} and {fasta_file}")
                    continue
                combined_results.update(paras_results)

    parsed_results = {}
    for aa34, paras_predictions in combined_results.items():
        parsed_results[aa34] = paras_predictions_to_nerpa_predictions(paras_predictions, monomer_names_helper)

    return parsed_results