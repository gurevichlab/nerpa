from typing import Dict, Optional
from pathlib import Path
from src.data_types import AA34, LogProb
from src.monomer_names_helper import MonomerResidue, MonomerNamesHelper
import pandas as pd
from Bio import SeqIO

PARAS_RESIDUE = str

def paras_residue_to_nerpa_residue(residue: PARAS_RESIDUE,
                                   monomer_names_helper: MonomerNamesHelper) -> MonomerResidue:


def extract_paras_results(paras_predictions: Path,
                          extended_signatures: Path) -> Dict[AA34, Dict[MonomerResidue, LogProb]]:
    # Step 1: Extract sequence ID to extended signature mapping using BioPython
    seq_id_to_signature = {
        record.id: str(record.seq)
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


def get_paras_results_all(paras_results: Path,
                          monomer_names_helper: MonomerNamesHelper) -> Dict[AA34, Dict[MonomerResidue, LogProb]]:
    combined_results = {}

    for folder in paras_results.rglob("*"):
        if folder.is_dir():
            fasta_file = folder / "run_extended_signatures.fasta"
            results_file = folder / "run_paras_results.txt"

            if fasta_file.exists() and results_file.exists():
                parsed_results = extract_paras_results(results_file, fasta_file)
                combined_results.update(parsed_results)

    return combined_results

