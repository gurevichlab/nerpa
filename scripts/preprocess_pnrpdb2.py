import sys

import pandas as pd
import json
import subprocess
from pathlib import Path
from typing import Dict

import yaml


def run_nerpa_results_on_pnrpdb2(nerpa_dir: Path,
                                 pnrpdb2_path: Path,
                                 output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    nerpa_script = nerpa_dir / "nerpa.py"

    # Construct the command
    command = [
        "python3", str(nerpa_script),
        "--antismash", str(nerpa_dir / "test_data" / "antismash"),
        "--smiles-tsv", str(pnrpdb2_path),
        "--col-id", "ID",
        "--output-dir", str(output_dir),
        "--force-output-dir",
        "--max-num-matches-per-nrp", '1',
        "--min-num-matches-per-nrp", '1',
        "--max-num-matches", '0',
        "--skip-molecule-drawing",
        "--threads", "8",
        "--fast-matching",
        "--dump-all-preprocessed",
        "--keep-intermediate-files",
    ]

    # Execute the command and capture output
    print('Executing command:', ' '.join(command))
    try:
        result = subprocess.run(command, check=True)
        print(f"Nerpa run completed successfully")
    except subprocess.CalledProcessError as e:
        print(f"Error running Nerpa:\n{e.stderr}")
        raise

def load_nrp_to_representative(nerpa_results: Path) -> Dict[str, str]:
    import yaml
    nrp_representatives_path = (nerpa_results
                                / 'preprocessed_input'
                                / 'NRP_representatives.yaml')
    with nrp_representatives_path.open() as f:
        nrp_representatives = yaml.safe_load(f)
    nrp_id_to_repr_id = {entry['nrp_id'][0]: entry['representative_id'][0]
                         for entry in nrp_representatives}
    return nrp_id_to_repr_id


def main():
    nerpa_dir = Path(__file__).resolve().parent.parent
    pnrpdb2_path = (nerpa_dir
                    / 'data'
                    / 'input'
                    / 'pnrpdb2.tsv')
    pnrpdb2_df = pd.read_csv(pnrpdb2_path, sep='\t')

    nerpa_results = (nerpa_dir
                     / 'nerpa_results'
                     / 'preprocess_pnrpdb2')

    if '--use-existing-results' not in sys.argv:  # for debugging purposes
        run_nerpa_results_on_pnrpdb2(nerpa_dir=nerpa_dir,
                                     pnrpdb2_path=pnrpdb2_path,
                                     output_dir=nerpa_results)

    nrp_id_to_repr_id = load_nrp_to_representative(nerpa_results)
    nrp_variants = yaml.safe_load(nerpa_results / 'preprocessed_input' / 'NRP_variants.yaml')
    rban_records = json.loads((nerpa_dir / 'intermediate_files' / 'rban.output.json').read_text())

    def is_mibig_norine(nrp_id: str) -> bool:
        return nrp_id.startswith('BGC') or nrp_id.startswith('NOR')

    def is_representative(nrp_id: str) -> bool:
        return nrp_id in nrp_id_to_repr_id.values()

    # 1. pnrpdb2_mibig_norine.tsv
    pnrpdb2_mibig_norine_path = (nerpa_dir
                                 / 'data'
                                 / 'input'
                                 / 'pnrpdb2_mibig_norine.tsv')
    is_mibig_norine_col = pnrpdb2_df['ID'].apply(is_mibig_norine)
    df_mibig_norine = pnrpdb2_df[is_mibig_norine_col]
    df_mibig_norine.to_csv(pnrpdb2_mibig_norine_path, sep='\t', index=False)

    # 2. pnrpdb2_deduplicated.tsv
    pnrpdb2_deduplicated_path = (nerpa_dir
                                 / 'data'
                                 / 'input'
                                 / 'pnrpdb2_deduplicated.tsv')
    is_representative_col = pnrpdb2_df['ID'].apply(is_representative)
    df_representatives = pnrpdb2_df[is_representative_col]
    df_representatives.to_csv(pnrpdb2_deduplicated_path, sep='\t', index=False)

    # 3. pnrpdb2_mibig_norine_deduplicated.tsv
    pnrpdb2_mibig_norine_deduplicated_path = (nerpa_dir
                                              / 'data'
                                              / 'input'
                                              / 'pnrpdb2_mibig_norine_deduplicated.tsv')
    df_mibig_norine_deduplicated = pnrpdb2_df[is_mibig_norine_col & is_representative_col]
    df_mibig_norine_deduplicated.to_csv(pnrpdb2_mibig_norine_deduplicated_path, sep='\t', index=False)

    # 4. Dump preprocessed NRP variants and rBAN records
    preprocessed_dir = nerpa_dir / 'data' / 'input' / 'preprocessed'
    preprocessed_dir.mkdir(exist_ok=True, parents=True)
    for table_name in (
            'pnrpdb2',
            'pnrpdb2_mibig_norine',
            'pnrpdb2_deduplicated',
            'pnrpdb2_mibig_norine_deduplicated',
    ):
        table = nerpa_dir / 'data' / 'input' / f'{table_name}.tsv'
        table_df = pd.read_csv(table, sep='\t')
        nrp_ids_in_table = set(table_df['ID'])

        nrp_variants_out_path = preprocessed_dir / f'{table_name}_nrp_variants.yaml'
        rban_records_out_path = preprocessed_dir / f'{table_name}_rban_records.json'

        nrp_variants_subset = [
            nrp_variant for nrp_variant in nrp_variants
            if nrp_variant['nrp_variant_id']['nrp_id'] in nrp_ids_in_table
        ]
        rban_records_subset = [
            rban_record for rban_record in rban_records
            if rban_record['id'] in nrp_ids_in_table
        ]

        with open(nrp_variants_out_path, 'w') as f:
            yaml.safe_dump(nrp_variants_subset, f)

        with open(rban_records_out_path, 'w') as f:
            json.dump(rban_records_subset, f, indent=2)


if __name__ == '__main__':
    main()
