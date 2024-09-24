import argparse
from pathlib import Path
from typing import List
from src.monomer_names_helper import MonomerResidue
from training.specificity_prediction.scripts.training_table_generator.load_data import load_data
from training.specificity_prediction.scripts.training_table_generator.process_data import generate_nerpa_table_entries
from training.specificity_prediction.scripts.training_table_generator.write_data import write_monomer_signatures, write_nerpa_table

SEPARATOR = "\t"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Process monomer signature data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Input group
    input_group = parser.add_argument_group("input", "Input options")
    input_group.add_argument("--num-monomers", type=int, default=None,
                             help="Number of most frequent monomers to consider. By default all monomers are included.")
    input_group.add_argument("--unknown-monomer-name", type=str, default="unknown",
                             help="Special name designating all rare (unsupported) monomers.")
    input_group.add_argument("--as-monomer-frequencies", type=str, default="data/core_frequency.tsv",
                             help="Path to a TSV file with all antismash core monomers and their frequencies sorted from the most common to the rarest.")
    input_group.add_argument("--norine-monomer-frequencies", type=str, default="data/norine_residues_freqs.tsv",
                             help="Path to a TSV file with all norine core monomers and their frequencies sorted from the most common to the rarest.")
    input_group.add_argument("--extended-signatures", type=str, default="data/extended_signatures.tsv",
                             help="Path to the extended signatures file.")
    input_group.add_argument("--monomers-table", type=str, default="data/monomers_unique.tsv",
                             help="Path to the monomer names conversion table.")

    # Output group
    output_group = parser.add_argument_group("output", "Output options")
    output_group.add_argument("--monomer-signatures-dict", type=str, default="/home/ilianolhin/git/nerpa2/training/specificity_prediction/output/monomer_signatures.yaml",
                              help="Path to a YAML file to store the dictionary with all aa10 and aa34 signatures per every supported core monomer.")
    output_group.add_argument("--nerpa-training-table", type=str, default="/home/ilianolhin/git/nerpa2/training/specificity_prediction/output/nerpa_scoring_table.tsv",
                              help="Path to a TSV file to store the main specificity training table.")

    return parser.parse_args()


def ensure_output_directories_exist(output_paths: List[Path]):
    for path in output_paths:
        # Convert the file path to a Path object
        parent_dir = Path(path).parent

        # Create the directory (and any missing parents) if it doesn't exist
        if not parent_dir.exists():
            parent_dir.mkdir(parents=True, exist_ok=True)
            print(f"Created directory: {parent_dir}")
        # else:
        #     print(f"Directory already exists: {parent_dir}")


def main():
    # Parse arguments
    args = parse_args()

    # reading
    input_data = load_data(as_monomer_frequencies_tsv=args.as_monomer_frequencies,
                           norine_monomer_frequencies_tsv=args.norine_monomer_frequencies,
                           monomers_table_tsv=args.monomers_table,
                           extended_signatures_tsv=args.extended_signatures,
                           num_monomers=args.num_monomers,
                           nerpa_unknown_monomer=MonomerResidue(args.unknown_monomer_name),
                           SEPARATOR=SEPARATOR)

    # processing
    nerpa_table_entries = generate_nerpa_table_entries(input_data)

    # writing
    ensure_output_directories_exist([args.monomer_signatures_dict, args.nerpa_training_table])
    write_monomer_signatures(args.monomer_signatures_dict,
                             input_data.substrate_to_aa10_codes,
                             input_data.substrate_to_aa34_codes,
                             input_data.nerpa_supported_monomers,
                             input_data.nerpa_unknown_monomer)
    write_nerpa_table(args.nerpa_training_table, nerpa_table_entries, SEPARATOR)


if __name__ == "__main__":
    main()
