#!/usr/bin/env python3
"""
Extends signatures.tsv with additional columns containing SVM (NRPys) predictions
Every column corresponds to a single amino acid supported by SVM (see SVM_SUBSTRATES)
and contains a 4-tuple with SVM prediction scores from single aa to small cluster to large cluster to phys.-chem. class
"""

# # You should add the antiSMASH root directory to sys.path as shown below
# # (or configure this via PyCharm: right-click on the antismash subdir inside the antismash dir,
# #                                 Mark Directory as -> Sources Root)
# antismash_root = '/path/to/antismash/source/code'
# import sys
# sys.path.append(antismash_root)
import antismash.modules.nrps_pks.nrpys as antismash_nrpys
from collections import OrderedDict
from pathlib import Path
import nrpys

# TODO: make it options?
# Get the absolute path of the directory where the script is located
script_dir = Path(__file__).resolve().parent
# Define the path to the 'data' directory relative to the script directory
data_dir = script_dir.parent / "data"

signatures_fpath = data_dir / "signatures.tsv"
SVM_models_dir = data_dir / "models"

output_fpath = data_dir / "extended_signatures.tsv"
OUTPUT_SEPARATOR = "\t"
MAX_LINES_TO_PROCESS = -1  # for debug only, -1 to process everything


SVM_SUBSTRATES = ["Arg", "Asp", "Glu", "Asn", "Lys", "Gln", "Orn", "ohOrn", "Aad",
                  "Ala", "Gly", "Val", "Leu", "Ile", "Abu", "Iva", "Ser", "Thr", "Hpg", "dHpg", "Cys", "Pro", "Pip",
                  "Phe", "Tyr", "Trp", "2,3-dohBza", "Sal", "Pgl", "R-ohTyr"]
# Note on column naming:
# A domain has (experimentally) validated ability to attract "all substrates"
# but actual NRP structure(s) contains only "nrp_substrates"
HEADER = ["aa10", "aa34", "substrates", "num_domains"] + SVM_SUBSTRATES


class SignatureEntry:
    def __init__(self, raw_entry_line):
        assert len(raw_entry_line.split()) == 5, "Error: corrupted line in the signatures file"
        self.aa10, self.aa34, all_substrates, nrp_substrates, sources = raw_entry_line.split()
        self.substrates = all_substrates.split('|')  # TODO: make an option to easily switch between all/nrp _substrates
        self.num_sources = sources.count('|') + 1

    def get_content(self):
        return [self.aa10, self.aa34, self.substrates, self.num_sources]


def run_nrpys(signature_entries):
    # some magic configuration copied from antiSMASH
    config = nrpys.Config()
    config.model_dir = SVM_models_dir
    config.stachelhaus_signatures = signatures_fpath
    config.skip_v1 = True
    config.skip_v3 = True

    names: list[str] = [f"{idx}" for idx in range(len(signature_entries))]
    signatures: list[str] = [signature_entry.aa34 for signature_entry in signature_entries]

    return [antismash_nrpys.PredictorSVMResult.from_nrpys(nrpys_res)
            for nrpys_res in nrpys.run(config, names, signatures)]


def main():
    signature_entries = list()

    # TODO: get rid of SignatureEntry and make reading simply with Pandas
    with open(signatures_fpath) as f:
        for idx, raw_entry in enumerate(f):
            if MAX_LINES_TO_PROCESS > 0 and idx > MAX_LINES_TO_PROCESS:
                break
            signature_entries.append(SignatureEntry(raw_entry))

    svm_results = run_nrpys(signature_entries)

    with open(output_fpath, 'w') as f:
        f.write(OUTPUT_SEPARATOR.join(HEADER) + "\n")
        for signature_entry, svm_result in zip(signature_entries, svm_results):
            svm_prediction_scores = OrderedDict((short_name, [0.0] * 4) for short_name in SVM_SUBSTRATES)
            for idx, svm_prediction in enumerate([svm_result.single_amino, svm_result.small_cluster,
                                                  svm_result.large_cluster, svm_result.physicochemical_class]):
                for substrate in svm_prediction.substrates:
                    svm_prediction_scores[substrate.short][idx] = svm_prediction.score
            f.write(OUTPUT_SEPARATOR.join(map(str, signature_entry.get_content() +
                                                   list(map(tuple, svm_prediction_scores.values()))))
                    + '\n')
    print(f"Done! See results in {output_fpath}")


if __name__ == "__main__":
    main()
