import yaml
from pathlib import Path
from typing import List, Dict, Tuple
from collections import OrderedDict
from src.monomer_names_helper import MonomerResidue
from src.write_results import write_yaml
from training.specificity_prediction.scripts.training_table_generator.load_data import (
    AA10_code,
    AA34_code,
    SubstrateCodesDict
)
from training.specificity_prediction.scripts.training_table_generator.process_data import NerpaTableEntry


def write_monomer_signatures(monomer_signatures_yaml: Path,
                             substrate_to_aa10_codes: SubstrateCodesDict,
                             substrate_to_aa34_codes: SubstrateCodesDict,
                             nerpa_supported_monomers: List[MonomerResidue],
                             nerpa_unknown_monomer: MonomerResidue):
    # Reformatting the two internal dicts into one codes_dict with lexicographical order of monomers
    codes_dict: Dict[MonomerResidue, Tuple[List[AA10_code], List[AA34_code]]] = OrderedDict()

    for substrate in nerpa_supported_monomers + [nerpa_unknown_monomer]:
        aa10_list = list(substrate_to_aa10_codes[substrate].keys())
        aa34_list = list(substrate_to_aa34_codes[substrate].keys())
        codes_dict[substrate] = (aa10_list, aa34_list)

    # Dumping codes_dict to YAML file
    write_yaml(codes_dict, monomer_signatures_yaml)
    print(f"Saved monomer aa10 and aa34 codes in {monomer_signatures_yaml}")


def write_nerpa_table(output_file: Path,
                      nerpa_table_entries: List[NerpaTableEntry],
                      SEPARATOR: str):
    with open(output_file, 'w') as f:
        f.write(SEPARATOR.join(NerpaTableEntry.header()) + "\n")
        f.write('\n'.join(table_row.to_str(SEPARATOR)
                          for table_row in nerpa_table_entries))
    print(f'Saved the Nerpa training table in {output_file}')

