import copy
from typing import Dict, List, Optional, Tuple, Set
from src.monomer_names_helper import (
    antiSMASH_MonomerName,
    MonomerResidue,
    MonomerNamesHelper,
    NorineMonomerName,
    UNKNOWN_RESIDUE
)
from src.write_results import write_yaml
from collections import Counter, defaultdict
import pandas as pd
import sys
from pathlib import Path
import yaml


def load_supported_monomers(monomer_frequencies_tsv: Path,
                            num_monomers: Optional[int],
                            SEPARATOR: str = '\t') -> List[MonomerResidue]:
    df = pd.read_csv(monomer_frequencies_tsv, sep=SEPARATOR)
    selected_monomers = df['core'].head(num_monomers) if num_monomers is not None else df['core']
    return sorted(selected_monomers)


def remove_unsupported_monomers(monomer_names_table: pd.DataFrame,
                                supported_monomers: List[MonomerResidue],
                                UNKNOWN_MONOMER: str) -> pd.DataFrame:
    new_monomer_names_table = copy.deepcopy(monomer_names_table)
    new_monomer_names_table['core'] = new_monomer_names_table['core'].apply(lambda x: x if x in supported_monomers else UNKNOWN_MONOMER)
    return new_monomer_names_table


AA10_code = str
AA34_code = str


def load_known_codes(signatures: pd.DataFrame) -> Tuple[Dict[antiSMASH_MonomerName, Set[AA10_code]], \
                                                        Dict[antiSMASH_MonomerName, Set[AA34_code]]]:
    signatures.columns = ['aa10', 'aa34', 'substrates_loose', 'substrates_strict', 'a_domains']

    # Split 'substrates_loose' into multiple rows, repeating 'aa10' and 'aa34' as necessary
    expanded = (signatures
                .assign(substrates_loose=signatures['substrates_loose'].str.split('|'))
                .explode('substrates_loose'))

    # Group by 'substrates_loose' and aggregate the 'aa10' and 'aa34' values into sets
    aa10_codes = expanded.groupby('substrates_loose')['aa10'].apply(set).to_dict()
    aa34_codes = expanded.groupby('substrates_loose')['aa34'].apply(set).to_dict()

    return aa10_codes, aa34_codes


def get_known_codes_only_cores(known_codes: Tuple[Dict[antiSMASH_MonomerName, Set[AA10_code]], \
                                                  Dict[antiSMASH_MonomerName, Set[AA34_code]]],
                               monomer_names_helper: MonomerNamesHelper) \
        -> Tuple[Dict[MonomerResidue, Set[AA10_code]],
                 Dict[MonomerResidue, Set[AA34_code]]]:
    aa10_codes, aa34_codes = known_codes
    aa10_codes_cores, aa34_codes_cores = defaultdict(set), defaultdict(set)
    for monomer_name, aa10_set in aa10_codes.items():
        aa34_set = aa34_codes[monomer_name]
        parsed_name = monomer_names_helper.parsed_name(monomer_name, 'antismash')
        aa10_codes_cores[parsed_name.residue].update(aa10_set)
        aa34_codes_cores[parsed_name.residue].update(aa34_set)
    return aa10_codes_cores, aa34_codes_cores


def get_core_frequencies(norine_graphs_table: pd.DataFrame,
                         monomer_names_helper: MonomerNamesHelper,
                         supported_monomers: List[MonomerResidue]) -> Dict[MonomerResidue, float]:
    norine_graphs_table.columns = ['id', 'graph', '_']

    def extract_monomer_names(graph: str) -> List[MonomerResidue]:
        return [monomer_names_helper.parsed_name(monomer_name, 'norine').residue
                for monomer_name in graph.split('@')[0].split(',')]

    counter = Counter(norine_graphs_table['graph']
                      .apply(extract_monomer_names)
                      .explode())

    for res in supported_monomers:
        if res not in counter:
            counter[res] = 1  # dirty fix: add monomers missing in norine graphs to the counter

    total = sum(counter.values())
    return {monomer: count / total for monomer, count in counter.items()}


def dump_known_codes_only_cores(aa10_codes: Dict[MonomerResidue, Set[AA10_code]],
                                aa34_codes: Dict[MonomerResidue, Set[AA34_code]],
                                output_file: Path):
    aa10_codes_lists = {monomer: list(codes) for monomer, codes in aa10_codes.items()}
    aa34_codes_lists = {monomer: list(codes) for monomer, codes in aa34_codes.items()}
    write_yaml({'aa10': aa10_codes_lists, 'aa34': aa34_codes_lists}, output_file)


def dump_monomer_names_table(monomer_names_table: pd.DataFrame, output_file: Path):
    monomer_names_table.to_csv(output_file, sep='\t', index=False)


SVM_SUBSTRATES = ["Arg", "Asp", "Glu", "Asn", "Lys", "Gln", "Orn", "ohOrn", "Aad", "Ala",
                  "Gly", "Val", "Leu", "Ile", "Abu", "Iva", "Ser", "Thr", "Hpg", "dHpg",
                  "Cys", "Pro", "Pip", "Phe", "Tyr", "Trp", "2,3-dohBza", "Sal", "Pgl", "R-ohTyr"]  # TODO: move to config

def update_config(nerpa_dir: Path,
                  core_frequencies: Dict[MonomerResidue, float],
                  monomer_names_helper: MonomerNamesHelper,
                  UNKNOWN_RESIDUE: MonomerResidue):
    config = yaml.safe_load((nerpa_dir / 'data' / 'old_config.yaml').read_text())
    config['specificity_prediction_config']['apriori_residue_prob'] = core_frequencies

    svm_cores = list(filter(lambda res: res != UNKNOWN_RESIDUE,
                            (monomer_names_helper.parsed_name(monomer_name, 'antismash').residue
                             for monomer_name in SVM_SUBSTRATES)))
    config['antismash_parsing_config']['SVM_SUBSTRATES'] = list(set(svm_cores))
    write_yaml(config, nerpa_dir / 'configs' / 'config.yaml')


def main():
    '''
    if len(sys.argv) != 2:
        print("Usage: python set_num_supported_monomers.py <num_supported_monomers>")
        sys.exit(1)

    num_supported_monomers = int(sys.argv[1])
    '''
    num_supported_monomers = 49

    nerpa_dir = Path(__file__).parent.parent.parent.resolve()
    supported_monomers = load_supported_monomers(nerpa_dir / 'data' / 'core_frequency.tsv', num_supported_monomers)

    monomer_names_table = pd.read_csv(nerpa_dir / 'data' / 'monomers_unique.tsv', sep='\t')
    monomer_names_table = remove_unsupported_monomers(monomer_names_table, supported_monomers, UNKNOWN_RESIDUE)
    monomer_names_helper = MonomerNamesHelper(monomer_names_table)

    known_codes = load_known_codes(pd.read_csv(nerpa_dir / 'data' / 'signatures.tsv',
                                               sep='\t', header=None))
    known_codes_only_cores = get_known_codes_only_cores(known_codes, monomer_names_helper)

    core_frequencies = get_core_frequencies(pd.read_csv(nerpa_dir / 'data' / 'norine_graphs.csv', sep=';'),
                                            monomer_names_helper,
                                            supported_monomers)

    dump_known_codes_only_cores(*known_codes_only_cores, nerpa_dir / 'configs' / 'aa_codes.yaml')
    dump_monomer_names_table(monomer_names_table, nerpa_dir / 'configs' / 'monomer_names_table.tsv')

    update_config(nerpa_dir, core_frequencies, monomer_names_helper, UNKNOWN_RESIDUE)


if __name__ == "__main__":
    main()