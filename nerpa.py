import sys
from typing import List
from src.pipeline.pipeline_helper import PipelineHelper
from src.pipeline.logger import NerpaLogger
from src.rban_parsing.rban_parser import Parsed_rBAN_Record
from src.monomer_names_helper import MonomerNamesHelper, UNKNOWN_RESIDUE
from pathlib import Path
import yaml


def compute_compound_similarity(rban_records: List[Parsed_rBAN_Record]):
    print('Computing compound similarity')
    rban_graphs = [(rban_record.compound_id, rban_record.to_nx_monomer_graph())
                   for rban_record in rban_records]
    isomorphism_table = {}
    for i, (compound_id1, graph1) in enumerate(rban_graphs):
        print(f'{i}/{len(rban_graphs)}')
        isomorphism_table[compound_id1] = []
        for compound_id2, graph2 in rban_graphs:
            if compound_id1 == compound_id2:
                continue
            if Parsed_rBAN_Record.graphs_isomorphic(graph1, graph2):
                isomorphism_table[compound_id1].append(compound_id2)
    with open('isomorphism_table.yaml', 'w') as f:
        yaml.dump(isomorphism_table, f)
    return isomorphism_table


def compute_compound_info_table(rban_records: List[Parsed_rBAN_Record],
                                monomer_helper: MonomerNamesHelper):
    # q: create pandas table with 3 columns: id, num of recognized monomers, num of nerpa2-supported monomers
    # q: how to count recognized monomers? recognized monomers are those which rban name doesn't start with 'X' or have ':' inside
    # q: how to count nerpa2-supported monomers? nerpa2-supported which residues are not UNKNOWN_RESIDUE
    print('Computing compound info table')
    table = []
    for rban_record in rban_records:
        recognized_monomers = 0
        nerpa2_supported_monomers = 0
        for monomer in rban_record.monomers.values():
            rban_name = monomer.name
            if not rban_name.startswith('X') and ':' not in rban_name:
                recognized_monomers += 1
            if monomer_helper.parsed_name(rban_name, name_format='norine').residue != UNKNOWN_RESIDUE:
                nerpa2_supported_monomers += 1
        table.append({'id': rban_record.compound_id,
                      'recognized_monomers': recognized_monomers,
                      'nerpa2_supported_monomers': nerpa2_supported_monomers})
    with open('compound_info_table.tsv', 'w') as f:
        f.write('id\trecognized_monomers\tnerpa2_supported_monomers\n')
        for row in table:
            f.write(f"{row['id']}\t{row['recognized_monomers']}\t{row['nerpa2_supported_monomers']}\n")


def main(log: NerpaLogger):  # log is passed as an argument to make it easier to write log in case of exception
    pipeline_helper = PipelineHelper(log)

    bgc_variants = pipeline_helper.pipeline_helper_antismash.get_bgc_variants()
    hmms = pipeline_helper.construct_hmms(bgc_variants)
    #for hmm in hmms:
    #    hmm.draw(Path(f'{hmm.bgc_variant.genome_id}.png'))

    nrp_variants, rban_records = pipeline_helper.get_nrp_variants_and_rban_records()
    # compute_compound_similarity(rban_records)
    # compute_compound_info_table(rban_records, pipeline_helper.monomer_names_helper)
    # exit(0)
    nrp_linearizations = pipeline_helper.get_nrp_linearizations(nrp_variants)

    if pipeline_helper.args.only_preprocessing:
        matches = []
    else:
        matches = pipeline_helper.get_matches(hmms, nrp_linearizations)
    pipeline_helper.write_results(matches, bgc_variants, nrp_variants, rban_records)


if __name__ == "__main__":
    log = NerpaLogger()
    try:
        main(log)
    except Exception as e:
        _, exc_value, _ = sys.exc_info()
        log.exception(exc_value)
    finally:
        # TODO: clean up: remove all intermediate files
        pass
