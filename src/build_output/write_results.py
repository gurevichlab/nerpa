from typing import (
    Callable,
    List,
    Optional, Set
)

from src.hmm.detailed_hmm import DetailedHMM
from src.matching.hmm_match import HMM_Match
from src.matching.match_type import Match
from src.config import OutputConfig
from src.antismash_parsing.bgc_variant_types import BGC_Variants_Info, BGC_Variant_ID
from src.monomer_names_helper import MonomerNamesHelper
from src.rban_parsing.nrp_variant_types import NRP_Variants_Info, NRP_Variant_ID
from src.generic.combinatorics import sort_groupby
from src.rban_parsing.rban_parser import Parsed_rBAN_Record
from src.build_output.html_reporter import create_html_report
from src.nerpa_ms.monomer_graph.draw_graph import draw_molecule, draw_monomer_graph
from src.nerpa_ms.monomer_graph.monomer_graph import MonomerGraph
from src.pipeline.logging.logger import NerpaLogger
from pathlib import Path
from io import StringIO
import csv
import yaml


def write_yaml(data, out_file: Path):
    # dirty hack to erase information about types and make output less verbose
    # https://github.com/yaml/pyyaml/issues/408
    yaml.emitter.Emitter.prepare_tag = lambda self, tag: ''

    # another hack (albeit less dirty) to forbid yaml creating references
    # https://stackoverflow.com/questions/13518819/avoid-references-in-pyyaml
    yaml.Dumper.ignore_aliases = lambda *args: True

    out_file.parent.mkdir(parents=True, exist_ok=True)
    with open(out_file, 'w') as out:
        yaml.dump(data, out,
                  default_flow_style=None, sort_keys=False)


def draw_hmms_with_optimal_paths(hmms: List[DetailedHMM],
                                 hmm_matches: List[HMM_Match],
                                 main_out_dir: Path):
    bgc_variant_to_hmm = {detailed_hmm.bgc_variant.bgc_variant_id: detailed_hmm
                          for detailed_hmm in hmms}
    for hmm_match in hmm_matches:
        detailed_hmm = bgc_variant_to_hmm[hmm_match.bgc_variant_id]
        for i, opt_path in enumerate(hmm_match.optimal_paths):
            fname = f'{hmm_match.bgc_variant_id.bgc_id.to_str_short()}_{hmm_match.nrp_id}_al{i}.png'
            detailed_hmm.draw(main_out_dir / 'HMMs_with_paths' / fname,
                              opt_path)


def build_report(matches: List[Match]) -> str:
    result = StringIO()
    csv_writer = csv.DictWriter(result,
                                fieldnames=('LogOdds_vs_avg_NRP', 'LogOdds_vs_avg_BGC', 'Raw_score', 'p_value',
                                            'NRP_ID', 'NRP_Variant_Idx',
                                            'Genome_ID', 'antiSMASH_input', 'BGC_idx', 'BGC_Variant_Idx'),
                                delimiter='\t')
    csv_writer.writeheader()
    csv_writer.writerows({'Score': match.log_odds_vs_avg_bgc,
                          'NRP_ID': match.nrp_variant_id.nrp_id,
                          'NRP_Variant_Idx': match.nrp_variant_id.variant_idx,
                          'antiSMASH_input': match.bgc_variant_id.bgc_id.antiSMASH_file,
                          'BGC_idx': match.bgc_variant_id.bgc_id.bgc_idx,
                          'Contid_idx': match.bgc_variant_id.bgc_id.contig_idx,
                          'BGC_Variant_Idx': match.bgc_variant_id.variant_idx}
                         for match in matches)
    return result.getvalue()


def write_matches_per_id(matches: List[Match],
                         output_dir: Path,
                         get_id: Callable[[Match], str]):
    for id_, id_matches in sort_groupby(matches, get_id):  # python sort is stable so groups will be sorted by score
        (output_dir / Path(id_)).write_text('\n\n'.join(map(str, id_matches)))


def write_nrp_variants(nrp_variants_info: NRP_Variants_Info,
                       nrp_ids_to_write: Set[NRP_Variant_ID],
                       output_cfg: OutputConfig,
                       rban_records: Optional[List[Parsed_rBAN_Record]] = None,
                       log: Optional[NerpaLogger] = None,
                       monomer_names_helper: Optional[MonomerNamesHelper] = None,):
    if rban_records is None:
        rban_records = []

    compound_ids_to_write = {nrp_id.nrp_id for nrp_id in nrp_ids_to_write}

    write_yaml([
        nrp_variant.to_dict()
        for nrp_variant in nrp_variants_info.nrp_variants
        if nrp_variant.nrp_variant_id in nrp_ids_to_write
    ],
        output_cfg.nrp_variants)

    write_yaml([
        {
            'nrp_id': nrp_id,
            'representative_id': repr_id
        }
        for nrp_id, repr_id in nrp_variants_info.nrp_id_to_repr_id.items()
        if repr_id in nrp_ids_to_write
    ],
        output_cfg.nrp_representatives)

    if rban_records:
        write_yaml([rban_record.to_compact_dict()
                    for rban_record in rban_records
                    if rban_record.compound_id in compound_ids_to_write],
                   output_cfg.rban_graphs)
        write_yaml([rban_record.to_dict(monomer_names_helper=monomer_names_helper)
                    for rban_record in rban_records
                    if rban_record.compound_id in compound_ids_to_write],
                   output_cfg.parsed_rban_records)
        if output_cfg.draw_molecules:
            for rban_record in filter(lambda r: r.compound_id in compound_ids_to_write,
                                      rban_records):
                #print(f'Drawing {rban_record.compound_id}', flush=True)
                monomer_graph = MonomerGraph.from_rban_record(rban_record)
                draw_monomer_graph(monomer_graph,
                                   with_rban_indexes=True,
                                   output_path=output_cfg.nrp_images_dir / f'graphs/{rban_record.compound_id}.svg',
                                   monomer_names_helper=monomer_names_helper)
                try:
                    draw_molecule(monomer_graph,
                                  rban_indexes=True,
                                  monomer_labels=True,
                                  output_file=output_cfg.nrp_images_dir / f'molecules/{rban_record.compound_id}.svg',
                                  monomer_names_helper=monomer_names_helper)
                except Exception as e:
                    if log is not None:
                        log.info(f'Failed to draw molecule for {rban_record.compound_id}: {e}')
                    raise
    else:
        if log is not None:
            log.info('rBAN records not provided, skipping writing rBAN graphs and drawing molecules')


def write_bgc_variants(bgc_variants_info: BGC_Variants_Info,
                       bgc_ids_to_write: Set[BGC_Variant_ID],
                       cfg: OutputConfig):
    write_yaml([
        bgc_variant.to_dict()
        for bgc_variant in sorted(bgc_variants_info.bgc_variants,
                                  key=lambda bgc_id_variant: bgc_id_variant.bgc_variant_id)
        if bgc_variant.bgc_variant_id in bgc_ids_to_write
    ],
        cfg.bgc_variants)

    write_yaml([
        {
            'bgc_id': bgc_id.to_dict(),
            'representative_id': repr_id.to_dict()
        }
        for bgc_id, repr_id in bgc_variants_info.bgc_id_to_repr_id.items()
        if repr_id in bgc_ids_to_write
    ],
        cfg.bgc_representatives)

def write_matches_details(matches: List[Match],
                          matches_details_output_dir: Path):
    matches_details_output_dir.mkdir()
    write_yaml([match.to_dict() for match in matches],
               matches_details_output_dir / 'matches.yaml')

    (matches_details_output_dir / Path('per_BGC')).mkdir()
    write_matches_per_id(matches, matches_details_output_dir / Path('per_BGC'),
                         get_id=lambda match: f'{match.bgc_variant_id.bgc_id.to_str_short()}')

    (matches_details_output_dir / Path('per_NRP')).mkdir()
    write_matches_per_id(matches, matches_details_output_dir / Path('per_NRP'),
                         get_id=lambda match: match.nrp_variant_id.nrp_id)


def write_results(matches: List[Match],
                  bgc_variants_info: BGC_Variants_Info,
                  nrp_variants_info: NRP_Variants_Info,
                  output_cfg: OutputConfig,
                  matches_details: bool = True,
                  html_report: bool = True,
                  debug_output: bool = False,
                  write_only_what_is_matched: bool = True,
                  log: Optional[NerpaLogger] = None,
                  monomer_names_helper: Optional[MonomerNamesHelper] = None,):
    output_cfg.report.write_text(build_report(matches))

    if write_only_what_is_matched:
        bgc_ids_to_write = {match.bgc_variant_id for match in matches}
        nrp_ids_to_write = {match.nrp_variant_id for match in matches}
    else:
        bgc_ids_to_write = bgc_variants_info.bgc_id_to_repr_id.keys()
        nrp_ids_to_write = nrp_variants_info.nrp_id_to_repr_id.keys()

    output_cfg.bgc_variants.parent.mkdir(exist_ok=True, parents=True)
    write_bgc_variants(bgc_variants_info,
                       bgc_ids_to_write,
                       output_cfg)

    output_cfg.nrp_variants.parent.mkdir(exist_ok=True, parents=True)
    write_nrp_variants(nrp_variants_info,
                       nrp_ids_to_write,
                       rban_records=nrp_variants_info.rban_records,
                       output_cfg=output_cfg,
                       log=log,
                       monomer_names_helper=monomer_names_helper)

    if matches_details:
        write_matches_details(matches, output_cfg.matches_details)

    if html_report:
        create_html_report(output_cfg, matches, bgc_variants_info, nrp_variants_info, debug_output)
