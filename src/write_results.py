from typing import (
    Callable,
    Iterable,
    List,
    TypeVar,
    Tuple,
    Union
)
from src.matching.alignment_types import Match
from src.data_types import BGC_Variant, NRP_Variant
from src.rban_parsing.rban_parser import Parsed_rBAN_Record
from src.build_output.html_reporter import create_html_report
from src.nerpa_ms.monomer_graph.draw_graph import draw_molecule, draw_monomer_graph
from src.nerpa_ms.monomer_graph.monomer_graph import MonomerGraph
from src.config import ConfigPaths

from pathlib import Path
from io import StringIO
import csv
import yaml
import json
from itertools import groupby

T = TypeVar('T')
U = TypeVar('U')

def sort_groupby(items: Iterable[T],
                 key: Callable[[T], U],
                 reverse: bool=False) -> Iterable[Tuple[U, Iterable[T]]]:
    return groupby(sorted(items, key=key, reverse=reverse), key=key)


def write_yaml(data, out_file: Path):
    # dirty hack to erase information about types and make output less verbose
    # https://github.com/yaml/pyyaml/issues/408
    yaml.emitter.Emitter.prepare_tag = lambda self, tag: ''

    # another hack (albeit less dirty) to forbid yaml creating references
    # https://stackoverflow.com/questions/13518819/avoid-references-in-pyyaml
    yaml.Dumper.ignore_aliases = lambda *args: True

    with open(out_file, 'w') as out:
        yaml.dump(data, out,
                  default_flow_style=None, sort_keys=False)


def build_report(matches: List[Match]) -> str:
    result = StringIO()
    csv_writer = csv.DictWriter(result,
                                fieldnames=('Score', 'NRP_ID', 'NRP_Variant_Idx', 'Genome_ID', 'BGC_ID', 'BGC_Variant_Idx'),
                                delimiter='\t')
    csv_writer.writeheader()
    csv_writer.writerows({'Score': match.normalized_score,
                          'NRP_ID': match.nrp_variant.nrp_id,
                          'NRP_Variant_Idx': match.nrp_variant.variant_idx,
                          'Genome_ID': match.bgc_variant.genome_id,
                          'BGC_ID': match.bgc_variant.bgc_idx,
                          'BGC_Variant_Idx': match.bgc_variant.variant_idx}
                         for match in matches)
    return result.getvalue()


def write_matches_per_id(matches: List[T],
                         output_dir: Path,
                         get_id: Callable[[T], str]):
    for id_, id_matches in sort_groupby(matches, get_id):  # python sort is stable so groups will be sorted by score
        (output_dir / Path(id_)).write_text('\n\n'.join(map(str, id_matches)))


def write_nrp_variants(nrp_variants: List[NRP_Variant],
                       config_paths: ConfigPaths,
                       rban_records: Union[List[Parsed_rBAN_Record], None] = None,
                       draw: bool = True):
    if rban_records is not None:
        write_yaml([rban_record.to_compact_dict() for rban_record in rban_records], config_paths.rban_graphs)
        for rban_record in rban_records:
            monomer_graph = MonomerGraph.from_rban_record(rban_record)
            if draw:
                draw_molecule(monomer_graph, config_paths.nrp_images_dir / Path(f'molecules/{rban_record.compound_id}.png'))
                draw_monomer_graph(monomer_graph,
                                   output_file=config_paths.nrp_images_dir / Path(f'graphs/{rban_record.compound_id}.png'))

    config_paths.nrp_variants_dir.mkdir()
    for nrp_id, nrp_id_variants in sort_groupby(nrp_variants, key=lambda nrp_variant: nrp_variant.nrp_id):
        write_yaml(list(nrp_id_variants), config_paths.nrp_variants_dir / Path(f'{nrp_id}.yaml'))


def write_bgc_variants(bgc_variants: List[BGC_Variant],
                       output_dir: Path):
    for (genome_id, bgc_id), bgc_id_variants in sort_groupby(bgc_variants, key=lambda bgc_variant: (bgc_variant.genome_id, bgc_variant.bgc_idx)):
        write_yaml(list(bgc_id_variants), output_dir / f'{genome_id}_{bgc_id}.yaml')


def write_matches_details(matches: List[Match],
                          matches_details_output_dir: Path):
    matches_details_output_dir.mkdir()
    write_yaml([match.to_dict_light() for match in matches],
               matches_details_output_dir / 'matches.yaml')

    (matches_details_output_dir / Path('per_BGC')).mkdir()
    write_matches_per_id(matches, matches_details_output_dir / Path('per_BGC'),
                         get_id=lambda match: f'{match.bgc_variant.genome_id}_{match.bgc_variant.bgc_idx}')

    (matches_details_output_dir / Path('per_NRP')).mkdir()
    write_matches_per_id(matches, matches_details_output_dir / Path('per_NRP'),
                         get_id=lambda match: match.nrp_variant.nrp_id)


def write_results(matches: List[Match],
                  config_paths: ConfigPaths,
                  bgc_variants: Union[List[BGC_Variant], None] = None,
                  nrp_variants: Union[List[NRP_Variant], None] = None,
                  rban_records: Union[List[Parsed_rBAN_Record], None] = None,
                  matches_details: bool = True,
                  draw_molecules: bool = True,
                  html_report: bool = True):
    config_paths.report.write_text(build_report(matches))

    json_report_path = config_paths.report.with_suffix('.json')
    with open(json_report_path, 'w') as json_file:
        json.dump([match.to_dict_light() for match in matches], json_file, indent=4)

    if bgc_variants is not None:
        config_paths.bgc_variants_dir.mkdir()
        write_bgc_variants(bgc_variants, config_paths.bgc_variants_dir)
    if nrp_variants is not None:
        write_nrp_variants(nrp_variants, config_paths, rban_records, draw=draw_molecules)

    if matches_details:
        write_matches_details(matches, config_paths.matches_details)

    if html_report:
        create_html_report(config_paths, json_report_path)  # TODO?: pass the matches directly instead of reading from 'report.tsv'

