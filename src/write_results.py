from typing import (
    Callable,
    Iterable,
    List,
    TypeVar,
    Tuple,
    Optional,
    Union
)
from src.matching.match_type import Match
from src.config import OutputConfig
from src.data_types import BGC_Variant, NRP_Variant
from src.generic.combinatorics import sort_groupby
from src.rban_parsing.rban_parser import Parsed_rBAN_Record
from src.build_output.html_reporter import create_html_report
from src.nerpa_ms.monomer_graph.draw_graph import draw_molecule, draw_monomer_graph
from src.nerpa_ms.monomer_graph.monomer_graph import MonomerGraph
from src.pipeline.logger import NerpaLogger
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


def build_report(matches: List[Match]) -> str:
    result = StringIO()
    csv_writer = csv.DictWriter(result,
                                fieldnames=('Score', 'NRP_ID', 'NRP_Variant_Idx', 'Genome_ID', 'BGC_ID', 'BGC_Variant_Idx'),
                                delimiter='\t')
    csv_writer.writeheader()
    csv_writer.writerows({'Score': match.score,
                          'NRP_ID': match.nrp_variant_id.nrp_id,
                          'NRP_Variant_Idx': match.nrp_variant_id.variant_idx,
                          'Genome_ID': match.bgc_variant_id.bgc_id.genome_id,
                          'BGC_ID': match.bgc_variant_id.bgc_id.bgc_idx,
                          'BGC_Variant_Idx': match.bgc_variant_id.variant_idx}
                         for match in matches)
    return result.getvalue()


def write_matches_per_id(matches: List[Match],
                         output_dir: Path,
                         get_id: Callable[[Match], str]):
    for id_, id_matches in sort_groupby(matches, get_id):  # python sort is stable so groups will be sorted by score
        (output_dir / Path(id_)).write_text('\n\n'.join(map(str, id_matches)))


def write_nrp_variants(nrp_variants: List[NRP_Variant],
                       matches: List[Match],
                       output_cfg: OutputConfig,
                       rban_records: Optional[List[Parsed_rBAN_Record]] = None,
                       log: Optional[NerpaLogger] = None):
    matched_nrp_ids = {match.nrp_variant_id.nrp_id for match in matches}
    if rban_records is not None:
        write_yaml([rban_record.to_compact_dict() for rban_record in rban_records],
                   output_cfg.main_out_dir / Path('rban_graphs.yaml'))
        if output_cfg.draw_molecules:
            for rban_record in rban_records:
                if rban_record.compound_id not in matched_nrp_ids:
                    continue
                monomer_graph = MonomerGraph.from_rban_record(rban_record)
                draw_monomer_graph(monomer_graph,
                                   output_file=output_cfg.nrp_images_dir / f'graphs/{rban_record.compound_id}.png')
                try:
                    draw_molecule(monomer_graph, output_cfg.nrp_images_dir / f'molecules/{rban_record.compound_id}.png')
                except Exception as e:
                    if log is not None:
                        log.info(f'Failed to draw molecule for {rban_record.compound_id}: {e}')

    output_cfg.nrp_variants_dir.mkdir()
    for nrp_id, nrp_id_variants in sort_groupby(nrp_variants, key=lambda nrp_variant: nrp_variant.nrp_variant_id.nrp_id):
        write_yaml(list(nrp_id_variants), output_cfg.nrp_variants_dir / f'{nrp_id}.yaml')


def write_bgc_variants(bgc_variants: List[BGC_Variant],
                       output_dir: Path):
    for bgc_id, bgc_id_variants in sort_groupby(bgc_variants, key=lambda bgc_variant: bgc_variant.bgc_variant_id.bgc_id):
        bgc_id_str = f'{bgc_id.genome_id}_{bgc_id.contig_idx}_{bgc_id.bgc_idx}'
        bgc_id_variants = list(bgc_id_variants)
        write_yaml([bgc_variant.to_dict() for bgc_variant in bgc_id_variants],
                   output_dir / f'{bgc_id_str}.yaml')


def write_matches_details(matches: List[Match],
                          matches_details_output_dir: Path):
    matches_details_output_dir.mkdir()
    write_yaml([match.to_dict() for match in matches],
               matches_details_output_dir / 'matches.yaml')

    (matches_details_output_dir / Path('per_BGC')).mkdir()
    write_matches_per_id(matches, matches_details_output_dir / Path('per_BGC'),
                         get_id=lambda match: f'{match.bgc_variant_id.bgc_id.genome_id}_{match.bgc_variant_id.get_antismash_id()}')

    (matches_details_output_dir / Path('per_NRP')).mkdir()
    write_matches_per_id(matches, matches_details_output_dir / Path('per_NRP'),
                         get_id=lambda match: match.nrp_variant_id.nrp_id)


def write_results(matches: List[Match],
                  output_cfg: OutputConfig,
                  bgc_variants: Optional[List[BGC_Variant]] = None,
                  nrp_variants: Optional[List[NRP_Variant]] = None,
                  rban_records: Optional[List[Parsed_rBAN_Record]] = None,
                  matches_details: bool = True,
                  html_report: bool = True,
                  debug_output: bool = False,
                  log: Optional[NerpaLogger] = None):
    output_cfg.report.write_text(build_report(matches))

    if bgc_variants is not None:
        output_cfg.bgc_variants_dir.mkdir()
        write_bgc_variants(bgc_variants, output_cfg.bgc_variants_dir)
    if nrp_variants is not None:
        write_nrp_variants(nrp_variants,
                           rban_records=rban_records,
                           matches=matches,
                           output_cfg=output_cfg,
                           log=log)

    if matches_details:
        write_matches_details(matches, output_cfg.matches_details)

    if html_report:
        create_html_report(output_cfg, matches, debug_output)

