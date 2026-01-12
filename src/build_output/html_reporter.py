import json
import math
import os
import shutil
from dataclasses import asdict
from pathlib import Path
from typing import List, Dict
from src.config import OutputConfig
from src.matching.match_type import Match
from src.antismash_parsing.bgc_variant_types import BGC_Variants_Info
from src.rban_parsing.nrp_variant_types import NRP_Variants_Info


def _create_match_dicts(matches: List[Match],
                        debug_output: bool,
                        default_score_field: str = 'log_odds_vs_avg_bgc') -> List[Dict]:
    # TODO: make this function less hard-coded, maybe do all stuff in AlignmentStep.to_dict()? See alignment_step_type.py
    def _postprocess_alignments_for_html(match_dict: Dict) -> Dict:
        for i, alignment in enumerate(match_dict['alignments']):
            cleaned_alignment = []
            for alignment_step in alignment:
                # omitting auxiliary HMM transitions unless explicilty asked to show everything (--debug)
                if alignment_step['A-domain_idx'] == alignment_step['rBAN_idx'] == '---' and not debug_output:
                    continue
                alignment_step['Modifying_domains'] = alignment_step['Modifying_domains'].replace('EPIMERIZATION', 'E')
                alignment_step['Modifying_domains'] = alignment_step['Modifying_domains'].replace('METHYLATION', 'MT')
                alignment_step['NRP_chirality'] = alignment_step['NRP_chirality'].replace('UNKNOWN', 'unk')
                alignment_step['NRP_modifications'] = 'MT' if alignment_step['NRP_methylated'] == 'True' else '---'  # possible values are 'True', 'False', '---'
                cleaned_alignment.append(alignment_step)
            # Update the original alignment in the match_dict
            match_dict['alignments'][i] = cleaned_alignment

        # Adding log-odds score and the default score column
        match_dict['score'] = match_dict[default_score_field]
        return match_dict

    # TODO: think if we need some minimal leading_zeros always, e.g., at least 3
    # Determine the number of digits needed based on the length of matches
    min_leading_zeros = 3
    leading_zeros = max(min_leading_zeros,
                        math.ceil(math.log10(len(matches))) if matches else 1)

    return [
        _postprocess_alignments_for_html({
            **match.to_dict(),
            'Match_ID': f'{i:0{leading_zeros}d}'
        })
        for i, match in enumerate(matches)
    ]


def _create_serializable_bgc_metadata(bgc_variants_info: BGC_Variants_Info):
    # Collect only the metadata we need for BGCs (JSON-serializable)
    bgc_metadata = []
    for v in bgc_variants_info.bgc_variants:
        meta = getattr(v, 'metadata', None)
        if not meta:
            continue
        bgc_vid = v.bgc_variant_id.to_dict()  # {'bgc_id': {...}, 'variant_idx': int}
        bgc_metadata.append({
            'bgc_id': bgc_vid['bgc_id'],  # keep as dict so JS can access structured fields
            'variant_idx': bgc_vid['variant_idx'],
            'metadata': v.metadata.to_dict()
        })

    return bgc_metadata


def _create_serializable_nrp_metadata(nrp_variants_info: NRP_Variants_Info):
    # Collect only the metadata we need for NRPs (JSON-serializable)
    nrp_metadata = []
    for v in nrp_variants_info.nrp_variants:
        meta = getattr(v, 'metadata', None)
        if not meta:
            continue
        nrp_metadata.append({
            'nrp_id': v.nrp_variant_id.nrp_id,
            'variant_idx': v.nrp_variant_id.variant_idx,
            'metadata': asdict(v.metadata)
        })

    return nrp_metadata


def _create_serializable_bgc_representatives(bgc_variants_info: BGC_Variants_Info):
    groups = {}
    for member_id, repr_id in bgc_variants_info.bgc_id_to_repr_id.items():
        groups.setdefault(repr_id, []).append(member_id)

    out = []
    for repr_id, member_ids in groups.items():
        # TODO: thnk of it, see below. Ensure representative itself is not included
        # member_ids = [m for m in member_ids if m != repr_id]
        out.append({
            "repr": repr_id.to_dict(),
            "members": [m.to_dict() for m in sorted(member_ids,
                                                    key=lambda x: (x.bgc_id.antiSMASH_file,
                                                                   x.bgc_id.contig_idx,
                                                                   x.bgc_id.bgc_idx,
                                                                   x.variant_idx))]
        })
    return out


def _create_serializable_nrp_representatives(nrp_variants_info: NRP_Variants_Info):
    """
    Reverse nrp_id_to_repr_id into JSON-serializable list of groups:
    [
        { "repr": {...}, "members": [{...}, ...] },
        ...
    ]
    Note: the representative itself is also a part of the members
    """
    groups = {}
    for member_id, repr_id in nrp_variants_info.nrp_id_to_repr_id.items():
        groups.setdefault(repr_id, []).append(member_id)

    out = []
    for repr_id, member_ids in groups.items():
        # TODO: think how it is better to do with or without representative
        #  Ensure representative itself is not included (robust to potential self-maps)
        # member_ids = [m for m in member_ids if m != repr_id]
        out.append({
            "repr": repr_id._asdict(),
            "members": [m._asdict() for m in sorted(
                member_ids, key=lambda x: (x.nrp_id, x.variant_idx)
            )]
        })
    return out


def _apply_substitutions(template: str, substitutions: Dict[str, str]) -> str:
    for placeholder, value in substitutions.items():
        template = template.replace(placeholder, value)
    return template


def create_html_report(output_cfg: OutputConfig,
                       matches: List[Match],
                       bgc_variants_info: BGC_Variants_Info,
                       nrp_variants_info: NRP_Variants_Info,
                       debug_output: bool = False,
                       default_score_field: str = 'log_odds_vs_avg_bgc'):
    current_dir = Path(__file__).resolve().parent
    template_main_report_path = current_dir / 'main_report_template.html'
    main_report_path = output_cfg.html_report

    with open(template_main_report_path, 'r') as file:
        main_report_html_template = file.read()

    # TODO: Maybe add these paths directly to config_paths?
    html_aux_dir = output_cfg.main_out_dir / 'html_aux'
    report_data_js_path = html_aux_dir / 'report_data.js'
    html_aux_dir.mkdir()
    match_dicts = _create_match_dicts(matches,
                                      debug_output,
                                      default_score_field=default_score_field)
    bgc_metadata = _create_serializable_bgc_metadata(bgc_variants_info)
    nrp_metadata = _create_serializable_nrp_metadata(nrp_variants_info)
    bgc_representatives = _create_serializable_bgc_representatives(bgc_variants_info)
    nrp_representatives = _create_serializable_nrp_representatives(nrp_variants_info)

    # the main (root) HTML report and associated JSON
    with open(report_data_js_path, 'w') as json_file:
        json_file.write('var data = ')
        json.dump(match_dicts, json_file, indent=4)
        json_file.write(';\n')

        json_file.write('var bgc_metadata = ')
        json.dump(bgc_metadata, json_file, indent=4)
        json_file.write(';\n')

        json_file.write('var nrp_metadata = ')
        json.dump(nrp_metadata, json_file, indent=4)
        json_file.write(';\n')

        json_file.write('var bgc_representatives = ')
        json.dump(bgc_representatives, json_file, indent=4)
        json_file.write(';\n')

        json_file.write('var nrp_representatives = ')
        json.dump(nrp_representatives, json_file, indent=4)
        json_file.write(';\n')

    path_substitutions = {
        '{{HTML_AUX_DIR}}': str(html_aux_dir.relative_to(output_cfg.main_out_dir)),
        '{{ANTISMASH_OUT_DIR}}': str(output_cfg.antismash_out_dir.relative_to(output_cfg.main_out_dir))
    }
    main_html_report = _apply_substitutions(main_report_html_template, path_substitutions)
    with open(main_report_path, 'w') as file:
        file.write(main_html_report)
    # copying logo to be embedded in the HTML report
    shutil.copy(output_cfg.logo, html_aux_dir)