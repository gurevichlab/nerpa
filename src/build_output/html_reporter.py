import json
import math
import os
import shutil
from pathlib import Path
from typing import List, Dict
from src.config import OutputConfig
from src.matching.matching_types_match import Match


def _create_match_dicts(matches: List[Match], debug_output: bool) -> List[Dict]:
    # TODO: make this function less hard-coded, maybe do all stuff in AlignmentStep.to_dict()? See matching_types_alignment_step.py
    def _postprocess_alignments_for_html(match_dict: Dict) -> Dict:
        for i, alignment in enumerate(match_dict['Alignments']):
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
            match_dict['Alignments'][i] = cleaned_alignment

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


def _apply_substitutions(template: str, substitutions: Dict[str, str]) -> str:
    for placeholder, value in substitutions.items():
        template = template.replace(placeholder, value)
    return template


def create_html_report(output_cfg: OutputConfig, matches: List[Match], debug_output: bool = False):
    current_dir = Path(__file__).resolve().parent
    template_main_report_path = current_dir / 'main_report_template.html'
    main_report_path = output_cfg.html_report

    with open(template_main_report_path, 'r') as file:
        main_report_html_template = file.read()

    # TODO: Maybe add these paths directly to config_paths?
    html_aux_dir = output_cfg.main_out_dir / 'html_aux'
    report_data_js_path = html_aux_dir / 'report_data.js'
    html_aux_dir.mkdir()
    match_dicts = _create_match_dicts(matches, debug_output)
    # the main (root) HTML report and associated JSON
    with open(report_data_js_path, 'w') as json_file:
        json_file.write('var data = ')
        json.dump(match_dicts, json_file, indent=4)

    path_substitutions = {
        '{{HTML_AUX_DIR}}': str(html_aux_dir.relative_to(output_cfg.main_out_dir)),
        '{{ANTISMASH_OUT_DIR}}': str(output_cfg.antismash_out_dir.relative_to(output_cfg.main_out_dir))
    }
    main_html_report = _apply_substitutions(main_report_html_template, path_substitutions)
    with open(main_report_path, 'w') as file:
        file.write(main_html_report)
    # copying logo to be embedded in the HTML report
    shutil.copy(output_cfg.logo, html_aux_dir)
