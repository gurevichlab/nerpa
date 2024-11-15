import json
import math
import os
import shutil
from pathlib import Path
from typing import List, Dict
from src.config import ConfigPaths
from src.matching.matching_types_match import Match


def _create_match_dicts(matches: List[Match]) -> List[Dict]:
    # TODO: think if we need some minimal leading_zeros always, e.g., at least 3
    # Determine the number of digits needed based on the length of matches
    min_leading_zeros = 3
    leading_zeros = max(min_leading_zeros,
                        math.ceil(math.log10(len(matches))) if matches else 1)

    return [
        {
            **match.to_dict(),
            'Match_ID': f'{i:0{leading_zeros}d}'
        }
        for i, match in enumerate(matches)
    ]


def _apply_substitutions(template: str, substitutions: Dict[str, str]) -> str:
    for placeholder, value in substitutions.items():
        template = template.replace(placeholder, value)
    return template


def create_html_report(config_paths: ConfigPaths, matches: List[Match]):
    current_dir = Path(__file__).resolve().parent
    template_main_report_path = current_dir / 'main_report_template.html'
    template_match_report_path = current_dir / 'match_report_template.html'
    main_report_path = config_paths.html_report

    with open(template_main_report_path, 'r') as file:
        main_report_html_template = file.read()

    # TODO/FIXME: It is better (but more difficult) to load everything in just one main HTML
    #  and skip the creation of zillions per match JSONs and HTML reports
    with open(template_match_report_path, 'r') as file:
        match_report_html_template = file.read()

    # TODO: Maybe add these paths directly to config_paths?
    html_aux_dir = config_paths.main_out_dir / 'html_aux'
    report_data_js_path = html_aux_dir / 'report_data.js'
    html_aux_dir.mkdir()
    match_dicts = _create_match_dicts(matches)
    # the main (root) HTML report and associated JSON
    with open(report_data_js_path, 'w') as json_file:
        json_file.write('var data = ')
        json.dump(match_dicts, json_file, indent=4)

    path_substitutions = {
        '{{HTML_AUX_DIR}}': str(html_aux_dir.relative_to(config_paths.main_out_dir))
    }
    main_html_report = _apply_substitutions(main_report_html_template, path_substitutions)
    with open(main_report_path, 'w') as file:
        file.write(main_html_report)
    # copying logo to be embedded in the HTML report
    shutil.copy(config_paths.logo, html_aux_dir)

    # per-match detailed HTML reports and associated JSONs
    # TODO?: limit the number of generate reports to the first (best) X matches? e.g., 1000
    for match_dict in match_dicts:
        match_data_js_path = html_aux_dir / f'match_{match_dict["Match_ID"]}.js'
        match_html_report_path = html_aux_dir / f'match_{match_dict["Match_ID"]}.html'
        with open(match_data_js_path, 'w') as json_file:
            json_file.write('var data = ')
            json.dump(match_dict, json_file, indent=4)

        path_substitutions = {
            # Note: the first one could be done with Path's .relative_to()
            # but it will fail for the rest two -- it can't properly form the "../" like output
            '{{DATA_JSON_PATH}}': os.path.relpath(match_data_js_path, start=html_aux_dir),
            '{{MAIN_REPORT_PATH}}': os.path.relpath(main_report_path, start=html_aux_dir),
            '{{ANTISMASH_OUT_DIR}}': os.path.relpath(config_paths.antismash_out_dir, start=html_aux_dir),
            '{{NRP_IMAGES_DIR}}': os.path.relpath(config_paths.nrp_images_dir, start=html_aux_dir)
        }
        match_html_report = _apply_substitutions(match_report_html_template, path_substitutions)
        with open(match_html_report_path, 'w') as file:
            file.write(match_html_report)
