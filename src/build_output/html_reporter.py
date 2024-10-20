import json
import csv
from pathlib import Path
from typing import List
from src.config import ConfigPaths
from src.matching.alignment_types import Match


def create_html_report(config_paths: ConfigPaths, matches: List[Match]):
    current_dir = Path(__file__).resolve().parent
    template_file = current_dir / 'report_template.html'
    html_report_file = config_paths.html_report

    # Read the template HTML
    with open(template_file, 'r') as file:
        html_template = file.read()

    # TODO: Maybe add these paths directly to config_paths?
    html_aux_dir = config_paths.main_out_dir / 'html_aux'
    report_data_js_path = html_aux_dir / 'report_data.js'
    html_aux_dir.mkdir()
    with open(report_data_js_path, 'w') as json_file:
        json_file.write('var data = ')
        json.dump([match.to_dict_light() for match in matches], json_file, indent=4)

    final_html = html_template.replace('{{REPORT_DATA_JS_PATH}}',
                                       str(report_data_js_path.relative_to(config_paths.main_out_dir)))

    # Write the final HTML file with data
    with open(html_report_file, 'w') as file:
        file.write(final_html)
