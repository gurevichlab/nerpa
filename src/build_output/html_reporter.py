import csv
from pathlib import Path
from src.config import ConfigPaths


def create_html_report(config_paths: ConfigPaths):
    current_dir = Path(__file__).resolve().parent
    template_file = current_dir / 'report_template.html'
    tsv_report_file = config_paths.report
    html_report_file = config_paths.html_report

    # Read the template HTML
    with open(template_file, 'r') as file:
        html_template = file.read()

    # FIXME: don't read from TSV and generate HTML table manually
    # Better generate JSON (report.json) and read/generate the table from it directly in JavaScript
    # (in report_template.html)

    # Read the TSV file and generate table headers and rows
    with open(tsv_report_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        headers = next(reader)
        headers_html = ''.join([f'<th>{header}</th>' for header in headers])
        rows_html = ''.join(
            ['<tr>' + ''.join([f'<td>{cell}</td>' for cell in row]) + '</tr>' for row in reader]
        )

    # Replace placeholders with actual headers and rows
    final_html = html_template.replace('{{table_headers}}', headers_html).replace('{{table_rows}}', rows_html)

    # Write the final HTML file with data
    with open(html_report_file, 'w') as file:
        file.write(final_html)
