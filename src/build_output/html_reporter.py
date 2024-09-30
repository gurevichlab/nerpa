import csv
from pathlib import Path


def create_html_report(output_dir: Path):
    current_dir = Path(__file__).resolve().parent
    template_file = current_dir / 'report_template.html'
    # FIXME: use filenames from the config instead of hard coding!
    tsv_report_file = output_dir / 'report.tsv'
    html_report_file = output_dir / 'report.html'

    # Read the template HTML
    with open(template_file, 'r') as file:
        html_template = file.read()

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
