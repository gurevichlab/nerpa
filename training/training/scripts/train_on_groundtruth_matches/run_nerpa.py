import pandas as pd
from pathlib import Path
import subprocess


def run_nerpa_on_all_pairs(matches_table: pd.DataFrame,
                           antismash_results: Path,
                           rban_results: Path,
                           nerpa_dir: Path,
                           output_dir: Path) -> Path:
    nrp_ids = matches_table['NRP variant'].unique()
    for i, nrp_id in enumerate(nrp_ids, start=1):
        bgc_id = nrp_id.split('.')[0]
        command = ' '.join([
            'python', str(nerpa_dir / 'nerpa.py'),
            '-a', str(antismash_results / bgc_id),
            '--rban-json', str(rban_results / f'{nrp_id}.json'),
            '-o', str(output_dir / nrp_id), '--force-existing-outdir'
        ])
        print(f'{i}/{len(nrp_ids)}: {nrp_id}')
        subprocess.run(command, shell=True, check=True)
    return output_dir

