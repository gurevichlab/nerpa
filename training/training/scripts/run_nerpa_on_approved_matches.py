from pathlib import Path
import yaml
import datetime
import subprocess

approved_matches_file = Path('/home/ilianolhin/git/nerpa2/training/training/matches_inspection_results/approved_matches.yaml')

current_time = str(datetime.datetime.now()).split('.')[0].replace('-', '_').replace(' ', '_').replace(':', '_')
#output_dir = Path(f'/home/ilianolhin/git/nerpa2/training/training/nerpa_results/new/{current_time}')
output_dir = Path(f'/home/ilianolhin/git/nerpa2/training/training/nerpa_results/new/temporary')

antismash_results_dir = Path('/home/ilianolhin/programs/antismash7/mibig_results')
rban_results_dir = Path('/home/ilianolhin/git/nerpa2/training/training/rban_records')

for match in yaml.safe_load(approved_matches_file.read_text())[:1]:
    nrp_id = match['NRP']
    bgc_id = nrp_id.split('.')[0]
    command = ' '.join([
        'python', '/home/ilianolhin/git/nerpa2/nerpa.py',
        '-a', str(antismash_results_dir / bgc_id),
        '--rban-json', str(rban_results_dir / f'{nrp_id}.json'),
        '-o', str(output_dir / nrp_id), '--force-existing-outdir'
    ])
    subprocess.run(command, shell=True, check=True)