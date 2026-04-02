# ARGS=""
import subprocess
import os
os.environ["RUST_BACKTRACE"] = "1"

subprocess.run(["cargo", "build"])
command = ' '.join(['cargo run --bin nerpa_ms_core --',
                    '--input fixtures/input_example.json',
                    '--max-edits 3',
                    '--num-variants-per-num-edits 5',
                    '--out output',
                    '--monomers-db-json data/monomers_db.json'])
subprocess.run(command.split())

