# ARGS=""
import subprocess
import os
os.environ["RUST_BACKTRACE"] = "1"

subprocess.run(["cargo", "build", "--bin", "create_monomers_db"])
command = ' '.join(['cargo run --bin create_monomers_db --',
                    '--parsed-rban-records-json fixtures/rban_records_list_example.json',
                    '--out output/monomers_db.json'])
subprocess.run(command.split())

