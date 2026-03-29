# ARGS=""
import subprocess
import os
os.environ["RUST_BACKTRACE"] = "1"

subprocess.run(["cargo", "build"])
command = ' '.join(['cargo run --',
                    '--input fixtures/input_example.json',
                    '--max-edits 3',
                    '--num-variants-per-num-edits 5',
                    '--out output/example_output.json'])
subprocess.run(command.split())

