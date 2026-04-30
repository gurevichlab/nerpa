import json
from pathlib import Path
import yaml
import sys

inp = Path(sys.argv[1])
out = Path(sys.argv[2]) if len(sys.argv) > 2 else inp.with_suffix(".json")
print(f"Converting {inp} to {out}...")

data = yaml.safe_load(inp.read_text())

out.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
