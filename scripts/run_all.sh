#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python "$root/scripts/download_data.py"
python "$root/Functional-Regions/process.py"
python "$root/Functional-Regions/process.py" \
  --include-noncoding \
  --regions-dir "$root/Functional-Regions/regions-noncoding-added"
python "$root/ATG/process.py"
python "$root/scripts/verify_outputs.py"
