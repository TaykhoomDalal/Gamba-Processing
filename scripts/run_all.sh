#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
args=()
if [[ "${1:-}" == "--skip-phylop" ]]; then
  args+=(--skip-phylop)
elif [[ $# -gt 0 ]]; then
  echo "usage: $0 [--skip-phylop]" >&2
  exit 2
fi

python "$root/scripts/download_data.py" "${args[@]}"
python "$root/Functional-Regions/process.py" "${args[@]}"
python "$root/Functional-Regions/process.py" \
  --include-noncoding \
  --regions-dir "$root/Functional-Regions/regions-noncoding-added" \
  "${args[@]}"
python "$root/ATG/process.py" "${args[@]}"
python "$root/scripts/verify_outputs.py" "${args[@]}"
