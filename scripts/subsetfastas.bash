#!/bin/bash
set -euo pipefail

usage() {
  echo "Usage: $(basename "$0") MIN MAX file1.fasta [file2.fasta ...]"
  echo "Example: $(basename "$0") 800 20000 sample1.fasta sample2.fasta"
  echo ".fasta, .fa, .fas, .faa, .fna, are accepted extensions"
  exit 1
}

# --- args ---
[[ $# -lt 3 ]] && usage
MIN="$1"; shift
MAX="$1"; shift

# basic validation
[[ ! "$MIN" =~ ^[0-9]+$ || ! "$MAX" =~ ^[0-9]+$ || "$MIN" -gt "$MAX" ]] && {
  echo "Error: MIN/MAX must be integers and MIN <= MAX." >&2
  usage
}

# require bioawk
command -v bioawk >/dev/null 2>&1 || {
  echo "Error: bioawk not found in PATH." >&2
  exit 1
}

# --- process each fasta ---
for f in "$@"; do
  if [[ ! -f "$f" ]]; then
    echo "Warning: '$f' does not exist. Skipping." >&2
    continue
  fi

  base="${f%.*}"                   # strip last extension
  out="${base}_subset.fasta"

  bioawk -c fastx -v min="$MIN" -v max="$MAX" '
    length($seq) >= min && length($seq) <= max {
      print ">"$name
      print $seq
    }
  ' "$f" > "$out"

  echo "Wrote: $out"
done
