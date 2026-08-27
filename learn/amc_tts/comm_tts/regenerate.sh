#!/bin/sh
# Regenerate J-coupled equations for scalar TTS commutators.
# Run from repo root:  learn/amc_tts/comm_tts/regenerate.sh
set -e
DIR="$(cd "$(dirname "$0")" && pwd)"
OUT="$DIR/output"
mkdir -p "$OUT"
for f in "$DIR"/input/*_unred.txt; do
  base=$(basename "$f" .txt)
  echo "AMC $base ..."
  python3 -m amc "$f" -o "$OUT/${base}.tex" --wet-convention wigner
done
echo "Done. See AMC_CHECK.md for packaging notes."
