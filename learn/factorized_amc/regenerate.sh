#!/bin/bash
# Regenerate every AMC output in learn/factorized_amc from its input.
#
#   learn/factorized_amc/regenerate.sh            # rewrite outputs in place
#   learn/factorized_amc/regenerate.sh --check    # regenerate to a temp dir and diff
#
# Flags follow how each stored output was produced:
#   <stem>.tex        --wet-convention wigner
#   <stem>_ninej.tex  --wet-convention wigner --collect-ninejs
#   leftover132/, nested_tts/  always --collect-ninejs (stored as <stem>.tex)
# Inputs without a stored output are generated in plain mode.
set -u
DIR="$(cd "$(dirname "$0")" && pwd)"
PY="${AMC_PYTHON:-python3}"
CHECK=0
[ "${1:-}" = "--check" ] && CHECK=1
TMP=$(mktemp -d)
fail=0

for f in "$DIR"/*/input/*.txt; do
  d=$(basename "$(dirname "$(dirname "$f")")")
  s=$(basename "$f" .txt)
  o="$DIR/$d/output"
  modes=()
  case "$d" in
    leftover132|nested_tts) modes=("$s.tex|--collect-ninejs") ;;
    *)
      [ -f "$o/$s.tex" ] && modes+=("$s.tex|")
      [ -f "$o/${s}_ninej.tex" ] && modes+=("${s}_ninej.tex|--collect-ninejs")
      [ ${#modes[@]} -eq 0 ] && modes=("$s.tex|")
      ;;
  esac
  for m in "${modes[@]}"; do
    name=${m%%|*}; flags=${m#*|}
    if [ $CHECK -eq 1 ]; then dest="$TMP/${d}__$name"; else dest="$o/$name"; fi
    if ! "$PY" -m amc "$f" -o "$dest" --wet-convention wigner $flags > "$TMP/log" 2>&1; then
      echo "AMC-FAIL $d/$name"; fail=1; continue
    fi
    if [ $CHECK -eq 1 ]; then
      if cmp -s "$dest" "$o/$name"; then echo "same     $d/$name"; else echo "DIFF     $d/$name"; fail=1; fi
    else
      echo "wrote    $d/$name"
    fi
  done
done
rm -rf "$TMP"
exit $fail
