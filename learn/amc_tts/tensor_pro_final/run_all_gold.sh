#!/usr/bin/env bash
# Smoke-run gold benches for tensor_pro_final (emax=1, λ=2 by default).
# Continues on failure; exits 1 if any failed.
# Usage: ./run_all_gold.sh [emax] [lambda]
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$ROOT"
export PYTHONPATH="${PYTHONPATH:-build}"
EMAX="${1:-1}"
LAM="${2:-2}"
# Per-test wall clock cap. A bench that blows past this is a bug (usually an
# argv-order slip feeding a huge emax into a deeply nested AMC sum), not a slow
# machine -- fail it loudly instead of hanging the whole pack.
TIMEOUT="${GOLD_TIMEOUT:-900}"
TIMEOUT_BIN="$(command -v timeout || command -v gtimeout || true)"

pass=0
fail=0
timedout=0
TIMES_FILE="$(mktemp -t gold_times)"
trap 'rm -f "$TIMES_FILE"' EXIT

run_one() {
  local label="$1"
  shift
  echo ""
  echo "======== $label ========"
  local t0 t1 elapsed
  t0=$(python3 -c 'import time; print(time.time())')
  if [ -n "$TIMEOUT_BIN" ]; then
    "$TIMEOUT_BIN" "$TIMEOUT" "$@"
  else
    "$@"
  fi
  local rc=$?
  t1=$(python3 -c 'import time; print(time.time())')
  elapsed=$(python3 -c "print(f'{$t1 - $t0:.2f}')")
  local status
  if [ $rc -eq 0 ]; then
    status=OK
    pass=$((pass + 1))
  elif [ $rc -eq 124 ]; then
    status=TIMEOUT
    timedout=$((timedout + 1))
    fail=$((fail + 1))
  else
    status="FAIL(exit $rc)"
    fail=$((fail + 1))
  fi
  echo "$status $label  [${elapsed}s]"
  printf '%s\t%s\t%s\n' "$elapsed" "$status" "$label" >> "$TIMES_FILE"
}

run_one "fI/fII/fIII mscheme" python3 -B run/test_tts_f_mscheme.py "$EMAX" "$LAM"
run_one "fII chi_beta" python3 -B run/test_chi_beta_mscheme.py "$EMAX" "$LAM"
run_one "fIIIa chi_gamma PathB" python3 -B run/test_chi_gamma_pathB_amc.py "$EMAX" "$LAM"
run_one "fIIIb" python3 -B run/test_tts_fIIIb.py "$EMAX" "$LAM"

run_one "GI chi_epsilon" python3 -B run/test_chi_epsilon_mscheme.py "$EMAX" "$LAM"
run_one "GI mscheme" python3 -B run/test_tts_GI_mscheme.py "$EMAX" "$LAM"

run_one "GII chi_zeta" python3 -B run/test_chi_zeta_mscheme.py "$EMAX" "$LAM"
run_one "GII pathB mscheme" python3 -B run/test_tts_GII_pathB_mscheme.py "$EMAX" "$LAM"

run_one "GIIIa chi_eta" python3 -B run/test_chi_eta_mscheme.py "$EMAX" "$LAM"
run_one "GIIIa ladder" python3 -B run/test_GIIIa_ladder_mscheme.py "$EMAX" "$LAM"
run_one "GIIIa Path B" python3 -B run/test_tts_GIIIa.py "$EMAX" "$LAM"

run_one "GIIIb Path B" python3 -B run/test_tts_GIIIb.py "$EMAX" "$LAM"

run_one "GIIIc chi_theta" python3 -B run/test_chi_theta_mscheme.py "$EMAX" "$LAM"
run_one "GIIIc mscheme" python3 -B run/test_tts_GIIIc_mscheme.py "$EMAX" "$LAM"

run_one "GIVa chi PathB" python3 -B run/test_chi_kappa_pathB_vs_direct.py "$EMAX" "$LAM"
run_one "GIVa full PathB" python3 -B run/test_G4a_pathB_mscheme.py "$EMAX" "$LAM" 15
run_one "GIVa ethS vs PathB" python3 -B run/test_tts_GIVa_eths_vs_pathB.py "$EMAX" "$LAM"

run_one "GIVc chi_lambda" python3 -B run/test_chi_lambda_mscheme.py "$EMAX" "$LAM"
run_one "GIVc mscheme" python3 -B run/test_tts_GIVc_mscheme.py "$EMAX" "$LAM"
run_one "GIVc PathB" python3 -B run/test_tts_GIVc_pathB.py "$EMAX" "$LAM"

echo ""
echo "======== TIMING (slowest first) ========"
sort -rn "$TIMES_FILE" | awk -F'\t' '{printf "%9.2fs  %-14s %s\n", $1, $2, $3}'
awk -F'\t' '{t += $1} END {printf "%9.2fs  TOTAL\n", t}' "$TIMES_FILE"

echo ""
echo "======== SUMMARY ========"
echo "pass=$pass fail=$fail timeout=$timedout"
exit $(( fail > 0 ? 1 : 0 ))
