#!/usr/bin/env python3
"""He6: neutral MR-EOM roots vs full-space diagonalization of the same Hamiltonian.

With two valence nucleons rho3 = rho4 = 0, so the MR-EOM kernels are exact for
the two-body H. The valence (P) space is exactly decoupled from the rest (Q),
and every EOM vector chi|Psi> lies in Q. Rayleigh-Ritz in a subspace of Q then
gives, for each J, the bound

    k-th MR-EOM root  >=  k-th non-valence FCI eigenvalue.

A violation means a bug; the size of the gap measures how incomplete the
two-body operator pool is.

Inputs
  he6_fci/log_He6_fci_hs_j{0,4}p.txt      KSHELL full-space logs (he6_fci.slurm)
  he6_fci/fci_hs_He6.snt                  for the zero-body term
  run2/he6_rho3/he6_kshell_e3.json        p-shell states (to remove from FCI)
  eom_result/He6_*_rho3_full.json         MR-EOM (mr_rho3_full.py, NUCLEUS=He6)

  PYTHONPATH=../build python3 compare_he6_fci.py
"""
from __future__ import annotations

import json
import os

from he8_io import PAPER, RUN2, eom_result_path
from plot_he8_mr_eom import load_kshell_jlog

FCI_DIR = os.path.join(PAPER, "he6_fci")
TOL_MATCH = 2e-3   # MeV, identify a FCI level with a p-shell level
TOL_BOUND = 1e-6   # MeV, allowed violation of the Ritz bound

# MR-EOM jobs whose roots have the given J (reference J x operator rank)
JOBS = {0: ("mr0gs", "mr0x", "t20", "t20x"), 2: ("t02", "t02x", "mr2gs", "mr2x")}


def zero_body():
    with open(os.path.join(FCI_DIR, "fci_hs_He6.snt")) as f:
        for line in f:
            if "Zero body term" in line:
                return float(line.split(":")[1])
    raise RuntimeError("no zero-body term in fci_hs_He6.snt")


def main():
    zb = zero_body()
    ks = json.load(open(os.path.join(RUN2, "he6_rho3", "he6_kshell_e3.json")))
    valence = {J: sorted(s["E_sm"] + ks["E_core"] for s in ks["he8"] if s["twoJ"] == 2 * J)
               for J in (0, 2)}
    runs = json.load(open(eom_result_path("He6", 16, 3, "rho3_full")))["runs"]

    n_viol = n_cmp = 0
    for J in (0, 2):
        log = os.path.join(FCI_DIR, f"log_He6_fci_hs_j{2 * J}p.txt")
        if not os.path.isfile(log):
            print(f"J={J}: missing {log}")
            continue
        fci = sorted(s["E"] + zb for s in load_kshell_jlog(log, j_want=J))
        q = list(fci)
        print(f"\nHe6 J={J}+   FCI levels: {len(fci)}   zero-body {zb:.6f}")
        for ev in valence[J]:
            hit = min(q, key=lambda e: abs(e - ev))
            ok = abs(hit - ev) < TOL_MATCH
            print(f"  p-shell {ev:10.5f}  FCI {hit:10.5f}  diff {hit - ev:+.1e}  "
                  f"{'removed' if ok else 'NOT FOUND (decoupling?)'}")
            if ok:
                q.remove(hit)
        print("  non-valence FCI:", " ".join(f"{e:8.3f}" for e in q[:10]))
        for name in JOBS[J]:
            rec = runs.get(f"{name}_rdm3")
            if not rec or "error" in rec:
                continue
            roots = sorted(rec["E_abs"])
            line = []
            for k, e in enumerate(roots):
                if k >= len(q):
                    break
                gap = e - q[k]
                bad = gap < -TOL_BOUND
                n_viol += int(bad)
                n_cmp += 1
                line.append(f"{e:8.3f}({gap:+.2f}{'!' if bad else ''})")
            print(f"  {name:6s} ref {rec['E_ref']:9.3f}:  " + " ".join(line))
    if n_cmp == 0:
        print("\nno MR-EOM root could be compared (FCI logs missing?)")
        return
    print(f"\nbound checks: {n_cmp}   violations: {n_viol}  ->  {'FAIL' if n_viol else 'PASS'}")
    print("entries: MR root (MR root - k-th non-valence FCI level); '!' marks a violation")


if __name__ == "__main__":
    main()
