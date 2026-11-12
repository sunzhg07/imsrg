#!/usr/bin/env python3
"""Tensor EOM checks on Mg24 (emax 4).

1. Block projector (ConstructProjectMatrix_blocks) vs the old single dense
   SVD of all of N: timing, P v agreement on the E1 seed and random vectors,
   and idempotency P P v = P v.
2. Metric sign at lambda=2: single particle-hole chi = (a^dag i)^2, core 0s1/2
   -> valence 0d5/2 / 0d3/2. Expect N = 1 - n_a (physical metric).

From run_paper/:  python3 check_projector.py      (PART=1 or PART=2 for one)
"""
from __future__ import annotations

import math
import os
import sys
import time

from he8_io import BUILD, RUN2, hs_path

if BUILD not in sys.path:
    sys.path.insert(0, BUILD)

from pyIMSRG import EOM, ModelSpace, Operator, ReadWrite, UnitTest  # noqa: E402

emax, hw = 4, 16
Hs_file = hs_path("Mg24", hw, emax, val="sd-shell")
E1_file = os.path.join(os.path.dirname(Hs_file), "E1c" + os.path.basename(Hs_file)[len("Hs"):])
ref_file = os.path.join(RUN2, "mg24", f"mg24_gs_e{emax}.ref")


def load():
    ms = ModelSpace(emax, "O16", "sd-shell")
    ms.SetHbarOmega(hw)
    rw = ReadWrite()
    Hs = Operator(ms, 0, 0, 0, 3)
    rw.ReadOperator(Hs, Hs_file)
    return ms, rw, Hs


def build(Hs, J2, parity):
    eom = EOM(Hs, ref_file, J2, parity, 0)
    eom.SetUseRdm3(False)
    eom.PrepareHamiltonianForArnoldi()
    eom.ConstructConfigs_tensor()
    t = time.time()
    eom.ConstructNormMatrix_tensor()
    print(f"ConstructNormMatrix_tensor: {time.time() - t:.1f} s", flush=True)
    return eom


def part1(ms, rw, Hs):
    eom = build(Hs, 1, 1)
    E1 = Operator(ms, 1, 0, 1, 2)
    rw.ReadOperator(E1, E1_file)
    ut = UnitTest(ms)
    ut.SetRandomSeed(5)
    vecs = [("E1 seed", eom.GetVSEOM_ladder_multiref(E1, 1))]
    for k in range(3):
        vecs.append((f"random{k}", eom.GetVSEOM_ladder_multiref(ut.RandomOp(ms, 1, 0, 1, 2, 1), 1)))

    def project(v):
        w = v * 1.0
        eom.ProjectOprator(w)
        return w

    t = time.time()
    eom.ConstructProjectMatrix()
    t_blocks = time.time() - t
    Pb = [project(v) for _, v in vecs]
    PPb = [project(w) for w in Pb]
    t = time.time()
    for _, v in vecs:
        project(v)
    t_apply_b = (time.time() - t) / len(vecs)

    eom.SetProjectorFullSVD(True)
    t = time.time()
    eom.ConstructProjectMatrix()
    t_full = time.time() - t
    t = time.time()
    Pf = [project(v) for _, v in vecs]
    t_apply_f = (time.time() - t) / len(vecs)

    print(f"\nprojector build: blocks {t_blocks:.1f} s   full SVD {t_full:.1f} s   "
          f"(x{t_full / max(t_blocks, 1e-9):.0f})", flush=True)
    print(f"ProjectOprator per call: blocks {t_apply_b:.3f} s   full {t_apply_f:.3f} s", flush=True)
    worst = 0.0
    for (lab, v), pb, ppb, pf in zip(vecs, Pb, PPb, Pf):
        d = (pb - pf).Norm() / max(pf.Norm(), 1e-300)
        idem = (ppb - pb).Norm() / max(pb.Norm(), 1e-300)
        worst = max(worst, d)
        print(f"  {lab:<8} |v|={v.Norm():.4e} |Pv|={pf.Norm():.4e}  "
              f"|P_blk v - P_full v|/|Pv|={d:.2e}  |PPv-Pv|/|Pv|={idem:.2e}  "
              f"N(Pv): blk {eom.ComputeNorm(pb, pb):+.6e}  full {eom.ComputeNorm(pf, pf):+.6e}",
              flush=True)
    print(f"PART1 {'PASS' if worst < 1e-8 else 'FAIL'}  (worst rel diff {worst:.2e})", flush=True)


def part2(ms, Hs):
    lam = 2
    eom = build(Hs, lam, 0)
    occ = {(0, 2, 5): 2.952 / 6, (0, 2, 3): 0.635 / 4}
    for (ni, li, ji), (na, la, ja) in (((0, 0, 1), (0, 2, 5)), ((0, 0, 1), (0, 2, 3))):
        i = ms.GetOrbitIndex(ni, li, ji, -1)
        a = ms.GetOrbitIndex(na, la, ja, -1)
        x = Operator(ms, lam, 0, 0, 2)
        x.SetOneBodyME(a, i, 1.0)
        chi = eom.GetVSEOM_ladder_multiref(x, 1)
        N = eom.ComputeNorm(chi, chi)
        n = occ[(na, la, ja)]
        exp = 1 - n   # physical metric (EOM.cc multiplies by (-1)^lam lam_hat)
        print(f"PART2 lam=2 core ({ni},{li},{ji}/2) -> val ({na},{la},{ja}/2): N={N:+.6e}  "
              f"expected 1-n={exp:+.6e}  ratio={N / exp:+.4f}", flush=True)


def main():
    ms, rw, Hs = load()
    part = os.environ.get("PART", "12")
    if "1" in part:
        part1(ms, rw, Hs)
    if "2" in part:
        part2(ms, Hs)


if __name__ == "__main__":
    main()
