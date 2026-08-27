#!/usr/bin/env python3
"""Does the tensor Pandya kernel reduce to the scalar one at lambda=0?

This is the precondition for having a single universal-lambda Pandya routine
instead of separate lambda==0 / lambda!=0 branches.

Scalar kernel (FactorizedDoubleCommutator.cc:1290-1293):
    Xbar -= (2J+1) * SixJ(ja,jb,Jcc,jc,jd,J) * X^{J}_{adcb}

Tensor kernel, two conventions live in FactorizedDoubleCommutator_eths.cc:
  GIIIa  pandya_eta      (:4109-4120)  hats * scale, scale=sqrt((2J1+1)/(2Jbra+1))
  GIVc   pandya_tensor   (:6831-6842)  hats only, no scale
Both are labelled "IMSRG tensor Pandya"; they differ by `scale`. At most one
can reduce to the scalar kernel.

Usage:
  PYTHONPATH=build python3 -B run/test_pandya_lambda0_reduction.py
"""

from __future__ import annotations

import math
import sys

from pyIMSRG import *


def hat2(J: float) -> float:
    return 2.0 * J + 1.0


def phase(n: int) -> float:
    return 1.0 if n % 2 == 0 else -1.0


def tri(a: float, b: float, c: float) -> bool:
    return abs(a - b) <= c <= a + b


def scalar_kernel(ja, jb, jc, jd, Jcc, J) -> float:
    """(2J+1) * SixJ(ja,jb,Jcc,jc,jd,J), with the leading minus factored out."""
    return hat2(J) * SixJ(ja, jb, Jcc, jc, jd, J)


def tensor_kernel(ja, jb, jc, jd, Jbra, Jket, J1, J2, lam, with_scale) -> float:
    """hats * [scale] * phase(jb+jd+Jket+J2) * NineJ(...), minus factored out."""
    ninej = NineJ(ja, jd, J1, jb, jc, J2, Jbra, Jket, lam)
    if abs(ninej) < 1e-14:
        return 0.0
    hats = math.sqrt(hat2(J1) * hat2(J2) * hat2(Jbra) * hat2(Jket))
    scale = math.sqrt(hat2(J1) / hat2(Jbra)) if with_scale else 1.0
    ph = phase(int(round(jb + jd)) + int(Jket) + int(J2))
    return hats * scale * ph * ninej


# Half-integer single-particle j values, as in a real model space.
jvals = [0.5, 1.5, 2.5, 3.5]

results = {True: {"n": 0, "max": 0.0, "worst": None, "ratios": {}},
           False: {"n": 0, "max": 0.0, "worst": None, "ratios": {}}}

for ja in jvals:
    for jb in jvals:
        for jc in jvals:
            for jd in jvals:
                # Jcc couples (a,d-bar) and (c,b-bar) in the CC channel.
                Jmin = int(max(abs(ja - jd), abs(jc - jb)))
                Jmax = int(min(ja + jd, jc + jb))
                for Jcc in range(Jmin, Jmax + 1):
                    # scalar sum index J couples (a,b) and (c,d)
                    for J in range(0, 8):
                        if not (tri(ja, jd, J) and tri(jc, jb, J)):
                            continue
                        sk = scalar_kernel(ja, jb, jc, jd, Jcc, J)
                        for with_scale in (True, False):
                            # lambda=0 forces J2=J1=J and Jket=Jbra=Jcc
                            tk = tensor_kernel(
                                ja, jb, jc, jd, Jcc, Jcc, J, J, 0, with_scale
                            )
                            if abs(sk) < 1e-13 and abs(tk) < 1e-13:
                                continue
                            r = results[with_scale]
                            r["n"] += 1
                            err = abs(sk - tk)
                            if err > r["max"]:
                                r["max"] = err
                                r["worst"] = (ja, jb, jc, jd, Jcc, J, sk, tk)
                            if abs(tk) > 1e-13:
                                key = round(sk / tk, 6)
                                r["ratios"][key] = r["ratios"].get(key, 0) + 1

tol = 1e-10
print("lambda=0 reduction: tensor Pandya kernel vs scalar Pandya kernel\n")
ok_any = False
for with_scale, label in ((True, "GIIIa convention (hats*scale)"),
                          (False, "GIVc convention (hats only) ")):
    r = results[with_scale]
    top = sorted(r["ratios"].items(), key=lambda kv: -kv[1])[:3]
    passed = r["max"] < tol
    ok_any = ok_any or passed
    print(f"{label}: n={r['n']} max|scalar-tensor|={r['max']:.3e} "
          f"=> {'REDUCES' if passed else 'DOES NOT REDUCE'}")
    print(f"    scalar/tensor ratios (top 3): {top}")
    if not passed and r["worst"]:
        ja, jb, jc, jd, Jcc, J, sk, tk = r["worst"]
        print(f"    worst: j=({ja},{jb},{jc},{jd}) Jcc={Jcc} J={J} "
              f"scalar={sk:.6e} tensor={tk:.6e}")
    print()

print("A universal-lambda Pandya is possible iff exactly one convention reduces.")
sys.exit(0 if ok_any else 1)
