#!/usr/bin/env python3
"""Tensor EOM H3 = Dcom: ethS 223 vs nested 3b gold, then ⟨·⟩_ρ / 2.

H_EOM = H2 + H3 with
  H2 = ⟨χ | Htc([H,χ]_2)⟩_N
  H3 = Dcom = ⟨[χ, [χ, H]_3]⟩_ρ / 2

Nested gold (tensor ω, scalar H, leftover scalar 0+):
  W = −comm223st(H, ω) = [ω, H]_3
  Z = comm231tts(ω,W) + comm232tts(ω,W) + comm132tts(ω,W)
  gold = ⟨Z⟩_ρ / 2

ethS (what DcomMultiref uses):
  223_231_st + 223_232 + 223_132_tts

emax=1 He4+p-shell so tensor 3b is cheap. RDM = he8.ref (p-shell).

From run2/:
  PYTHONPATH=../build python3 test_dcom_tts.py
"""
import os
import sys

import numpy as np
from pyIMSRG import *

emax = int(os.environ.get("EMAX", "1"))
rank_j = int(os.environ.get("JRANK", "2"))
parity = int(os.environ.get("PARITY", "0"))
rank_Tz = int(os.environ.get("ITZ", "0"))
tol = float(os.environ.get("DCOM_TOL", "1e-6"))
ref_file = os.environ.get("REF_FILE", "he8.ref")
seed = int(os.environ.get("SEED", "17"))

eth = Commutator.FactorizedDoubleCommutator_eths
cm = Commutator


def make_Z():
    Z = Operator(ms, 0, 0, 0, 2)
    Z.SetHermitian()
    Z *= 0.0
    return Z


def make_W3():
    W = Operator(ms, rank_j, rank_Tz, parity, 3)
    W.SetHermitian()
    W.ThreeBody.SetMode("pn")
    W *= 0.0
    return W


def leftover_eths(omega, H):
    Z = make_Z()
    eth.SetUse_1b_Intermediates(True)
    eth.SetUse_2b_Intermediates(True)
    eth.comm223_231_st(omega, H, Z)
    eth.comm223_232(omega, H, Z)
    eth.comm223_132_tts(omega, H, Z)
    return Z


def leftover_nested(omega, H):
    W = make_W3()
    cm.comm223st(H, omega, W)
    W *= -1.0
    Z = make_Z()
    cm.comm231tts(omega, W, Z)
    cm.comm232tts(omega, W, Z)
    cm.comm132tts(omega, W, Z)
    return Z, W


def rel_err(a, b):
    return abs(a - b) / max(1.0, abs(a), abs(b))


ms = ModelSpace(emax, "He4", "p-shell")
ms.SetHbarOmega(16)
if not os.path.isfile(ref_file):
    raise FileNotFoundError(ref_file)

Hs = Operator(ms, 0, 0, 0, 2)
unt = UnitTest(ms)
unt.SetRandomSeed(seed)
Hs = unt.RandomOp(ms, 0, 0, 0, 2, 1)

eom = EOM(Hs, ref_file, rank_j, parity, rank_Tz)
eom.SetUseRdm3(False)
eom.SetIncludeConfigs(True, True, True, True, True)
eom.ConstructConfigs_tensor()
eom.ConstructNormMatrix_tensor()
eom.SetArnoldiUseProjection(False)

print(f"MS=He4+p-shell emax={emax}  JπTz={rank_j} {parity} {rank_Tz}  "
      f"EOM dim={eom.eom_dims}  seed={seed}")

unt.SetRandomSeed(seed)
h_rand = unt.RandomOp(ms, rank_j, rank_Tz, parity, 2, 1)
chi = eom.GetVSEOM_ladder_multiref(h_rand, 1)
eta = eom.GetVSEOM_ladder_multiref(chi, -1)

# Rank-3 copy so comm132tts can touch X3 (then unused).
eta3 = Operator(ms, rank_j, rank_Tz, parity, 3)
eta3.SetAntiHermitian()
eta3.ThreeBody.SetMode("pn")
eta3 *= 0.0
eom.UnflattenOperator(eta3, eom.FlattenOperator(eta))
# The tensor commutators (comm223st, comm*tts) take Ω as a reduced RME for
# every λ. At J=0 the IMSRG scalar convention is unreduced; ethS converts
# on entry, the nested gold needs it done here.
if rank_j == 0 and not eta3.IsReduced():
    eta3.MakeReduced()

print(f"  ||η1||={eta.OneBodyNorm():.4g}  ||η2||={eta.TwoBodyNorm():.4g}")

print("\n=== 1. leftover Z: nested 3b vs ethS ===")
Z_nested, W = leftover_nested(eta3, Hs)
Z_eths = leftover_eths(eta, Hs)
n_nested = Z_nested.Norm()
n_eths = Z_eths.Norm()
dZ = (Z_nested - Z_eths).Norm()
print(f"  ||Z_nested||={n_nested:.6e}  ||W3||={W.ThreeBodyNorm():.4g}")
print(f"  ||Z_eths||  ={n_eths:.6e}  ||ΔZ||={dZ:.3e}")
ok1 = dZ <= tol * max(1.0, n_nested)
print(f"  leftover match: {'OK' if ok1 else 'FAIL'}")

print("\n=== 2. DcomMultiref vs ⟨Z⟩_ρ / 2 ===")
gold = eom.GetVSEOM_Overlap_multiref(Z_nested) / 2.0
code = eom.DcomMultiref(Hs, chi)
eths_ov = eom.GetVSEOM_Overlap_multiref(Z_eths) / 2.0
print(f"  nested ⟨Z⟩/2     = {gold:.8e}")
print(f"  ethS  ⟨Z⟩/2      = {eths_ov:.8e}")
print(f"  DcomMultiref     = {code:.8e}")
print(f"  Dcom/nested      = {code / gold if abs(gold) > 1e-14 else np.nan:.6f}")
ok2 = rel_err(code, gold) <= tol
ok2b = rel_err(code, eths_ov) <= tol
print(f"  Dcom vs nested: {'OK' if ok2 else 'FAIL'}")
print(f"  Dcom vs ethS ov: {'OK' if ok2b else 'FAIL'}")

print("\n=== 3. polarization  H3(a,b) = ½(Dcom(a+b)−Dcom(a)−Dcom(b))  ===")
unt.SetRandomSeed(seed + 1)
chi_b = eom.GetVSEOM_ladder_multiref(
    unt.RandomOp(ms, rank_j, rank_Tz, parity, 2, 1), 1)
d_a = eom.DcomMultiref(Hs, chi)
d_b = eom.DcomMultiref(Hs, chi_b)
d_ab = eom.DcomMultiref(Hs, chi + chi_b)
cross = 0.5 * (d_ab - d_a - d_b)
d_ba = eom.DcomMultiref(Hs, chi_b + chi)
cross2 = 0.5 * (d_ba - d_a - d_b)
print(f"  Dcom(a)={d_a:.6e}  Dcom(b)={d_b:.6e}")
print(f"  H3_ab={cross:.6e}  H3_ba={cross2:.6e}  |ab−ba|={abs(cross-cross2):.3e}")
ok3 = rel_err(cross, cross2) <= tol
print(f"  H3_ab = H3_ba: {'OK' if ok3 else 'FAIL'}")

print("\n=== 4. H2 + H3 Rayleigh pieces (one vector) ===")
eom.SetArnoldiUseProjection(False)
h2 = eom.ComputeNorm(chi, eom.HtcMultiref(Hs, chi))
n0 = eom.ComputeNorm(chi, chi)
ray_h2 = h2 / n0 if abs(n0) > 1e-14 else np.nan
ray_h3 = d_a / n0 if abs(n0) > 1e-14 else np.nan
print(f"  N={n0:.6e}  H2={h2:.6e}  H3={d_a:.6e}")
print(f"  H2/N={ray_h2:.6f}  H3/N={ray_h3:.6f}  (H2+H3)/N={ray_h2+ray_h3:.6f}")

all_ok = ok1 and ok2 and ok2b and ok3
print("\nPASS" if all_ok else "\nFAIL")
sys.exit(0 if all_ok else 1)
