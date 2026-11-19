#!/usr/bin/env python3
"""Benchmark ConstructOverlapMatrix vs ConstructNormMatrix.

Hermitian RDM: overlap kernel must match N[ρ]; Nkernel itself unchanged.
Same-ref OverlapAB must match ComputeNorm.

TDM: N[ρ] still matches ConstructNormMatrix on the TDM. The reverse
contraction is N[ρ†] of the same operator, and (Qa†Qb)·ρ = (Qb†Qa)·ρ†.
"""
from __future__ import annotations

import os
import sys

import numpy as np

from he8_io import BUILD
from he8_mr_eom import JOBS, Hs_file, build_hs, resolve_ref

if BUILD not in sys.path:
    sys.path.insert(0, BUILD)

from pyIMSRG import EOM  # noqa: E402

PAPER = os.path.dirname(os.path.abspath(__file__))
TRANSIT = os.path.join(PAPER, "kshell_he8", "transit.log")
SNT = os.path.join(PAPER, "kshell_he8", "he8.snt")
GS_RITZ = os.path.join(os.path.dirname(PAPER), "eom_result", "ritz_live",
                       "he8_mr0gs_ritz")


def maxabs(a, b=None):
    if b is None:
        return float(np.max(np.abs(a)))
    return float(np.max(np.abs(a - b)))


def main():
    ms, Hs = build_hs()
    gs_ref = resolve_ref(JOBS["mr0gs"])["ref"]
    eom = EOM(Hs, gs_ref, 0, 0, 0)
    eom.SetUseRdm3(False)
    eom.ConstructConfigs()
    eom.ConstructNormMatrix(True)
    N_herm = np.array(eom.Nkernel, copy=True)
    print(f"N[ρ] fermi-on  shape={N_herm.shape}  |N|={maxabs(N_herm):.6e}",
          flush=True)

    eom.ConstructNormMatrix(False)
    N_fermi_off = np.array(eom.Nkernel, copy=True)
    eom.ConstructNormMatrix(True)
    dN = maxabs(np.array(eom.Nkernel), N_herm)
    print(f"reload fermi-on after fermi-off  |ΔN|={dN:.3e}", flush=True)
    if dN > 1e-12:
        raise RuntimeError("ConstructNormMatrix(True) not reproducible")

    eom.ConstructOverlapMatrix(True)
    N_ov_rdm = np.array(eom.Nkernel_overlap, copy=True)
    d_keep = maxabs(np.array(eom.Nkernel), N_herm)
    d_match = maxabs(N_ov_rdm, N_herm)
    print(f"ConstructOverlapMatrix(True) on ρ:  |Nkernel−N_herm|={d_keep:.3e}  "
          f"|N_ov−N[ρ]|={d_match:.3e}", flush=True)
    if d_keep > 1e-12:
        raise RuntimeError("ConstructOverlapMatrix overwrote Nkernel")
    if d_match > 1e-10:
        raise RuntimeError("same-ref overlap kernel must match N[ρ] (ComputeNorm)")

    from pyIMSRG import ReadWrite, Operator
    rw = ReadWrite()
    qs = []
    for k in range(5):
        q = Operator(ms, 0, 0, 0, 2)
        rw.ReadOperator(q, f"{GS_RITZ}_{k}.op")
        qs.append(q)
    dS = 0.0
    for i in range(5):
        for j in range(5):
            n = eom.ComputeNorm(qs[i], qs[j])
            o = eom.OverlapAB(qs[i], qs[j])
            dS = max(dS, abs(n - o))
    print(f"max |OverlapAB − ComputeNorm| on 5 Ritz = {dS:.3e}", flush=True)
    if dS > 1e-10:
        raise RuntimeError("same-ref OverlapAB must match ComputeNorm")

    eom.ConstructOverlapMatrix(False)
    d_off = maxabs(np.array(eom.Nkernel_overlap), N_fermi_off)
    print(f"ConstructOverlapMatrix(False) vs N fermi-off  |Δ|={d_off:.3e}",
          flush=True)
    if d_off > 1e-10:
        raise RuntimeError("fermi-off overlap kernel must match N fermi-off")
    eom.ConstructNormMatrix(True)

    probe = EOM(Hs, 0, 0, 0)
    probe.SetUseRdm3(False)
    tdm = probe.ReadKshellRdm(TRANSIT, SNT, 1, 2)
    v3344 = tdm.TwoBody.GetTBME_J_norm(0, 0, 2, 2, 3, 3)
    v4433 = tdm.TwoBody.GetTBME_J_norm(0, 0, 3, 3, 2, 2)
    print(f"TDM 3344={v3344:.8f}  4433={v4433:.8f}", flush=True)

    gold = EOM(Hs, 0, 0, 0)
    gold.SetUseRdm3(False)
    gold.LoadKshellRdm(TRANSIT, SNT, 1, 2)
    gold.ConstructConfigs()
    gold.ConstructNormMatrix(True)
    N_old = np.array(gold.Nkernel, copy=True)

    eom.ConstructOverlapKernel(tdm)
    N_new = np.array(eom.Nkernel_overlap, copy=True)
    N_dag = np.array(eom.Nkernel_overlap_dag, copy=True)
    d_keep2 = maxabs(np.array(eom.Nkernel), N_herm)
    d_tdm = maxabs(N_new, N_old)
    print(f"TDM N[ρ] vs ConstructNormMatrix on TDM  |Δ|={d_tdm:.6e}  "
          f"|Nkernel still herm|={d_keep2:.3e}", flush=True)
    if d_keep2 > 1e-12:
        raise RuntimeError("ConstructOverlapKernel overwrote Nkernel")
    if d_tdm > 1e-10:
        raise RuntimeError("TDM N[ρ] must match ConstructNormMatrix")

    d_fwd_rev = 0.0
    d_avg = 0.0
    for i in range(min(3, len(qs))):
        for j in range(min(3, len(qs))):
            s_ab = eom.OverlapAB(qs[i], qs[j])
            s_ba = eom.OverlapBA(qs[j], qs[i])
            s_h = eom.OverlapHermitian(qs[i], qs[j])
            d_fwd_rev = max(d_fwd_rev, abs(s_ab - s_ba))
            d_avg = max(d_avg, abs(s_h - 0.5 * (s_ab + s_ba)))
    print(f"max |(Qa†Qb)·ρ − (Qb†Qa)·ρ†| = {d_fwd_rev:.3e}", flush=True)
    print(f"max |OverlapHermitian − ½(fwd+rev)| = {d_avg:.3e}", flush=True)
    if d_fwd_rev > 1e-10:
        raise RuntimeError("cross-ref product must satisfy (Qa†Qb)·ρ = (Qb†Qa)·ρ†")
    if d_avg > 1e-12:
        raise RuntimeError("OverlapHermitian is not the average of ρ and ρ†")

    asym_new = maxabs(N_new - N_new.T)
    asym_old = maxabs(N_old - N_old.T)
    print(f"|N[ρ]−N[ρ]ᵀ|={asym_new:.6e}  |N_partner−N_partnerᵀ|={asym_old:.6e}  "
          f"|N[ρ†]|={maxabs(N_dag):.6e}", flush=True)
    print("PASS", flush=True)


if __name__ == "__main__":
    os.chdir(PAPER)
    main()
