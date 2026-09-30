#!/usr/bin/env python3
"""He8 0+ MR-EOM cross-reference overlap.

States: Q_a^i |0_gs⟩  (mr0gs)  and  Q_b^j |0_1⟩  (mr0x).
Live Arnoldi Ritz vectors Q = V y are N-orthonormal on the
same reference: ComputeNorm(Q_i, Q_j) = δ_ij.

Cross block is ⟨a|Q_a† Q_b|b⟩ from the same Wick diagrams as
ComputeNorm, including A1/A2, with ρ replaced by the kshell TDM
(transit.log states 1↔2).
The leftover commutator CrossRefOverlap is not this product.

The full S_aa / S_bb blocks use the Hermitian N[ρ].
Cross S_ab uses TDM ⟨gs|·|exc⟩, S_ba uses TDM ⟨exc|·|gs⟩.
With N[ρ]ᵀ = N[ρᵀ] these are adjoints: S_ba ≈ S_abᵀ.

"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

from he8_io import BUILD, RUN2, eom_result_path
from he8_mr_eom import JOBS, Hs_file, build_hs, nucleus, resolve_ref

if BUILD not in sys.path:
    sys.path.insert(0, BUILD)

from pyIMSRG import EOM, Operator, ReadWrite  # noqa: E402

TRANSIT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "kshell_he8", "transit.log")
SNT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "kshell_he8", "he8.snt")
GS_RITZ = os.path.join(RUN2, "he8_mr_ritz")
X_RITZ = os.path.join(RUN2, "he8", "he8_exc0_mr_ritz")


def n_states():
    return int(os.environ.get("N_STATES", "5"))


def load_ritz(ms, prefix, n):
    rw = ReadWrite()
    out = []
    for k in range(n):
        path = f"{prefix}_{k}.op"
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
        q = Operator(ms, 0, 0, 0, 2)
        rw.ReadOperator(q, path)
        out.append(q)
        print(f"  read {path}", flush=True)
    return out


def collect_ritz(run_result, n):
    ritz = list(run_result.arnoldi.ritz)
    if len(ritz) < n:
        raise RuntimeError(f"only {len(ritz)} Ritz vectors, need {n}")
    return ritz[:n]


def pack_ha(eom, qs):
    ha, aa = [], []
    for q in qs:
        ha.append(eom.GetVSEOM_ladder_multiref(q, 1))
        aa.append(eom.GetVSEOM_ladder_multiref(q, -1))
    return ha, aa


def block_same(eom, qs, label, use_compute_norm=True):
    n = len(qs)
    S = np.zeros((n, n))
    name = "ComputeNorm" if use_compute_norm else "NormMultiref"
    for i in range(n):
        for j in range(n):
            S[i, j] = (eom.ComputeNorm(qs[i], qs[j]) if use_compute_norm
                       else eom.NormMultiref(qs[i], qs[j]))
    print(f"\n=== same-ref {label}  {name} ===", flush=True)
    print(np.array2string(S, precision=6, suppress_small=True), flush=True)
    print(f"max |S-I| = {np.max(np.abs(S - np.eye(n))):.3e}", flush=True)
    return S


def block_cross(eom, qs_a, qs_b):
    na, nb = len(qs_a), len(qs_b)
    S = np.zeros((na, nb))
    for i in range(na):
        for j in range(nb):
            S[i, j] = eom.OverlapAB(qs_a[i], qs_b[j])
            print(f"  S[{i},{j}] = {S[i, j]:+.8f}", flush=True)
    return S


def fmt_mat(M, names):
    w = 10
    head = " " * 12 + "".join(f"{n:>{w}}" for n in names)
    lines = [head]
    for i, n in enumerate(names):
        row = "".join(f"{M[i, j]:{w}.4f}" for j in range(len(names)))
        lines.append(f"{n:<12}{row}")
    return "\n".join(lines)


def main():
    n = n_states()
    load = os.environ.get("LOAD_RITZ", "0") != "0"
    print(f"Hs={Hs_file}  N_STATES={n}  LOAD_RITZ={int(load)}", flush=True)
    print(f"TDM={TRANSIT}", flush=True)

    ms, Hs = build_hs()

    gs_ref = resolve_ref(JOBS["mr0gs"])["ref"]
    x_ref = resolve_ref(JOBS["mr0x"])["ref"]

    eom_gs = EOM(Hs, gs_ref, 0, 0, 0)
    eom_gs.SetUseRdm3(False)
    eom_x = EOM(Hs, x_ref, 0, 0, 0)
    eom_x.SetUseRdm3(False)

    if load:
        print("\nloading saved Ritz operators", flush=True)
        qs_gs = load_ritz(ms, os.environ.get("GS_RITZ", GS_RITZ), n)
        qs_x = load_ritz(ms, os.environ.get("X_RITZ", X_RITZ), n)
        eom_gs.ConstructConfigs()
        eom_gs.ConstructNormMatrix(True)
        eom_x.ConstructConfigs()
        eom_x.ConstructNormMatrix(True)
        e_gs = []
        e_x = []
    else:
        print("\nrunning MR-EOM mr0gs / mr0x", flush=True)
        max_iter = int(os.environ.get("MAX_ITER", "50"))
        use_h3 = os.environ.get("USE_H3", "1") != "0"
        use_prj = os.environ.get("USE_PROJECTION", "1") != "0"
        eom_gs.SetArnoldiUseH3(use_h3)
        eom_gs.SetArnoldiUseProjection(use_prj)
        res_gs = eom_gs.Run(max_iter, n)
        qs_gs = collect_ritz(res_gs, n)
        eom_x.SetArnoldiUseH3(use_h3)
        eom_x.SetArnoldiUseProjection(use_prj)
        res_x = eom_x.Run(max_iter, n)
        qs_x = collect_ritz(res_x, n)
        e_gs = [float(x) for x in res_gs.arnoldi.energies[:n]]
        e_x = [float(x) for x in res_x.arnoldi.energies[:n]]
        print(f"mr0gs ΔE = {e_gs}", flush=True)
        print(f"mr0x  ΔE = {e_x}", flush=True)
        rw = ReadWrite()
        outdir = os.path.join(os.path.dirname(eom_result_path(
            nucleus, 16, int(os.environ.get("EMAX", "3")), "mr0_overlap")),
            "ritz_live")
        os.makedirs(outdir, exist_ok=True)
        for tag, qs in (("gs", qs_gs), ("x", qs_x)):
            for k, q in enumerate(qs):
                path = os.path.join(outdir, f"he8_mr0{tag}_ritz_{k}.op")
                rw.WriteOperator(q, path)
                print(f"  wrote {path}", flush=True)

    # Same-ref Arnoldi metric: must be I on live Ritz
    S_aa = block_same(eom_gs, qs_gs, "0_gs", True)
    S_bb = block_same(eom_x, qs_x, "0_x", True)
    S_aa_leftover = block_same(eom_gs, qs_gs, "0_gs leftover", False)
    S_bb_leftover = block_same(eom_x, qs_x, "0_x leftover", False)

    # Sanity: four-commutator on the same ρ recovers NormMultiref
    print("\n=== sanity  CrossRef(ρ_a, zb=1) vs NormMultiref  (0_gs) ===",
          flush=True)
    Ha_g, Aa_g = pack_ha(eom_gs, qs_gs)
    max_d = 0.0
    for i in range(min(2, n)):
        for j in range(min(2, n)):
            c = eom_gs.CrossRefOverlapHA(Ha_g[i], Aa_g[i], Ha_g[j], Aa_g[j], True)
            m = S_aa[i, j]
            d = c - m
            max_d = max(max_d, abs(d))
            print(f"  ({i},{j})  Cross={c:+.8f}  N={m:+.8f}  Δ={d:+.3e}",
                  flush=True)
    print(f"max |Cross−N| on 2×2 = {max_d:.3e}", flush=True)

    # TDM 1→2 and 2→1. ComputeNorm on eom_gs / eom_x still uses ρ.
    # OverlapAB uses a separate N[tdm] kernel.
    # Keep two probes: ReadKshellRdm rebuilds rdm_modelspace, which the
    # TDM Operator points at.
    probe_ab = EOM(Hs, 0, 0, 0)
    probe_ba = EOM(Hs, 0, 0, 0)
    tdm_ab = probe_ab.ReadKshellRdm(TRANSIT, SNT, 1, 2)
    tdm_ba = probe_ba.ReadKshellRdm(TRANSIT, SNT, 2, 1)
    v3344 = tdm_ab.TwoBody.GetTBME_J_norm(0, 0, 2, 2, 3, 3)
    v4433 = tdm_ab.TwoBody.GetTBME_J_norm(0, 0, 3, 3, 2, 2)
    print("\nTDM(1,2) TBTD J=0  3344={:.8f}  4433={:.8f}  (must differ)"
          .format(v3344, v4433), flush=True)

    print("\n=== cross  ⟨0_gs| Q_a† Q_b |0_x⟩  TDM(1,2) diagrams ===",
          flush=True)
    eom_gs.ConstructOverlapKernel(tdm_ab)
    S_ab = block_cross(eom_gs, qs_gs, qs_x)

    print("\n=== cross  ⟨0_x| Q_b† Q_a |0_gs⟩  TDM(2,1) diagrams ===",
          flush=True)
    eom_x.ConstructOverlapKernel(tdm_ba)
    S_ba = block_cross(eom_x, qs_x, qs_gs)

    names = [f"a{i}" for i in range(n)] + [f"b{j}" for j in range(n)]
    S = np.zeros((2 * n, 2 * n))
    S[:n, :n] = S_aa
    S[n:, n:] = S_bb
    S[:n, n:] = S_ab
    S[n:, :n] = S_ba

    print("\n=== full S (aa, ab; ba, bb) ===", flush=True)
    print(fmt_mat(S, names), flush=True)
    asym = np.max(np.abs(S - S.T))
    print(f"\nmax |S − S^T| = {asym:.3e}", flush=True)
    print(f"max |S_ba − S_ab^T| = {np.max(np.abs(S_ba - S_ab.T)):.3e}",
          flush=True)

    # correlation (normalize by own-metric diagonals)
    d = np.sqrt(np.clip(np.abs(np.diag(S)), 1e-30, None))
    Shat = S / np.outer(d, d)
    print("\n=== Ŝ_ij = S_ij / sqrt(|S_ii S_jj|) ===", flush=True)
    print(fmt_mat(Shat, names), flush=True)

    out = {
        "names": names,
        "S": S.tolist(),
        "S_hat": Shat.tolist(),
        "S_aa": S_aa.tolist(),
        "S_bb": S_bb.tolist(),
        "S_aa_leftover": S_aa_leftover.tolist(),
        "S_bb_leftover": S_bb_leftover.tolist(),
        "S_ab": S_ab.tolist(),
        "S_ba": S_ba.tolist(),
        "max_aa_off_I": float(np.max(np.abs(S_aa - np.eye(n)))),
        "max_bb_off_I": float(np.max(np.abs(S_bb - np.eye(n)))),
        "max_asym": float(asym),
        "sanity_cross_vs_N": float(max_d),
        "tdm_3344": float(v3344),
        "tdm_4433": float(v4433),
        "load_ritz": load,
        "E_gs": e_gs,
        "E_x": e_x,
    }
    path = eom_result_path(nucleus, 16, int(os.environ.get("EMAX", "3")),
                           "mr0_overlap")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nwrote {path}", flush=True)


if __name__ == "__main__":
    main()
