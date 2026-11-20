#!/usr/bin/env python3
"""Level schemes: MR-EOM (rho3 metric) vs full-space FCI, absolute energies.

  NUCLEUS=He8|Be8 python3 plot_mr_rho3.py

Reads eom_result/{Nucleus}_*_rho3_full.json (from mr_rho3_full.py) and the
KSHELL FCI logs log_{Nucleus}_fci_heb_{hf,hs}_n03_{Nucleus}_j{0,4}p.txt.
Columns: HF FCI (unevolved H), Hs FCI (IMSRG-evolved H), p-shell SM,
SR-EOM (He8 only), MR-EOM. FCI energies get the .snt zero-body term added.
Same line/marker conventions as plot_he8_mr_eom.py; energies are absolute.
"""
from __future__ import annotations

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from he8_io import EOM_RESULT, PAPER, eom_result_path
from plot_he8_mr_eom import (FCI_DIR, classify_ph, draw_fci_column, hlines, href,
                             load_kshell_jlog)

nucleus = os.environ.get("NUCLEUS", "He8")
emax, hw = 3, 16
E_WINDOW = float(os.environ.get("E_WINDOW", "18"))  # MeV above the SM 0+ gs

JOBS = {
    0: (("mr0gs", r"MR $0^+_{\mathrm{gs}}{\to}0^+$", "#c1121f"),
        ("mr0x",  r"MR $0^+_2{\to}0^+$",             "#e09f3e"),
        ("t20",   r"MR $2^+_1{\to}0^+$",             "#1d3557"),
        ("t20x",  r"MR $2^+_2{\to}0^+$",             "#48cae4")),
    2: (("mr2gs", r"MR $2^+_1{\to}2^+$",             "#6a4c93"),
        ("mr2x",  r"MR $2^+_2{\to}2^+$",             "#c77dff"),
        ("t02",   r"MR $0^+_{\mathrm{gs}}{\to}2^+$", "#2d6a4f"),
        ("t02x",  r"MR $0^+_2{\to}2^+$",             "#74c69d")),
}


def fci_zero_body(ham):
    """KSHELL does not add the .snt zero-body term; read it from the header."""
    with open(os.path.join(PAPER, f"fci_heb_{ham}_n03_{nucleus}.snt")) as f:
        for line in f:
            if "Zero body term" in line:
                return float(line.split(":")[1])
    raise RuntimeError("no zero-body term in FCI .snt")


def fci(ham, J):
    name = f"log_{nucleus}_fci_heb_{ham}_n03_{nucleus}_j{2 * J}p.txt"
    zb = fci_zero_body(ham)
    for d in (FCI_DIR, PAPER):
        p = os.path.join(d, name)
        if os.path.isfile(p):
            sts = load_kshell_jlog(p, j_want=J)
            for st in sts:
                st["E"] += zb  # absolute energy
            if nucleus != "He8":  # p-h rank relative to (0s)^4 (0p)^4
                for st in sts:
                    n0p = st["n_n0p"] + st["n_p0p"]
                    st["n_ph"] = (4.0 - st["n_0s"]) + max(0.0, 4.0 - n0p)
                    st["kind"], st["ph_lab"] = classify_ph(st["n_ph"])
            if sts:
                return sts
    return []


def main():
    with open(eom_result_path(nucleus, hw, emax, "rho3_full")) as f:
        data = json.load(f)
    ks, runs = data["kshell"], data["runs"]
    e_core = ks["E_core"]
    sm = {J: [s["E_sm"] + e_core for s in ks["he8"] if s["twoJ"] == 2 * J] for J in (0, 2)}
    e_gs = sm[0][0]
    sr = {}
    if nucleus == "He8":
        for J, nm in ((0, "sr0"), (2, "sr2")):
            p = eom_result_path(nucleus, hw, emax, nm)
            if os.path.isfile(p):
                sr[J] = json.load(open(p))

    fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.6), sharey=True)
    lo = e_gs
    for ax, J in zip(axes, (0, 2)):
        cols = []  # (label, draw function)
        hf = [s for s in fci("hf", J)]
        hs = [s for s in fci("hs", J)]
        top = e_gs + E_WINDOW
        x = 0.0
        xs, labels = [], []
        for lab, sts, color in ((rf"HF FCI ${J}^+$" "\n(unevolved)", hf, "#868e96"),
                                (rf"Hs FCI ${J}^+$" "\n(IMSRG)", hs, "#212529")):
            if not sts:
                continue
            keep = [s for s in sts if s["E"] <= top]
            draw_fci_column(ax, x, keep, [s["E"] for s in keep], color)
            lo = min(lo, keep[0]["E"]) if keep else lo
            xs.append(x); labels.append(lab); x += 1.0
        hlines(ax, x, [e for e in sm[J] if e <= top], "#212529", lw=2.6)
        xs.append(x); labels.append("p-shell SM\n" rf"${J}^+$"); x += 1.0
        if J in sr:
            rec = sr[J]
            href(ax, x, rec["eref"], "#2a6f97")
            hlines(ax, x, [e for e in rec["E_abs"] if e <= top], "#2a6f97", lw=2.15)
            xs.append(x); labels.append("SR-EOM\n" rf"${J}^+$"); x += 1.0
        for name, lab, color in JOBS[J]:
            rec = runs.get(f"{name}_rdm3")
            if not rec or "error" in rec:
                continue
            if not any(e <= top for e in rec["E_abs"]):
                continue  # nothing inside the plotted window
            href(ax, x, rec["E_ref"], color)
            res = rec.get("residuals") or [0.0] * len(rec["E_abs"])
            for e, r in zip(rec["E_abs"], res):
                if e > top:
                    continue
                ok = abs(r) < 0.1
                ax.hlines(e, x - 0.32, x + 0.32, colors=color, linewidth=2.15 if ok else 1.1,
                          linestyles="-" if ok else (0, (1, 1)), zorder=3)
            xs.append(x); labels.append(lab); x += 1.0
        ax.set_xticks(xs)
        ax.set_xticklabels(labels, fontsize=7.6)
        ax.set_xlim(-0.6, x - 0.4)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.axhline(e_gs, color="0.8", lw=0.6, zorder=1)
        ax.set_title(rf"$J^\pi={J}^+$", fontsize=10)
    axes[0].set_ylabel("E (MeV)")
    axes[0].set_ylim(lo - 1.0, e_gs + E_WINDOW + 0.5)
    fig.legend(handles=[
        Line2D([0], [0], color="0.15", lw=2.8, marker="D", ms=4, label=r"FCI $0\mathrm{p}0\mathrm{h}$ (p-shell)"),
        Line2D([0], [0], color="0.15", lw=1.7, marker="o", ms=3.2, markerfacecolor="white", label=r"FCI $1\mathrm{p}1\mathrm{h}$"),
        Line2D([0], [0], color="0.15", lw=1.6, marker="s", ms=3.0, markerfacecolor="white", label=r"FCI $2\mathrm{p}2\mathrm{h}$"),
        Line2D([0], [0], color="0.15", lw=1.4, linestyle=(0, (2.2, 1.4)), marker="^", ms=3.0, markerfacecolor="white", label=r"FCI $3\mathrm{p}3\mathrm{h}+$"),
        Line2D([0], [0], color="0.35", lw=1.25, linestyle=(0, (3, 2)), label=r"reference: SR $E_0$ / MR $\langle H\rangle_\rho$"),
        Line2D([0], [0], color="0.35", lw=2.15, label="EOM root, residual < 0.1 MeV"),
        Line2D([0], [0], color="0.35", lw=1.1, linestyle=(0, (1, 1)), label="EOM root, not converged"),
    ], frameon=False, loc="lower center", fontsize=8.0, ncol=7)
    A = {"He8": r"$^{8}$He", "Be8": r"$^{8}$Be"}[nucleus]
    fig.suptitle(rf"{A}  $e_{{\max}}=3$  EM1.8/2.0  $\hbar\omega=16$ MeV   "
                 r"MR-EOM (H2+H3, $\rho_3$ in norm) vs full-space FCI", fontsize=11)
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    stem = os.path.join(EOM_RESULT, f"{nucleus}_EM1.8_2.0_hw16_emax3_rho3_levels")
    fig.savefig(stem + ".png", dpi=180)
    fig.savefig(stem + ".pdf")
    print("wrote", stem + ".png")
    print("wrote", stem + ".pdf")

    # numeric comparison: nearest FCI level to each converged MR root
    for J in (0, 2):
        for ham in ("hs", "hf"):
            sts = fci(ham, J)
            if not sts:
                continue
            print(f"\n{nucleus} {J}+  nearest {ham.upper()} FCI level to each MR root")
            for name, lab, _ in JOBS[J]:
                rec = runs.get(f"{name}_rdm3")
                if not rec or "error" in rec:
                    continue
                for e, r in list(zip(rec["E_abs"], rec.get("residuals", [])))[:4]:
                    st = min(sts, key=lambda s: abs(s["E"] - e))
                    print(f"  {name:6s} E={e:9.3f} (resid {r:.3f})  FCI #{st['idx']:2d} "
                          f"{st['E']:9.3f} [{st['kind']}]  diff {e - st['E']:+.3f}")


if __name__ == "__main__":
    main()
