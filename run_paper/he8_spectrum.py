#!/usr/bin/env python3
"""He8 emax=3 spectrum figure: SR / MR scalar / MR tensor vs p-shell FCI.

Reads JSON from eom_result/{nucleus}_{LECs}_hw{hw}_emax{emax}_{job}.json
with fallbacks to run_paper/results and run2/he8_e3_results.

From run_paper/:
  PYTHONPATH=../build python3 he8_spectrum.py
"""
from __future__ import annotations

import json
import os

from he8_io import EOM_RESULT, PAPER, RUN2, eom_result_path

HE8 = os.path.join(RUN2, "he8")
FALLBACKS = (
    os.path.join(PAPER, "results"),
    os.path.join(RUN2, "he8_e3_results"),
)

emax = int(os.environ.get("EMAX", "3"))
hw = 16
nucleus = "He8"

FCI_0 = [-23.004707, -13.300021]
JOBS = {
    "sr0": {"lab": r"SR $0^+$", "kind": "sr0"},
    "sr2": {"lab": r"SR $2^+$", "kind": "sr2"},
    "mr0gs": {"lab": r"MR $0^+_{\mathrm{gs}}{\to}0^+$", "kind": "mr0"},
    "mr0x": {"lab": r"MR $0^+_1{\to}0^+$", "kind": "mr0x"},
    "mr2gs": {"lab": r"MR $2^+_{\mathrm{gs}}{\to}2^+$", "kind": "mr2"},
    "mr2x": {"lab": r"MR $2^+_1{\to}2^+$", "kind": "mr2x"},
    "t02": {"lab": r"MR $0^+_{\mathrm{gs}}{\to}2^+$", "kind": "t02"},
    "t20": {"lab": r"MR $2^+_{\mathrm{gs}}{\to}0^+$", "kind": "t20"},
}


def job_json(name):
    path = eom_result_path(nucleus, hw, emax, name)
    if os.path.isfile(path):
        return path
    for root in FALLBACKS:
        cand = os.path.join(root, f"{name}.json")
        if os.path.isfile(cand):
            return cand
    return path


def load_rec(name):
    path = job_json(name)
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        return json.load(f)


def load_fci_0():
    path = os.path.join(HE8, "he8_0p_fci.json")
    if os.path.isfile(path):
        with open(path) as f:
            d = json.load(f)
        return [float(x) for x in d.get("E_abs", [])]
    return list(FCI_0)


def load_fci_2():
    path = os.path.join(HE8, "he8_2p_fci.json")
    if not os.path.isfile(path):
        return [-17.731203, -11.727562]
    with open(path) as f:
        d = json.load(f)
    return [float(x) for x in d.get("E_abs", [])]


def hlines(ax, x, energies, color, lw=2.2, half=0.32):
    for e in energies:
        ax.hlines(e, x - half, x + half, colors=color, linewidth=lw, zorder=3)


def href(ax, x, e, color, half=0.36):
    ax.hlines(e, x - half, x + half, colors=color, linewidth=1.3,
              linestyles=(0, (3, 2)), zorder=2)


def plot():
    import matplotlib.pyplot as plt

    os.makedirs(EOM_RESULT, exist_ok=True)
    fci0 = load_fci_0()
    fci2 = load_fci_2()
    z0 = fci0[0] if fci0 else FCI_0[0]
    colors = {
        "fci0": "#222222",
        "fci2": "#6c757d",
        "sr0": "#2a6f97",
        "sr2": "#90e0ef",
        "mr0": "#c1121f",
        "mr0x": "#e09f3e",
        "mr2": "#6a4c93",
        "mr2x": "#c77dff",
        "t02": "#2d6a4f",
        "t20": "#95d5b2",
    }

    fig, ax = plt.subplots(figsize=(13.2, 5.8))
    columns = [("fci", None)] + [(n, JOBS[n]) for n in
                                 ("sr0", "sr2", "mr0gs", "mr0x",
                                  "mr2gs", "mr2x", "t02", "t20")]
    xs, labels = [], []
    t02_note = ""
    for i, (name, spec) in enumerate(columns):
        x = float(i)
        xs.append(x)
        if name == "fci":
            labels.append("FCI\n(p-shell)")
            hlines(ax, x, [e - z0 for e in fci0], colors["fci0"], lw=2.5,
                   half=0.22)
            if fci2:
                hlines(ax, x, [e - z0 for e in fci2], colors["fci2"], lw=2.5,
                       half=0.22)
            continue
        rec = load_rec(name)
        lab = spec["lab"]
        if rec is not None and rec.get("steps"):
            if name in ("t02", "t20"):
                lab = spec["lab"] + f"\n({rec['steps']} steps)"
            if name == "t02":
                reason = rec.get("stop_reason", "")
                t02_note = f"tensor $0^+\\to 2^+$ {rec['steps']} steps"
                if reason:
                    t02_note += f" ({reason.replace('_', ' ')})"
        labels.append(lab)
        if rec is None:
            ax.text(x, 0.5, "—", ha="center", va="bottom", fontsize=9,
                    color="0.6")
            continue
        color = colors[spec["kind"]]
        href(ax, x, rec["eref"] - z0, color)
        hlines(ax, x, [e - z0 for e in rec["E_abs"]], color, lw=2.15)

    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel(r"$E_x$ vs p-shell FCI $0^+_{\mathrm{gs}}$ (MeV)")
    title = (r"$^8$He  $e_{\mathrm{max}}=3$  EM1.8/2.0  "
             r"$\hbar\omega=16$  NO2B" + "\n5 EOM states")
    if t02_note:
        title += "; " + t02_note
    ax.set_title(title)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(-0.6, len(columns) - 0.4)
    ax.axhline(0.0, color="0.75", lw=0.6, zorder=1)
    ax.legend(
        handles=[
            plt.Line2D([0], [0], color=colors["fci0"], lw=2.4, label=r"FCI $0^+$"),
            plt.Line2D([0], [0], color=colors["fci2"], lw=2.4, label=r"FCI $2^+$"),
            plt.Line2D([0], [0], color=colors["sr0"], lw=2, label=r"SR-EOM $0^+$"),
            plt.Line2D([0], [0], color=colors["sr2"], lw=2, label=r"SR-EOM $2^+$"),
            plt.Line2D([0], [0], color=colors["mr0"], lw=2,
                       label=r"MR scalar $0^+{\to}0^+$"),
            plt.Line2D([0], [0], color=colors["mr2"], lw=2,
                       label=r"MR scalar $2^+{\to}2^+$"),
            plt.Line2D([0], [0], color=colors["t02"], lw=2, label=r"MR tensor"),
            plt.Line2D([0], [0], color="0.35", lw=1.3, linestyle=(0, (3, 2)),
                       label="SR $E_0$ / MR $\\langle H\\rangle_\\rho$"),
        ],
        frameon=False, loc="upper left", fontsize=7.5, ncol=2,
    )
    fig.tight_layout()
    png = os.path.join(EOM_RESULT, "he8_e3_spectrum.png")
    pdf = os.path.join(EOM_RESULT, "he8_e3_spectrum.pdf")
    fig.savefig(png, dpi=160)
    fig.savefig(pdf)
    print(f"wrote {png}")
    print(f"wrote {pdf}")
    missing = [n for n in JOBS if load_rec(n) is None]
    if missing:
        print("missing jobs:", ", ".join(missing))


if __name__ == "__main__":
    plot()
