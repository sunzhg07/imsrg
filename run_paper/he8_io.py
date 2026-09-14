"""Shared He8 paper I/O: Hs cache names from calculation parameters.

  hs_cache/Hs_{nucleus}_{LECs}_hw{hw}_emax{emax}.dat          closed-shell SR
  hs_cache/Hs_{nucleus}_{val}_{LECs}_hw{hw}_emax{emax}.dat    VS (val ≠ nucleus)

  eom_result/{nucleus}_{LECs}_hw{hw}_emax{emax}_{job}.json
"""
from __future__ import annotations

import os

PAPER = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(PAPER, ".."))
RUN2 = os.path.join(ROOT, "run2")
BUILD = os.path.join(ROOT, "build")
HS_CACHE = os.path.join(ROOT, "hs_cache")
EOM_RESULT = os.path.join(ROOT, "eom_result")

LECS_DEFAULT = "EM1.8_2.0"


def lecs():
    return os.environ.get("LECS", LECS_DEFAULT)


def _val_tag(val: str | None) -> str | None:
    if not val:
        return None
    return val.replace("-", "").replace("_", "")


def hs_filename(nucleus: str, hw: float | int, emax: int, val: str | None = None) -> str:
    """Hs name from nucleus, (optional valence), LECs, ħω, emax."""
    parts = ["Hs", nucleus]
    tag = _val_tag(val)
    if tag and tag != nucleus:
        parts.append(tag)
    parts += [lecs(), f"hw{int(hw)}", f"emax{int(emax)}"]
    return "_".join(parts) + ".dat"


def hs_path(nucleus: str, hw: float | int, emax: int, val: str | None = None) -> str:
    os.makedirs(HS_CACHE, exist_ok=True)
    return os.path.join(HS_CACHE, hs_filename(nucleus, hw, emax, val))


def eom_result_filename(nucleus: str, hw: float | int, emax: int, job: str) -> str:
    return f"{nucleus}_{lecs()}_hw{int(hw)}_emax{int(emax)}_{job}.json"


def eom_result_path(nucleus: str, hw: float | int, emax: int, job: str) -> str:
    os.makedirs(EOM_RESULT, exist_ok=True)
    return os.path.join(EOM_RESULT, eom_result_filename(nucleus, hw, emax, job))


def fci_snt_filename(nucleus: str, hw: float | int, emax: int,
                     tag: str | None = None) -> str:
    """Vacuum FCI .snt after UndoNormalOrdering.

    tag=None → Hs_{nuc}_FCI_...          (VS / MR)
    tag='SR' → Hs_{nuc}_SR_FCI_...       (closed-shell SR)
    """
    parts = ["Hs", nucleus]
    if tag:
        parts.append(tag)
    parts += ["FCI", lecs(), f"hw{int(hw)}", f"emax{int(emax)}"]
    return "_".join(parts) + ".snt"


def fci_snt_path(nucleus: str, hw: float | int, emax: int,
                 tag: str | None = None) -> str:
    os.makedirs(HS_CACHE, exist_ok=True)
    return os.path.join(HS_CACHE, fci_snt_filename(nucleus, hw, emax, tag))


def pshell_snt_filename(nucleus: str, hw: float | int, emax: int) -> str:
    """He4-core p-shell effective interaction for KSHELL."""
    return f"Hs_{nucleus}_pshell_{lecs()}_hw{int(hw)}_emax{int(emax)}.snt"


def pshell_snt_path(nucleus: str, hw: float | int, emax: int) -> str:
    os.makedirs(HS_CACHE, exist_ok=True)
    return os.path.join(HS_CACHE, pshell_snt_filename(nucleus, hw, emax))
