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

# EM1.8/2.0 hw=16 interaction files: (path, file truncations, 3N storage mode).
# Local (Mac) files first; otherwise the OLCF proj-shared copies (riker/andes),
# which carry larger file truncations and half-precision NO2B 3N.
_ME_LOCAL = "/Users/wolf/work/srg_io/input"
_ME_OLCF = os.path.expanduser("~/proj-share/nuccor_interactions/me_takayuki")
_INTERACTIONS = (
    (os.path.join(_ME_LOCAL, "TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax14_e2max28.me2j.gz"),
     (14, 28, 14),
     os.path.join(_ME_LOCAL, "NO2B_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms18_36_18.stream.bin"),
     (18, 36, 18), "no2b"),
    (os.path.join(_ME_OLCF, "me2j", "TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16_emax16_e2max32.me2j.gz"),
     (16, 32, 16),
     os.path.join(_ME_OLCF, "me3j", "NO2B_half_ThBME_EM1.8_2.0_3NFJmax15_IS_hw16_ms18_36_28.stream.bin"),
     (18, 36, 28), "no2bhalf"),
)


def read_em18_hw16(H, rw):
    """Read 2B + NO2B 3N (EM1.8/2.0, hw=16) into H from whichever copy exists."""
    for f2b, t2b, f3b, t3b, mode in _INTERACTIONS:
        if os.path.isfile(f2b) and os.path.isfile(f3b):
            print(f"2B: {f2b} {t2b}\n3N: {f3b} {t3b} mode={mode}", flush=True)
            rw.ReadBareTBME_Darmstadt(f2b, H, *t2b)
            H.ThreeBody.SetMode(mode)
            H.ThreeBody.ReadFile([f3b], list(t3b))
            return
    raise FileNotFoundError("EM1.8/2.0 hw16 interaction files not found: "
                            + ", ".join(i[0] for i in _INTERACTIONS))


# KSHELL executables: env override, else the repo-local kshell_he8/, else the
# riker install.
_KSHELL_DIRS = (os.path.join(PAPER, "kshell_he8"),
                os.path.expanduser("~/codes/kshell-20230714-ver4/bin"))


def kshell_path(name: str) -> str:
    """kshell.exe / transit.exe / gen_partition.py location."""
    env = {"kshell.exe": "KSHELL_EXE", "transit.exe": "TRANSIT_EXE"}.get(name)
    if env and os.environ.get(env):
        return os.environ[env]
    if os.environ.get("KSHELL_BIN"):
        return os.path.join(os.environ["KSHELL_BIN"], name)
    for d in _KSHELL_DIRS:
        if os.path.isfile(os.path.join(d, name)):
            return os.path.join(d, name)
    return os.path.join(_KSHELL_DIRS[-1], name)


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
