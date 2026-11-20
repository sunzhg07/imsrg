#!/usr/bin/env python3
"""MR-EOM norm kernel vs commutator gold (scalar and tensor, rho3 off and on).

For random amplitudes confined to a config sub-block, compare
  kernel:     ComputeNorm(a, b)  = a^T Nkernel b   (ConstructNormMatrix[_tensor])
  commutator: NormMultiref(a, b) = <[S+_a, S-_b]>_rho / 2  (comm223ss/tts for 3b)
on the FULL values (rho3 off and rho3 on), for a == b and for a != b
(bilinear, catches antisymmetric cross terms a quadratic form cannot see).
Blocks include the tensor "swapped" configs (creator pair in the stored ket).

  python3 test_rho3_norm.py EMAX REF_FILE [CORE VALENCE E3MAX J2 PARITY ITZ [DTZ]]
  DTZ: signed nucleus Tz change for itz != 0 (default -itz, n -> p)
"""
import os
import sys
import time

import numpy as np

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "build"))
from pyIMSRG import EOM, ModelSpace, Operator  # noqa: E402

emax = int(sys.argv[1])
ref = sys.argv[2]
core = sys.argv[3] if len(sys.argv) > 3 else "O16"
val = sys.argv[4] if len(sys.argv) > 4 else "sd-shell"
e3max = int(sys.argv[5]) if len(sys.argv) > 5 and int(sys.argv[5]) > 0 else None
J2, par, itz = (int(x) for x in sys.argv[6:9]) if len(sys.argv) > 8 else (0, 0, 0)
dtz = int(sys.argv[9]) if len(sys.argv) > 9 else -itz
tensor = (J2, par, itz) != (0, 0, 0)

ms = ModelSpace(emax, core, val)
if e3max is not None:
    ms.SetE3max(e3max)
Hs = Operator(ms, 0, 0, 0, 2)
eom = EOM(Hs, ref, J2, par, itz)
if itz:
    eom.SetDeltaTz(dtz)
if tensor:
    eom.ConstructConfigs_tensor()
else:
    eom.ConstructConfigs()
confs = eom.eom_confs


def build_norm():
    if tensor:
        eom.ConstructNormMatrix_tensor()
    else:
        eom.ConstructNormMatrix()


def cvq(o):
    return ms.GetOrbit(o).cvq


def kind(k):
    return "".join("cvq"[cvq(o)] for o in (k.p, k.q))


def sub(start, end, dim, test):
    """configs in [start,end] with test(creator kind, annihilator kind, swapped)"""
    if dim == 0:
        return []
    out = []
    for i in range(start, end + 1):
        b, k, chb, chk = confs[i]
        if not tensor:
            chk = chb
        kb = ms.GetTwoBodyChannel(chb).GetKet(b)
        kk = ms.GetTwoBodyChannel(chk).GetKet(k)
        o = eom.StoredExcitationOrientation(chb, chk, b, k) if tensor else 1
        cre, ann = (kind(kb), kind(kk)) if o >= 0 else (kind(kk), kind(kb))
        if test(cre, ann, o < 0):
            out.append(i)
    return out


def is_qv(s):
    return sorted(s) == ["q", "v"]


P = (eom.pphv_start, eom.pphv_end, eom.pphv_dim)
V = (eom.ppvv_start, eom.ppvv_end, eom.ppvv_dim)
H = (eom.pphh_start, eom.pphh_end, eom.pphh_dim)
blocks = {
    "pphv vv<-cv":    sub(*P, lambda c, a, s: c == "vv" and a in ("cv", "vc")),
    "pphv main":      sub(*P, lambda c, a, s: not s),
    "pphv swapped":   sub(*P, lambda c, a, s: s),
    "ppvv qv<-vv":    sub(*V, lambda c, a, s: is_qv(c)),
    "ppvv vv<-vv":    sub(*V, lambda c, a, s: c == "vv"),
    "ppvv swapped":   sub(*V, lambda c, a, s: s),
    "pphh swapped":   sub(*H, lambda c, a, s: s),
    "all": list(range(eom.eom_dims)),
}
blocks = {"all": blocks["all"]}  # quick: one random vector pair over all configs
print("J2=%d parity=%d itz=%d dtz=%d  eom_dims %d" % (J2, par, itz, dtz if itz else 0, eom.eom_dims),
      {k: len(v) for k, v in blocks.items()}, flush=True)

rng = np.random.default_rng(3)
vecs = {}
for name, ids in blocks.items():
    vs = []
    for _ in range(2):
        v = np.zeros(eom.eom_dims)
        v[ids] = rng.normal(size=len(ids))
        vs.append(v)
    vecs[name] = vs


def as_op(v):
    op = Operator(ms, J2, itz, par, 2)
    eom.UnflattenOperator(op, v)
    return eom.GetVSEOM_ladder_multiref(op, 1)


res = {}
for flag in (False, True):
    t0 = time.time()
    eom.SetUseRdm3(flag)
    build_norm()
    for name, vs in vecs.items():
        ops = [as_op(v) for v in vs]
        pairs = [(0, 0), (0, 1), (1, 0)]
        for (s, t) in pairs:
            res[(flag, name, s, t)] = (eom.ComputeNorm(ops[s], ops[t]),
                                       eom.NormMultiref(ops[s], ops[t]))
    print("use_rdm3=%s done in %.1fs" % (flag, time.time() - t0), flush=True)

worst = 0.0
for name in vecs:
    for (s, t) in [(0, 0), (0, 1), (1, 0)]:
        line = []
        for flag in (False, True):
            c, m = res[(flag, name, s, t)]
            d = c - m
            scale = max(abs(c), abs(m), 1e-12)
            rel = abs(d) / scale if scale > 1e-10 else 0.0
            worst = max(worst, rel)
            line.append("rho3=%d kernel %+.8f commut %+.8f rel %.1e" % (flag, c, m, rel))
        print("%-14s v%d,v%d  %s" % (name, s, t, " | ".join(line)))
print("max relative diff (full values): %.2e  -> %s" % (worst, "PASS" if worst < 1e-8 else "FAIL"))
