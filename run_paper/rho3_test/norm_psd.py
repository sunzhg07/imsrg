#!/usr/bin/env python3
"""Spectrum of the MR-EOM norm kernel, rho3 off and on.
  psd.py EMAX REF CORE VAL E3MAX J2 PAR ITZ [DTZ]"""
import os
import sys, os
import numpy as np
ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "build"))
from pyIMSRG import EOM, ModelSpace, Operator
emax=int(sys.argv[1]); ref=sys.argv[2]; core=sys.argv[3]; val=sys.argv[4]; e3=int(sys.argv[5])
J2,par,itz=(int(x) for x in sys.argv[6:9]); dtz=int(sys.argv[9]) if len(sys.argv)>9 else -itz
tensor=(J2,par,itz)!=(0,0,0)
ms=ModelSpace(emax,core,val)
if e3>0: ms.SetE3max(e3)
Hs=Operator(ms,0,0,0,2)
eom=EOM(Hs,ref,J2,par,itz)
if itz: eom.SetDeltaTz(dtz)
(eom.ConstructConfigs_tensor if tensor else eom.ConstructConfigs)()
for flag in (False,True):
    eom.SetUseRdm3(flag)
    (eom.ConstructNormMatrix_tensor if tensor else eom.ConstructNormMatrix)()
    N=np.array(eom.Nkernel.todense() if hasattr(eom.Nkernel,"todense") else eom.Nkernel)
    asym=abs(N-N.T).max() if N.size else 0.0
    w=np.linalg.eigvalsh(0.5*(N+N.T)) if N.size else np.zeros(0)
    tol=1e-8*max(1.0,abs(w).max() if w.size else 1.0)
    print("PSD J2=%d P=%d itz=%d rho3=%d dim=%d min=%+.3e max=%.3e nneg=%d nnull=%d asym=%.1e -> %s"%(
        J2,par,itz,flag,len(w),w.min() if w.size else 0,w.max() if w.size else 0,
        int((w<-tol).sum()),int((abs(w)<=tol).sum()),asym,"PSD" if (w>=-tol).all() else "NOT PSD"),flush=True)
