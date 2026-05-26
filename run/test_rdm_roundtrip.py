#!/usr/bin/env python3
"""
Round-trip test: read he8.ref with EOM(H, file), write it back with WriteTdm,
read the written file again, and verify all 3b matrix elements are identical.
Usage:  python test_rdm_roundtrip.py
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build'))
from pyIMSRG import *

os.chdir(os.path.dirname(os.path.abspath(__file__)))
TDM_FILE   = 'he8.ref'
COPY_FILE  = '/tmp/rdm_roundtrip_copy.dat'
TOL = 1e-8

# ----------------------------------------------------------------
# 1. Build a minimal Hs so EOM constructor is satisfied
#    (the string constructor builds its own rdm_modelspace from the file)
# ----------------------------------------------------------------
ms_hs = ModelSpace(1, 'He4', 'p-shell')
ms_hs.SetHbarOmega(20.0)
H = Operator(ms_hs, 0, 0, 0, 2)

# ----------------------------------------------------------------
# 2. Read he8.ref -> eom1.rdm gets its own rdm_modelspace from file
# ----------------------------------------------------------------
eom1 = EOM(H, TDM_FILE, 0, 0, 0)
op1  = eom1.rdm
ms1  = op1.GetModelSpace()
norb = ms1.GetNumberOrbits()
nch3 = ms1.GetNumberThreeBodyChannels()
print(f'Read {TDM_FILE}: {norb} orbits, {nch3} 3-body channels')

# ----------------------------------------------------------------
# 3. Write op1 back to file using WriteTdm
# ----------------------------------------------------------------
eom1.WriteTdm(op1, COPY_FILE)
print(f'Written to {COPY_FILE}')

# ----------------------------------------------------------------
# 4. Read the written file back
# ----------------------------------------------------------------
eom2 = EOM(H, COPY_FILE, 0, 0, 0)
op2  = eom2.rdm
ms2  = op2.GetModelSpace()

# ----------------------------------------------------------------
# 5. Compare every injected element
# ----------------------------------------------------------------
# 5. Compare all 3b elements (ch,ibra,iket) in native order
# ----------------------------------------------------------------
max_diff = 0.0
n_compared = 0
n_mismatch = 0

for ch_bra in range(nch3):
    TBC_bra = ms1.GetThreeBodyChannel(ch_bra)
    nk_bra  = TBC_bra.GetNumber3bKets()
    for ch_ket in range(ch_bra, nch3):
        TBC_ket = ms1.GetThreeBodyChannel(ch_ket)
        nk_ket  = TBC_ket.GetNumber3bKets()
        iket0   = 0
        for ibra in range(nk_bra):
            if ch_bra == ch_ket:
                iket0 = ibra
            for iket in range(iket0, nk_ket):
                v1   = op1.ThreeBody.GetME_pn_ch(ch_bra, ch_ket, ibra, iket)
                v2   = op2.ThreeBody.GetME_pn_ch(ch_bra, ch_ket, ibra, iket)
                diff = abs(v1 - v2)
                max_diff = max(max_diff, diff)
                n_compared += 1
                if diff > TOL:
                    n_mismatch += 1
                    if n_mismatch <= 5:
                        kb = TBC_bra.GetKet(ibra)
                        kk = TBC_ket.GetKet(iket)
                        print(f'  MISMATCH ch_bra={ch_bra} ch_ket={ch_ket} '
                              f'ib={ibra} ik={iket}  '
                              f'bra=({kb.p},{kb.q},{kb.r},J={kb.Jpq}) '
                              f'ket=({kk.p},{kk.q},{kk.r},J={kk.Jpq})  '
                              f'orig={v1:.8f} back={v2:.8f} diff={diff:.2e}')

print(f'Compared {n_compared} 3b elements')
if n_mismatch == 0:
    print(f'  ALL match  (max_diff={max_diff:.2e})  PASS')
else:
    print(f'  {n_mismatch}/{n_compared} mismatches, max_diff={max_diff:.2e}  FAIL')
    sys.exit(1)
