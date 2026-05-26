import sys
sys.path.insert(0, '../build')
from pyIMSRG import *

# Build modelspace from p-shell orbits (same as he8.ref)
ms = ModelSpace()
ms.SetEmax3Body(1000)
orbits = [(0,1,3,-1),(0,1,1,-1),(0,1,3,1),(0,1,1,1)]
for n,l,j2,tz2 in orbits:
    ms.AddOrbit(n, l, j2, tz2, 0.0, 1)
emax = max(2*n+l for n,l,j2,tz2 in orbits)
ms.SetEmax(emax)
ms.SetE2max(2*emax)
ms.SetE3max(3*emax)
ms.SetupKets()
ms.Setup3bKets()

nch = ms.GetNumberThreeBodyChannels()
print(f"Channels: {nch}")

# Build operator, set each element to a unique known value
op = Operator(ms, 0, 0, 0, 3)
op.ThreeBody.SetMode("pn")

ref = {}  # (ich, ib, ik) -> value
idx = 0
for ich in range(nch):
    ch = ms.GetThreeBodyChannel(ich)
    nk = ch.GetNumber3bKets()
    for ib in range(nk):
        for ik in range(ib, nk):
            val = float(idx + 1) * 0.123456
            op.ThreeBody.SetME_pn_ch(ich, ich, ib, ik, val)
            ref[(ich, ib, ik)] = val
            idx += 1

print(f"Elements set: {idx}")

# Test 1: GetME_pn_ch upper triangle
errors = 0
for (ich, ib, ik), expected in ref.items():
    got = op.ThreeBody.GetME_pn_ch(ich, ich, ib, ik)
    if abs(got - expected) > 1e-10:
        print(f"  FAIL upper: ch={ich} ib={ib} ik={ik}  expected={expected:.6f}  got={got:.6f}")
        errors += 1
print(f"Test 1 (upper triangle): {'PASS' if errors==0 else f'FAIL ({errors} errors)'}")

# Test 2: GetME_pn_ch lower triangle (symmetry: M[ik,ib] == M[ib,ik] for scalar hermitian)
errors2 = 0
for (ich, ib, ik), expected in ref.items():
    if ib == ik:
        continue
    got = op.ThreeBody.GetME_pn_ch(ich, ich, ik, ib)  # swapped
    if abs(got - expected) > 1e-10:
        print(f"  FAIL lower: ch={ich} ib={ib} ik={ik}  expected={expected:.6f}  got={got:.6f}")
        errors2 += 1
print(f"Test 2 (lower triangle symmetry): {'PASS' if errors2==0 else f'FAIL ({errors2} errors)'}")

# Test 3: GetME_pn with orbit labels - iterate all channels and kets,
# look up orbit indices and call GetME_pn, compare to GetME_pn_ch
errors3 = 0
for ich in range(nch):
    ch = ms.GetThreeBodyChannel(ich)
    nk = ch.GetNumber3bKets()
    twoJ = ch.twoJ
    for ib in range(nk):
        ket_b = ch.GetKet(ib)
        a, b, c = ket_b.p, ket_b.q, ket_b.r
        Jab = ket_b.Jpq   # integer Jab
        for ik in range(ib, nk):
            ket_k = ch.GetKet(ik)
            d, e, fk = ket_k.p, ket_k.q, ket_k.r
            Jde = ket_k.Jpq
            expected = op.ThreeBody.GetME_pn_ch(ich, ich, ib, ik)
            got = op.ThreeBody.GetME_pn(Jab, Jde, twoJ, a, b, c, d, e, fk)
            if abs(got - expected) > 1e-10:
                print(f"  FAIL GetME_pn: ch={ich} ({a},{b},{c},Jab={Jab}) ({d},{e},{fk},Jde={Jde}) twoJ={twoJ}")
                print(f"    expected={expected:.6f}  got={got:.6f}")
                errors3 += 1
print(f"Test 3 (GetME_pn vs GetME_pn_ch): {'PASS' if errors3==0 else f'FAIL ({errors3} errors)'}")

print()
print("ALL PASS" if errors+errors2+errors3==0 else "SOME TESTS FAILED")
