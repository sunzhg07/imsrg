import sys, math
sys.path.insert(0, '../build')
import pyIMSRG as imsrg

# Build minimal Hs on p-shell modelspace
ms = imsrg.ModelSpace()
ms.SetHbarOmega(20.0); ms.SetTargetMass(8); ms.SetTargetZ(2)
for tz2 in [-1,+1]:
    ms.AddOrbit(0,1,3,tz2,1.0,1)
    ms.AddOrbit(0,1,1,tz2,1.0,1)
ms.SetEmax(1); ms.SetE2max(2); ms.FindEFermi(); ms.SetupKets()
Hs = imsrg.Operator(ms,0,0,0,2)

eom = imsrg.EOM(Hs, 'he8.ref', 0, 0, 0)
rdm = eom.rdm

# Parse he8.ref 3body lines
entries = []
with open('he8.ref') as f:
    for line in f:
        tok = line.split()
        if tok and tok[0]=='TRBTD':
            a,b,c,d,e,f = [int(x)-1 for x in tok[1:7]]
            two_jab = int(tok[7])
            two_jde = int(tok[8])
            two_tot = int(tok[9])
            val     = float(tok[11])
            entries.append((a,b,c,d,e,f, two_jab//2, two_jde//2, two_tot, val))

print(f"he8.ref has {len(entries)} 3-body entries")
ok = 0; fail = 0
for (a,b,c,d,e,f,Jab,Jde,two_tot,ref) in entries:
    srg = rdm.ThreeBody.GetME_pn(Jab, Jde, two_tot, a, b, c, d, e, f)
    diff = abs(srg - ref)
    tol  = 1e-6 * max(1, abs(ref))
    if diff > tol:
        ratio = srg/ref if abs(ref)>1e-12 else float('inf')
        print(f"FAIL: ({a+1},{b+1},{c+1})^{Jab} ({d+1},{e+1},{f+1})^{Jde} 2J={two_tot}"
              f"  ref={ref:.8f}  srg={srg:.8f}  ratio={ratio:.4f}")
        fail += 1
    else:
        ok += 1

print(f"\nPASS: {ok}/{len(entries)}   FAIL: {fail}/{len(entries)}")
