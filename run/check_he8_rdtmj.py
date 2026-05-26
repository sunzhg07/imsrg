import sys, os, math
sys.path.insert(0, '../build')
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from pyIMSRG import *

# --- load rdtmj.dat (all orderings, 1-based orbits) ---
rdtmj = []
with open('rdtmj.dat') as f:
    for line in f:
        t = line.split()
        if not t or t[0] != 'TRBTD': continue
        a,b,c   = int(t[1])-1, int(t[2])-1, int(t[3])-1   # 0-based
        d,e,ff  = int(t[4])-1, int(t[5])-1, int(t[6])-1
        jab,jde = int(t[7]),   int(t[8])                   # 2*Jab, 2*Jde
        twoJ    = int(t[9])
        val     = float(t[11])
        rdtmj.append((a,b,c,jab, d,e,ff,jde, twoJ, val))

print(f"rdtmj entries: {len(rdtmj)}")

# --- build EOM from he8.ref ---
ms  = ModelSpace(1, "He8", "p-shell")
ms.SetHbarOmega(20)
ms.SetTargetMass(8)
H   = Operator(ms, 0, 0, 0, 2)
eom = EOM(H, 'he8.ref', 0, 0, 0)
rdm = eom.rdm

# --- categorize results ---
stats = {'can/can': [0,0], 'can/non': [0,0], 'non/can': [0,0], 'non/non': [0,0]}

for (a,b,c,jab, d,e,ff,jde, twoJ, expected) in rdtmj:
    cb = (a<=b and b<=c)
    ck = (d<=e and e<=ff)
    cat = ('can' if cb else 'non') + '/' + ('can' if ck else 'non')
    got = rdm.ThreeBody.GetME_pn(jab//2, jde//2, twoJ, a, b, c, d, e, ff)
    if abs(got - expected) > 1e-5:
        ratio = expected/got if abs(got) > 1e-8 else float('nan')
        if stats[cat][1] < 2:   # show first 2 mismatches per category
            print(f"  {cat}  abc=({a+1},{b+1},{c+1}) Jab={jab//2}"
                  f"  def=({d+1},{e+1},{ff+1}) Jde={jde//2}  2J={twoJ}"
                  f"  got={got:.6f}  exp={expected:.6f}  ratio={ratio:.3f}")
        stats[cat][1] += 1
    else:
        stats[cat][0] += 1

print()
print(f"{'category':10s}  {'match':>6}  {'mismatch':>8}")
for k,(ok,miss) in stats.items():
    print(f"{k:10s}  {ok:6d}  {miss:8d}")
