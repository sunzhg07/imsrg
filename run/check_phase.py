import sys, math, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/Users/wolf/work/imsrg/build')
import pyIMSRG as imsrg

# ── 1. Read he8.ref via EOM ─────────────────────────────────────────────────
with open('he8.ref') as f:
    lines = f.readlines()
norb = int(lines[1].strip())
ms0 = imsrg.ModelSpace()
for i in range(norb):
    tok = lines[2+i].split()
    ms0.AddOrbit(int(tok[1]), int(tok[2]), int(tok[3]), int(tok[4]), 0.0, 1)
Hs = imsrg.Operator(ms0, 0, 0, 0, 2)
eom = imsrg.EOM(Hs, 'he8.ref', 0, 0, 0)
rdm = eom.rdm
rdm_ms = rdm.GetModelSpace()
print(f"norb={rdm_ms.GetNumberOrbits()}  3b_channels={rdm_ms.GetNumberThreeBodyChannels()}  allocated={rdm.ThreeBody.IsAllocated()}")

# ── 2. Write rdm back out ────────────────────────────────────────────────────
eom.WriteTdm(rdm, '/tmp/he8_out.ref')

# ── 3. Read both files and compare 3b elements ──────────────────────────────
def read_3b(fname):
    """Returns dict (a,b,c,d,e,f,2Jab,2Jde,2J) -> value (1-based orbit indices)"""
    d = {}
    for line in open(fname):
        tok = line.split()
        if not tok or tok[0] != 'TRBTD':
            continue
        a,b,c,dd,e,f = tok[1],tok[2],tok[3],tok[4],tok[5],tok[6]
        two_jab,two_jde,two_tot = tok[7],tok[8],tok[9]
        # skip tok[10] if it is 2*Tz (he8.ref has it, WriteTdm does not)
        # value is always the last token
        val = float(tok[-1])
        key = (a,b,c,dd,e,f,two_jab,two_jde,two_tot)
        d[key] = val
    return d

orig = read_3b('he8.ref')
out  = read_3b('/tmp/he8_out.ref')

print(f"\nhe8.ref  3b elements : {len(orig)}")
print(f"he8_out  3b elements : {len(out)}")

nerr = 0
for key, val_orig in orig.items():
    val_out = out.get(key)
    if val_out is None:
        print(f"MISSING in output: {key}")
        nerr += 1
    elif abs(val_orig - val_out) > 1e-6 * max(1, abs(val_orig)):
        print(f"MISMATCH {key}: orig={val_orig:.8f}  out={val_out:.8f}")
        nerr += 1
if nerr == 0:
    print("ALL MATCH: read→write roundtrip is exact")

# ── 4. Verify GetME_pn recoupling: (ab)^Jab c  vs  (ba)^Jab c ───────────────
print("\n── Phase check: (ab)^J vs (ba)^J ──")
ncheck = nerr2 = 0
for key, val in orig.items():
    a,b,c,d,e,f,two_jab,two_jde,two_tot = [int(x) for x in key]
    Jab = two_jab // 2
    Jde = two_jde // 2
    # 0-based
    ai,bi,ci,di,ei,fi = a-1,b-1,c-1,d-1,e-1,f-1
    v_ab  = rdm.ThreeBody.GetME_pn(Jab, Jde, two_tot, ai, bi, ci, di, ei, fi)
    v_ba  = rdm.ThreeBody.GetME_pn(Jab, Jde, two_tot, bi, ai, ci, di, ei, fi)
    # expected phase for bra swap: (-1)^(ja/2 + jb/2 - Jab + 1) in imsrg convention
    ja2 = rdm_ms.GetOrbit(ai).j2
    jb2 = rdm_ms.GetOrbit(bi).j2
    phase = 1 if ((ja2 + jb2)//2 - Jab + 1) % 2 == 0 else -1
    expected_ba = phase * v_ab
    diff = abs(v_ba - expected_ba)
    ncheck += 1
    if diff > 1e-8:
        print(f"FAIL (ab)↔(ba): orbs=({a},{b},{c}|{d},{e},{f}) Jab={Jab} Jde={Jde} J={two_tot}  v_ab={v_ab:.6f}  v_ba={v_ba:.6f}  expected={expected_ba:.6f}")
        nerr2 += 1
if nerr2 == 0:
    print(f"ALL {ncheck} phase checks PASS")

# ── 5. Map he8.ref elements to rdtmj.dat (all orderings) ────────────────────
print("\n── he8.ref vs rdtmj.dat (GetME_pn recoupling) ──")
rdtmj = []
for line in open('rdtmj.dat'):
    tok = line.split()
    if not tok or tok[0] != 'TRBTD': continue
    a,b,c,d,e,f = [int(x)-1 for x in tok[1:7]]
    two_jab,two_jde,two_tot = int(tok[7]),int(tok[8]),int(tok[9])
    val = float(tok[-1])
    rdtmj.append((a,b,c,d,e,f,two_jab,two_jde,two_tot,val))

print(f"rdtmj.dat entries: {len(rdtmj)}")
nerr3 = 0
for (a,b,c,d,e,f,two_jab,two_jde,two_tot,val_ref) in rdtmj:
    Jab = two_jab // 2
    Jde = two_jde // 2
    srg = rdm.ThreeBody.GetME_pn(Jab, Jde, two_tot, a, b, c, d, e, f)
    diff = abs(srg - val_ref)
    if diff > 1e-5 * max(1e-10, abs(val_ref)):
        ratio = srg/val_ref if abs(val_ref) > 1e-12 else float('inf')
        print(f"MISMATCH ({a+1},{b+1},{c+1}|{d+1},{e+1},{f+1}) Jab={Jab} Jde={Jde} 2J={two_tot}: srg={srg:.8f}  ref={val_ref:.8f}  ratio={ratio:.4f}")
        nerr3 += 1
if nerr3 == 0:
    print(f"ALL {len(rdtmj)} entries MATCH")
else:
    print(f"{nerr3}/{len(rdtmj)} MISMATCHES")
