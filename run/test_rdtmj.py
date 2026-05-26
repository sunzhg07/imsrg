import sys, os
sys.path.insert(0, '/Users/wolf/work/imsrg/build')
os.chdir('/Users/wolf/work/imsrg/run')
from pyIMSRG import *

# Step 1: read he8.ref, write it back out
ms  = ModelSpace(1, "He4", "He4")
Hs  = Operator(ms, 0, 0, 0, 2)
eom = EOM(Hs, 'he8.ref', 0, 0, 0)
eom.WriteTdm(eom.rdm, 'he8_out.ref')
print("Written he8_out.ref")

# Step 2: compare TRBTD lines between he8.ref and he8_out.ref
# key = (a,b,c,d,e,f,2Jab,2Jde,2Jtot) — ignore 2*Tz column if present
def parse_3b(fname):
    elems = {}
    with open(fname) as f:
        for line in f:
            tok = line.split()
            if not tok or tok[0] != 'TRBTD':
                continue
            a,b,c,d,e,f = tok[1],tok[2],tok[3],tok[4],tok[5],tok[6]
            jab,jde,jtot = tok[7],tok[8],tok[9]
            key = (a,b,c,d,e,f,jab,jde,jtot)
            val = float(tok[-1])   # last column always the value
            elems[key] = val
    return elems

orig = parse_3b('he8.ref')
out  = parse_3b('he8_out.ref')

print(f"he8.ref    TRBTD count: {len(orig)}")
print(f"he8_out.ref TRBTD count: {len(out)}")

errors = 0
for key, val_out in out.items():
    if key in orig:
        diff = abs(val_out - orig[key])
        if diff > 1e-6:
            print(f"VALUE MISMATCH: {key}  orig={orig[key]:+.8f}  out={val_out:+.8f}  diff={diff:.2e}")
            errors += 1
    else:
        print(f"KEY NOT IN ORIG: {key}  val={val_out:+.8f}")
        errors += 1

if errors == 0:
    print("PASS: all written elements match he8.ref")
