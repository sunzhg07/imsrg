# GIVc odd-π T1 — fix plan (revised 2026-09-19)

## Diagnosis (locked)

| Claim | Status |
|---|---|
| Parity / opp-π CC filter broken | **false** |
| Reduced/unreduced / mid \(\hat J\) broken on odd | **false** (even + odd T2 PASS) |
| Old production = **ring** \(\chi_{ibal}\Omega_{ajkb}\) via adcb | **true** |
| m-gold is \(\chi_{ialb}\Omega_{bjak}\) | **true** (`fact_GIVc_which`) |
| Ring ≡ gold iff T1 χ AS (\(n_l=n_b\)) | **true**; fails odd-π T1 only |

## Solution (locked 2026-09-19)

**Gold Path B for all π** — AMC Pandya schemes that land on \((il)|(ab)\):

| Fill | AMC scheme | Reads |
|---|---|---|
| `PandyaChiIalb` | `((1,-3),(2,-4))` on \(\chi_{ialb}\) | \(\chi(i,a,l,b)\) |
| `PandyaOmegaBjak` | `((3,-1),(4,-2))` on \(\Omega_{bjak}\) | \(\Omega(b,j,a,k)\) |

Same ring mid \(w=(-1)^{J_p}/\hat J_p\cdot\hat\lambda^{-1}(-1)^{J_{ab}+\lambda}\)
and corrected inv. Even: gold ≡ adcb. Odd T1: ≡ m/CG.
T1-only → `AddToTBMENonHerm`.

AMC: `learn/amc_tts/factored_GIV/input/G4c_gold_il_ab_pandya.txt`.

## What failed (history)

| Attempt | Why |
|---|---|
| Ordinary JJ pack \((il)|(ab)\) of χ_m | ≢ IMSRG/AMC Pandya |
| Old `PandyaChiIalb` (adcb ph on wrong 9j) | corr ∼0.3–0.5 |
| Blind mid retune | no global fac |
| Path A ninej fold on odd | ≢ m |
| Odd CG leftover (hybrid) | correct but slow — replaced |

## Benches

| Run | Result |
|---|---|
| `PARITY=1 GIVC_CHI=1 TERMS=GIVc … 1 3` | **PASS** |
| `PARITY=0 TERMS=GIVc … 1 2` | **PASS** |
| `PARITY=1 GIVC_CHI=0 TERMS=GIVc … 1 3` | **PASS** |
