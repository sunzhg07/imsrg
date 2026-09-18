# GIVc leftover Path B — long debug lesson (2026-09-18)

Production leftover is **only** Pandya → DGEMM → inv (tts_ring).
CG recouple was gold during the lock, then removed.

Bench that locked the mid-J identity:
`run/test_givc_midj_bare_vs_cg.py` (bare \(K\), FoldAS off, T1, λ=3).

---

## What finally matched

m-scheme / CG gold for the bare ring:

\[
K=\sum_{ab}\mathrm{CG}(\lambda\mu,\lambda{-}\mu;00)\,
\chi^\lambda_{ialb}(m)\,\Omega_{bjak}(m)
\]

is **exactly** the J-scheme **tts_ring** Path A

\[
X_{pqsr}
=\sum_{ab}
\chi_{p\,\bar a\,r\,b}\,\Omega_{a\,q\,s\,b}
\qquad(p,q,s,r)=(i,j,k,l)
\]

(WE-reduced χ/Ω; \(X\) reduced; store \(Z_{\mathrm{unred}}=X/\hat J\)).

Path B speed form of **that same ring** (not of printed `G4c_from_chi`):

1. IMSRG **adcb** Pandya of χ and Ω into 2n×2n CC  
   with map \(\mathrm{AMC}\,\bar O(i,j,k,l)\equiv\mathrm{IMSRG\,adcb}(i,l,k,j)\).
2. Mid DGEMM  
   \(\bar G^{J_p}\mathrel{+}=
   \frac{(-1)^{J_p}}{\hat J_p}\,
   \hat\lambda^{-1}(-1)^{J_{ab}+\lambda}\,
   \bar\chi^{J_p J_{ab}}\,\bar\Omega^{J_{ab} J_p}\).
3. Inv Pandya with **no** AMC-sample overall minus  
   \(Z_{\mathrm{red}}=\hat J_0\sum_{J_p}\hat J_p\,
   \{j_r\,j_s\,J_0;\,j_q\,j_p\,J_p\}\,\bar G^{J_p}\).
4. Fermionic \(\tfrac12(1-P)^2\) on the leftover.

Locked: Path A ≡ Path B ≡ CG/m for λ=1,2,3 (T1/T2/both) at emax=1.

---

## Why we did not have it before (traps that burned days)

### 1. Trusted the wrong AMC print as mid-J gold

Printed `G4c_from_chi*_ninej.tex` (χ on `ialb`, Ω on `bjak` with that 6j/9j set)
is **wrong vs m** (~1–3%, non-Hermitian under bra↔ket). NOTES already said so
(2026-07-29), but Path B wiring still chased that print.

**Correct mid-J** was already locked for Ω×Ω as **tts_ring**
(`learn/amc_tts/factored_fIIIa/tts_ring.md`, `run/test_z_ring_mscheme_sign.py`).
GIVc leftover is the same ring with χ in place of the first Ω.

### 2. Same-label AMC Pandya ≠ IMSRG adcb without an explicit map

For **two AS Ω’s**, same-label AMC Pandya happens to work (up to the known
overall minus). For **occupancy-AS χ^λ** (no fermionic AS), same-label AMC
Pandya of labels `(p,b,a,r)` does **not** equal the adcb transform of
`χ(p,b,a,r)`.

Must use:

| AMC mid labels | IMSRG adcb call |
|---|---|
| \(\bar\chi(p,b,a,r)\) | `adcb(p,r,a,b)` ← reads `χ(p,b,a,r)` |
| \(\bar\Omega(a,q,s,b)\) | `adcb(a,b,s,q)` ← reads `Ω(a,q,s,b)` |

Do **not** invent `PandyaChiIalb` / `(ia)(lb)` helpers for production.

### 3. Missing mid factor \((-1)^{J_p}/\hat J_p\)

DGEMM had \(\hat\lambda^{-1}(-1)^{J_{ab}+\lambda}\) but omitted the
channel factor from tts_ring / `test_z_ring`. Ratios looked “almost packaging”
(\(-\sqrt{3}\), \(\hat J^{\pm1}\)) and invited blind hat retunes — those are
usually packaging, not angular bugs; here the factor was a real missing piece
of the ring Path B formula.

### 4. Compared Term1 alone to full AS leftover

`½(1-P)^2` has four sectors. Locking **bare** \(K\) / `FoldAS=false` first
avoids that class of fake FAIL.

### 5. Reduced vs unreduced / √2 identical (again)

`GetTBME_J` = unreduced tilde (identical √2 restored). AMC `reduce=true` =
\(S/\hat J\). Dividing √2 again when comparing to `GetTBME_J` fakes 0.5 / 0.707.
See `learn/amc_tts/REDUCED_UNREDUCED.md`.

### 6. 2n pack: unnormalized vs normalized

`GetTBME_J` / ChiTab = **unnormalized**. `GetMatrix` = normalized.
Pandya 2n packs unnormalized; ÷√2 only on `AddToTBME`.

---

## Playbook if a similar T×T→S leftover fails again

1. Name packaging on both sides (reduced / unreduced / m).
2. Strip \(P_{ij}P_{kl}\); lock **one bare ME** mid-J vs m/CG.
3. Prefer a **locked ring cousin** (`tts_ring`, `omega_cross`) over a new AMC
   fold print until that print passes m.
4. For Path B: lock fwd Pandya alone, then mid product, then inv — separately.
5. Use IMSRG adcb + explicit AMC↔adcb map; do not mix same-label ninejs into
   adcb slots.
6. Do not retune \(\hat J\)/√2/phases until (1)–(5) are clean.

---

## Code

- ethS: `comm223_232_GIVc` → Pandya+DGEMM+inv only
  (`src/FactorizedDoubleCommutator_eths.cc`).
- Packaging table: this directory’s `EQUATION.md`.
