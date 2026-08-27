# RC χ^ι — two paths that must match

χ^ι is **neither Hermitian nor AS** → Pandya/RC storage is always **2n×2n**
(direct + exchange). Never fill with \(h_\chi\); pack is topology only.

Goal: lock
\[
\mathrm{RC}[\chi]_{\text{from normal}}
\;=\;
\mathrm{RC}[\bar\chi]_{\text{from Pandya}}
\]
(same labels / same \(2n\) layout), then use either in ×Ω̄ → Inv.

---

## Path 1 — normal χ → RC (AMC `03b`, pack by hand)

**Input:** ordinary \(\chi^\iota\) (Path B invPlus / AMC direct / m).  
**AMC:** same-label scheme change only → `((1,-3),(4,-2))`.

### Scalar (`03b_normal_to_CC_scalar.tex`)

\[
\overline{\overline\chi}{}^{\iota\,J}_{ijkl}
=
(-1)^{j_k+j_l}
\sum_{J'}(-1)^{J'}\,\hat J'^2
\begin{Bmatrix} j_j & j_l & J \\ j_k & j_i & J' \end{Bmatrix}
\chi^{\iota\,J'}_{ijkl}.
\]

### Tensor (`03b_normal_to_CC_tensor.tex`)

\[
\overline{\overline\chi}{}^{\iota\,J_0 J_1\lambda}_{ijkl}
=
(-1)^{J_0+j_i+j_k+\lambda}\,
\hat J_0\hat J_1
\sum_{J_2 J_3 j_0}
(-1)^{J_2+J_3}\,
\hat J_2\hat J_3\,\hat j_0^2
\begin{Bmatrix}\lambda & J_3 & J_2 \\ j_i & j_j & j_0\end{Bmatrix}
\begin{Bmatrix}j_i & j_k & J_0 \\ j_l & j_0 & J_3\end{Bmatrix}
\begin{Bmatrix}J_1 & J_0 & \lambda \\ j_0 & j_j & j_l\end{Bmatrix}
\chi^{\iota\,J_2 J_3\lambda}_{ijkl}.
\]

### Pack (hand — not in AMC)

\[
\mathrm{RC}_{ab,ik}
=
\overline{\overline\chi}{}_{abik}
-h_Z\,\overline{\overline\chi}{}_{ikab}.
\]

This is the **official angular path** in `EQUATION_CHAIN.md`. Status: AMC OK;
full \(Z\) via 03b→04 still WIP (wiring vs Factorized dual).

---

## Path 2 — Pandya χ̄ → RC (Factorized **code**, not raw AMC `03a`)

**Input:** \(\bar\chi^\iota\) in Pandya (`bar_CHI_V`, Path B locked).  
**Gold identity:** `FactorizedDoubleCommutator.cc` L1833–1841.

### Scalar / λ=0 (code — locked in full \(Z\))

\[
\mathrm{RC}[\bar\chi^\iota]^{J}_{ab,cd}
=
\sum_{J'}(2J'+1)\,(-1)^{j_b+j_c+J'}
\begin{Bmatrix} j_a & j_b & J \\ j_c & j_d & J' \end{Bmatrix}
\Bigl(
  \bar\chi^{\iota\,J'}_{ad,bc}
  -h_Z\,\bar\chi^{\iota\,J'}_{bc,ad}
\Bigr).
\]

Notes:
- Index **rewire** \(ab,cd\leftarrow ad,bc\) is part of the definition (IIb leg recouple).
- Pack is inside the sum (ι: relative \(-h_Z\); η uses \(+\) and overall \(-\)).
- Needs **2n** reads of \(\bar\chi\) (ad/bc and bc/ad, including exchange offsets).

### AMC `03a` skeleton (angular only — **not** drop-in for Path 2)

Same-label Pandya→RC, pack stripped:

| | formula |
|---|---|
| scalar | \(\overline{\overline\chi}{}^J_{ijkl}=(-1)^{J+j_k+j_l}\sum(-1)^{J'}\hat J'^2\{j_l j_j J;\,j_k j_i J'\}\bar\chi^{J'}_{ijkl}\) |
| tensor | ninej continuous (`03a_RC_scheme_change_tensor.tex`) |

This is **Oc=Op twin** of scheme change only. Bench (`test_chi_iota_rc.py`):
AMC⊗pack ≠ code RC on many samples. **Do not implement Path 2 from raw `03a`.**

### Tensor Path 2 — equation we actually need

Continuous of the **code** (pack + rewire), with hats so \(\lambda\to 0\) recovers
\((2J'+1)\times 6j\):

\[
\mathrm{RC}[\bar\chi^\iota]^{J_0 J_1\lambda}_{ab,cd}
=
\sum_{J_3 J_4}
(-1)^{j_b+j_c+J_4}\,
\hat J_0\hat J_1\hat J_3\hat J_4
\cdot
\frac{\hat J_3}{\hat J_0}
\cdot
\begin{Bmatrix}
j_a & j_b & J_0 \\
j_c & j_d & J_1 \\
J_3 & J_4 & \lambda
\end{Bmatrix}
\Bigl(
  \bar\chi^{\iota\,J_3 J_4\lambda}_{ad,bc}
  -h_Z\,\bar\chi^{\iota\,J_3 J_4\lambda}_{bc,ad}
\Bigr).
\]

Identity check at \(\lambda=0\), \(J_0=J_1=J\), \(J_3=J_4=J'\):

\[
\hat J^2\hat J'^2\cdot\frac{\hat J'}{\hat J}
\cdot
\frac{\{j_a j_b J;\,j_c j_d J'\}}{\hat J\,\hat J'}
=
\hat J'^2\,\{j_a j_b J;\,j_c j_d J'\}
=(2J'+1)\,6j,
\]

matches Factorized. (NineJ layout above is the standard
`NineJ(ja,jb,J0, jc,jd,J1, J3,J4,λ)` continuous of `SixJ(ja,jb,J,jc,jd,J')`.)

**Status:** scalar Path 2 locked (full \(Z\) @λ=0). Tensor continuous above is the
candidate to lock vs Path 1 (and vs λ=0 code). Prior full-\(Z\) attempts failed when
exchange-sector / fill rules were wrong — fix storage first, then compare RC MEs.

---

## How they should match

```
χ_normal  ──03b+pack──►  RC₁
    ▲                      ║ must ≡
    │ invPlus              ║
χ̄_Pandya ──code RC──────►  RC₂
```

| Compare | What |
|---|---|
| λ=0 | RC₂ ≡ Factorized (known). RC₁ ≡ RC₂ is the open dual lock |
| λ≠0 | Build RC₂ from tensor code continuous; RC₁ from 03b tensor + pack; compare in 2n CC |

Pack is the same topology on both sides: \((ad,bc)-(bc,ad)\) on Pandya side ↔
\((abik)-(ikab)\) after normal→CC (labels differ by scheme).

**Production Factorized flow (Pandya Ω × RC χ → Z), scalar + tensor continuous:**
[`FACTORIZED_RC_FLOW.md`](FACTORIZED_RC_FLOW.md).

---

## What not to confuse

| Object | Role |
|---|---|
| AMC `03b` | Path 1 angular (normal→CC) — **keep** |
| Factorized L1833 | Path 2 gold (Pandya→RC + pack) — **keep** |
| AMC `03a` | Angular skeleton only; **not** Path 2 implementation |
| Analyze fold | Bypass RC entirely (locked full \(Z\) any λ) |
| \((1-P)^2\) | Only on final InvPandya / fold to \(Z\), never inside RC |

---

## Bench status (`run/test_chi_iota_rc_dual.py`)

```
PYTHONPATH=build python3 -B run/test_chi_iota_rc_dual.py [emax] [lambda] [nsamp]
```

| Check | Result (emax=1, λ=0) |
|---|---|
| Path 2 Factorized RC builds | OK (nonzero samples) |
| Path 1 `pack_abik` ≡ Path 2 | **FAIL** — mostly 0 (wrong triangles on Pandya slot) |
| Path 1 `dbar(a,d,b,c)` pack ≡ Path 2 | **FAIL** — best ~6/25 at ±1 (`neg_adbc_nop`); ratios mixed |

**Channel note:** Path 2 lives in Pandya triangles \((j_a j_b),(j_c j_d)\).
03b scheme `((1,-3),(4,-2))` has triangles \((j_i j_k),(j_l j_j)\). Matching
the Path 2 slot requires free labels \((i,j,k,l)=(a,d,b,c)\), i.e.
\(\overline{\overline\chi}_{adbc}-h_Z\overline{\overline\chi}_{bcad}\), not
raw `pack_abik` on \((a,b,c,d)\). Even with that remapping, ME dual is open.

Locked gold (not dual RC):
- `test_chi_iota_m_vs_amc.py` PASS
- `test_chi_iota_pathB_vs_direct.py` PASS
- `test_G4b_factorized_fullZ.py` PASS @λ=0
- `test_G4b_pathB_fold_mscheme.py` PASS (any λ via fold)

## Next lock steps

1. λ=0: close Path 1 ≡ Path 2 after channel remapping (phase / 6j orientation / hats).
2. If ME dual stays open: compare after ×Ω̄ (same \(W\)) rather than raw RC MEs.
3. λ≠0: implement tensor Path 2 candidate; compare to Path 1 tensor.

Refs: `output/steps/EQUATION_CHAIN.md`, `03_pathB/RC_EQUATION.md`,
`FactorizedDoubleCommutator.cc` L1833–1841,
`run/test_chi_iota_rc.py`, `run/test_G4b_factorized_fullZ.py`.
