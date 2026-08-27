# Γ^{IV_b} Factorized flow: Pandya Ω × RC χ → Z

**Schemes (this project):**

| Name | Coupling | Example |
|---|---|---|
| **Pandya** | `((1,-4),(3,-2))` | \(\bar O_{a\bar b\,c\bar d}\leftarrow O_{adcb}\) |
| **RC** (here) | leg-recouple of Pandya χ for DGEMM | \(\overline{\overline\chi}_{a\bar b\,c\bar d}\leftarrow\bar\chi_{a\bar d\,b\bar c}\) (+ pack) |

Arxiv has **two** double-bar wirings (CC #1 = \(j\bar l\,k\bar i\), CC #2 = \(i\bar k\,l\bar j\)).  
IV_b contracts in the **CC #1 / Pandya \(j\bar l\,k\bar i\)** slot.  
Factorized keeps **Pandya \(\Omega\) only** (never RC \(\Omega\)) and builds RC on \(\chi^\iota\).

**Gold:** m-scheme / AMC bare χ. **No TTS.**  
**χ^ι** is neither Hermitian nor AS → Pandya/RC storage is **2n×2n**; never fill with \(h_\chi\).  
Only \(\Omega\) may be swapped / transposed (with reduced phase).

Refs: `FactorizedDoubleCommutator.cc` L1289–1350, L1833–2013;  
`learn/factorized_code_analyze.tex` §code-GIVb; `NOTES.md` §Γ^{IV_b}.

---

## M-scheme target

\[
Z_{ijkl}
=(1-P_{ij})(1-P_{kl})
\sum_{ab}
\bigl(
\chi^\iota_{aibk}\,\Omega_{jbla}
-\chi^\iota_{akbi}\,\Omega_{jalb}
\bigr).
\]

Trusted χ (η-analog; spectator \(j\)):

\[
\chi^\iota_{ijkl}
=\sum_{ab}
\bigl(\bar n_a n_b\bar n_k+n_a\bar n_b n_k\bigr)
\,\Omega_{bjka}\,\Gamma_{iabl}.
\]

---

## Scalar Factorized flow (λ = 0) — **LOCKED**

```
Γ,Ω ──Pandya──► Γ̄, Ω̄
                  │
                  ▼
            χ̄ = Γ̄·(occ⊙Ω̄)          [Pandya]
                  │
                  ▼
            RC[χ̄]  (leg + pack)      [not Pandya]
                  │
                  ▼
            W = Ω̄ · RC               [Pandya slot j̄l k̄i]
                  │
                  ▼
            Inv_(1−P)²(W) → Z_ijkl
```

Bench: `run/test_G4b_factorized_fullZ.py` — **PASS**.

### 1. Pandya `((1,-4),(3,-2))`

\[
\bar O^{J}_{a\bar b\,c\bar d}
=-\sum_{J'}(2J'+1)
\begin{Bmatrix} j_a & j_b & J \\ j_c & j_d & J' \end{Bmatrix}
O^{J'}_{a\,d\,c\,b}.
\]

\(O=\Omega\) or \(\Gamma\). Hermiticity / AS fills **only** for \(\Omega\) (\(h_\Omega\)).

### 2. Pandya \(\bar\chi^\iota\)

\[
\bar\chi^{\iota\,J}_{i\bar l\,k\bar j}
=\sum_{ab}
\mathrm{occ}_{A\bar BC}(a,b,k)\,
\bar\Gamma^{J}_{i\bar l\,a\bar b}\,
\bar\Omega^{J}_{a\bar b\,k\bar j}.
\]

Code: `bar_CHI_V = bar_Gamma * nnnbar_Eta`. No \(h_\Omega\) phase, no transpose.

### 3. Pandya χ → RC χ (code L1833–1836)

OUT: bra \((a,b)\), ket \((c,d)\), channel \(J\).

\[
\overline{\overline\chi}{}^{\iota\,J}_{a\bar b\,c\bar d}
=
\sum_{J'}(2J'+1)\,(-1)^{j_b+j_c+J'}
\begin{Bmatrix} j_a & j_b & J \\ j_c & j_d & J' \end{Bmatrix}
\Bigl(
  \bar\chi^{\iota\,J'}_{a\bar d\,b\bar c}
  -h_Z\,\bar\chi^{\iota\,J'}_{b\bar c\,a\bar d}
\Bigr).
\]

| Pack term | Fold topology |
|---|---|
| \(\bar\chi_{a\bar d\,b\bar c}\) | \(\chi_{aibk}\Omega_{jbla}\) |
| \(-h_Z\,\bar\chi_{b\bar c\,a\bar d}\) | \(\chi_{akbi}\Omega_{jalb}\) |

### 4. DGEMM — Pandya Ω × RC χ

\[
W^{J}_{j\bar l\,k\bar i}
=
\sum_{ab}
\bar\Omega^{J}_{j\bar l\,a\bar b}\,
\overline{\overline\chi}{}^{\iota\,J}_{a\bar b\,k\bar i}.
\]

Code: `CHI_V_final = bar_Eta * bar_CHI_V_RC`.  
\(W\) sits in Pandya / CC #1 slot \(j\bar l\,k\bar i\).

### 5. Inv → \(Z_{ijkl}\) (L1884–2013)

Read \(W\) at four exchanges; assemble \((1-P_{ij})(1-P_{kl})\). Schematic:

\[
\begin{aligned}
C_{jikl}
&=
-\sum_{J'}(2J'+1)\,(-1)^{J'+j_i+j_k}
\begin{Bmatrix} j_j & j_i & J \\ j_k & j_l & J' \end{Bmatrix}
W^{J'}_{j\bar l\,i\bar k},
\\[0.5em]
C_{ijlk}
&=
-\sum_{J'}(2J'+1)\,(-1)^{J'+j_j+j_l}
\begin{Bmatrix} j_j & j_i & J \\ j_k & j_l & J' \end{Bmatrix}
W^{J'}_{i\bar k\,j\bar l},
\end{aligned}
\]

and the \((jk,il)\) / \((il,jk)\) pair with
\(\mathrm{SixJ}(j_i,j_j,J;\,j_k,j_l,J')\). Then

\[
\begin{aligned}
Z^{J}_{ijkl}
&=
(-1)^{J+j_i+j_j}
\Big[
  C_{jikl}
  -(-1)^{(j_i+j_j)/2-J}\,C_{ijkl}
  -(-1)^{(j_l+j_k)/2-J}\,C_{jilk}
  +(-1)^{(j_k+j_l+j_i+j_j)/2}\,C_{ijlk}
\Big]
\end{aligned}
\]

(\(\sqrt2\) store norms when \(i=j\) or \(k=l\).)

**\((1-P)^2\) only here** — not inside χ, RC, or GEMM.

---

## Tensor continuous (λ ≠ 0) — candidate

Same flow; rectangular mid-\(J\). \(\lambda\to 0\) must recover the scalar formulas above.  
**Status:** χ̄ Path B locked all λ; full RC×Ω̄×Inv twin **not** locked yet (fold path is gold for λ≠0).  
NineJ layouts below are the **code continuous** (same pack/rewire as scalar); validate \(\lambda\to0\) numerically before trusting as production.

### 1′. Tensor Pandya Ω (Path B locked)

\[
\bar\Omega^{J_0 J_1\lambda}_{a\bar b\,c\bar d}
=
-\sum_{J_1' J_2'}
\hat J_1'\hat J_2'\hat J_0\hat J_1
\,(-1)^{j_b+j_d+J_1+J_2'}
\begin{Bmatrix}
j_a & j_d & J_1' \\
j_b & j_c & J_2' \\
J_0 & J_1 & \lambda
\end{Bmatrix}
\Omega^{J_1' J_2'\lambda}_{a\,d\,c\,b}.
\]

(IMSRG / Path B; **no** extra \(\sqrt{(2J_1'+1)/(2J_0+1)}\) beyond the hats in the print you lock against.)  
Scalar Γ Pandya unchanged (equal-\(J\)).

### 2′. Rectangular \(\bar\chi^\iota\) (Path B locked)

\[
\bar\chi^{\iota\,J_0 J_1\lambda}_{i\bar l\,k\bar j}
=\sum_{ab}
\mathrm{occ}_{A\bar BC}(a,b,k)\,
\bar\Gamma^{J_0}_{i\bar l\,a\bar b}\,
\bar\Omega^{J_0 J_1\lambda}_{a\bar b\,k\bar j}.
\]

No \(h_\Omega\), no \((-1)^{J_0+J_1}\), no transpose (unlike χ^κ VI_II).

### 3′. RC continuous (Path 2 — code pack + rewire)

IMSRG `NineJ(ja,jb,J, jc,jd,J, Jp,Jp,0)` reduces to
\((-1)^{j_b+j_c+J+J'}/(\hat J\hat J')\,\mathrm{SixJ}(ja,jb,J;\,jd,jc,J')\) —
**not** Factorized’s \(\mathrm{SixJ}(ja,jb,J;\,jc,jd,J')\).
To recover scalar §3 under `NineJ`, use the **middle-row swap** \((j_d,j_c)\)
and phase \(\Phi=(-1)^{j_b+j_d+J_0+J_3}\):

\[
\overline{\overline\chi}{}^{\iota\,J_0 J_1\lambda}_{a\bar b\,c\bar d}
=
\sum_{J_3 J_4}
(-1)^{j_b+j_c+J_4}\,
\Phi\,
\hat J_0\hat J_1\hat J_3\hat J_4
\cdot\frac{\hat J_3}{\hat J_0}
\cdot
\begin{Bmatrix}
j_a & j_b & J_0 \\
j_d & j_c & J_1 \\
J_3 & J_4 & \lambda
\end{Bmatrix}
\Bigl(
  \bar\chi^{\iota\,J_3 J_4\lambda}_{a\bar d\,b\bar c}
  -h_Z\,\bar\chi^{\iota\,J_4 J_3\lambda}_{b\bar c\,a\bar d}
\Bigr),
\]

with \(\Phi=(-1)^{j_b+j_d+J_0+J_3}\) and NineJ =
`NineJ(ja,jb,J0, jd,jc,J1, J3,J4,λ)`.  
Transpose pack uses **swapped** mid ranks \((J_4,J_3)\) on \(\bar\chi_{bc,ad}\).

**λ→0 check** (\(J_0=J_1=J\), \(J_3=J_4=J'\)): recovers \((2J'+1)\,6j\) of scalar §3
(numerically exact on IMSRG `NineJ`).  
**λ≠0:** still a candidate — full \(Z\) vs m **FAIL** in
`run/test_G4b_factorized_tensor.py` (fold remains gold).

### 4′. Rectangular DGEMM

\[
W^{J_b J_k\lambda}_{j\bar l\,k\bar i}
=
\sum_{ab\,J_m}
\bar\Omega^{J_b J_m\lambda}_{j\bar l\,a\bar b}\,
\overline{\overline\chi}{}^{\iota\,J_m J_k\lambda}_{a\bar b\,k\bar i}.
\]

Same Pandya / CC #1 free-label slot; mid \(J\) summed.

### 5′. Inv continuous → scalar \(Z^{J}_{ijkl}\)

\(Z\) is scalar (rank 0). Same middle-row swap + \(\Phi\) as §3′ so λ→0
matches Factorized SixJ. Sector 1
(\(\mathrm{SixJ}(j_j,j_i,J;\,j_k,j_l,J')\rightarrow W_{j\bar l\,i\bar k}\)):

\[
\begin{aligned}
C_{jikl}
&=
-\sum_{J_3 J_4}
(-1)^{J_3+j_i+j_k}\,
\Phi_1\,
\hat J\,\hat J\,\hat J_3\hat J_4
\cdot\frac{\hat J_3}{\hat J}
\cdot
\begin{Bmatrix}
j_j & j_i & J \\
j_l & j_k & J \\
J_3 & J_4 & \lambda
\end{Bmatrix}
W^{J_3 J_4\lambda}_{j\bar l\,i\bar k},
\end{aligned}
\]

with \(\Phi_1=(-1)^{j_i+j_l+J+J_3}\) and NineJ =
`NineJ(jj,ji,J, jl,jk,J, J3,J4,λ)` (\(J_3\) on bra \((j,l)\), \(J_4\) on ket \((i,k)\)).  
Partner read \(W_{i\bar k\,j\bar l}\) at \((J_4,J_3)\).  
Sector 2: `NineJ(ji,jj,J, jl,jk,J, J3,J4,λ)` with
\(\Phi_2=(-1)^{j_j+j_l+J+J_3}\) and reads \(W_{j\bar k\,i\bar l}\) / \(W_{i\bar l\,j\bar k}\).  
Assemble \(Z\) with the same \((1-P)^2\) phase pattern as scalar §5.

**λ→0:** recovers \((2J'+1)\,6j\). **λ≠0:** same open status as §3′.

---

## What not to confuse

| Object | Role |
|---|---|
| Pandya `((1,-4),(3,-2))` | \(\bar\Omega\), \(\bar\Gamma\), \(\bar\chi\), and final \(W\) storage slot |
| RC (Factorized) | leg-recouple + ι pack on χ only |
| Arxiv CC #1 vs #2 | IV_b uses #1 slot; do **not** pick #2 |
| AMC `03a` / `03b` | angular skeletons; **not** drop-in for Factorized RC |
| Analyze fold | dual of RC×Ω̄ (gold any λ); no RC |
| \((1-P)^2\) | **only** Inv → \(Z\) |

---

## Status

| Piece | λ=0 | λ≠0 |
|---|---|---|
| χ̄ Path B ≡ AMC ≡ m | **PASS** | **PASS** |
| RC §3 / §3′ | **PASS** (full \(Z\)) | λ→0 NineJ fixed; rectangular **FAIL** vs m |
| Ω̄·RC → Inv → \(Z\) | **PASS** `test_G4b_factorized_fullZ.py` / `test_G4b_factorized_tensor.py` | **FAIL** same; use fold gold |
| Fold \(W_1-W_2\to(1-P)^2\to Z\) | **PASS** | **PASS** λ=0,1,2 — **WRAP-UP** m≡AMC≡Path B |

Related: [`RC_DUAL.md`](RC_DUAL.md), [`EQUATION_CHAIN.md`](EQUATION_CHAIN.md) (AMC alternate), [`../NOTES.md`](../NOTES.md).
