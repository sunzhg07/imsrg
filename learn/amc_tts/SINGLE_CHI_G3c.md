# Single-χ fold for \(\Gamma^{\mathrm{III}_c}\) / G3c

## M-scheme (AMC input)

\[
\Gamma^{\mathrm{III}_c}_{ijkl}
=-\tfrac12\sum_{abcd}f(a,b,c)\,(1-P_{ij})(1-P_{kl})
\Bigl(
\Omega_{abcl}\Omega_{idab}\Gamma_{cjkd}
+\Omega_{icab}\Omega_{abdl}\Gamma_{djkc}
\Bigr)
\]

with \(f(a,b,c)=\bar n_a\bar n_b n_c+n_a n_b\bar n_c\).

## Fold Term2 into Term1 (relabel \(c\leftrightarrow d\))

Term2 after \(c\leftrightarrow d\):

\[
\sum_{abcd}f(a,b,d)\,\Omega_{idab}\Omega_{abcl}\Gamma_{cjkd}.
\]

So (before external \(P\))

\[
\text{Term1}+\text{Term2}
=\sum_{cd}\chi_{id,cl}\,\Gamma_{cjkd},
\qquad
\chi_{id,cl}
=\sum_{ab}\bigl[f(a,b,c)+f(a,b,d)\bigr]\Omega_{idab}\Omega_{abcl}.
\]

- \(f(\cdot,c)\): occupation on **ket-first** of \(\chi\) \(\to\chi_k\)
- \(f(\cdot,d)\): occupation on **bra-second** of \(\chi\) \(\to\chi_j\)

Hence one intermediate

\[
\chi_{ij,kl}=\chi_k(ij,kl)+\chi_j(ij,kl)
\quad\text{(same slots)}.
\]

Index map \(a\leftrightarrow d\), \(b\leftrightarrow c\):

\[
\Gamma^{\mathrm{III}_c}
=-\tfrac12(1-P_{ij})(1-P_{kl})\sum_{ab}\chi_{iabl}\Gamma_{bjka},
\]

same skeleton as scalar Factorized IIe (one \(\chi^\theta\), then Pandya\(\times\bar\Gamma\)).

## What is *not* needed for the fold

- Antisymmetry \(P_{ij}\), \(P_{kl}\) stay outside (already in AMC / Path A).
- Hermiticity of \(\Gamma\) / antihermiticity of \(\Omega\) are **not** required for this rewrite; they matter for relating this \(\chi\) to Factorized’s outer-leg \(\mathrm{Om}\,\mathrm{Occ}_c+(\mathrm{Om}\,\mathrm{Occ}_d)^T\) packaging (different occ placement).

## Failed merges (for contrast)

| Merge | Meaning |
|---|---|
| \(\chi_k(ij,kl)+\chi_j(kl,ij)\) | Factorized-style transpose merge — wrong for AMC fold |
| \(\chi_k+\chi_l\) (occ on both ket legs) | Misread of Term2 after relabel |

## Code

- `ChiThetaToScalarOperator`: phys = \(\chi_k+\chi_j\) same slots → `MakeNotReduced`
- Path B: that Op → Pandya×Γ̄
- `SetUse_TypeGIIIc_single_chi(True)`: Path A Term1 sixj only with \(\chi_k+\chi_j\) (J-scheme check)
- Bench: `run/test_tts_GIIIc_single_chi.py`

## Numerical result (emax=1 He4, random ops)

| Construction | vs TTS |
|---|---|
| Path A dual-strip \(\chi_k\) Term1 + \(\chi_j\) Term2 | **PASS** (~1e-13) |
| Path A Term1 sixj only, \(x=\chi_k+\chi_j\) | **FAIL** (~0.4 ‖Z‖) |
| Path A Term1 only, \(x=\chi_k\) (no Term2) | **FAIL** (~0.6 ‖Z‖) |
| Path B Pandya with \(\chi_k+\chi_j\) same slots | **FAIL** (~1.8× ‖Z‖) |

## Why the m-scheme fold does not give Factorized-style Path A/B

1. **Antisymmetry \((1-P_{ij})(1-P_{kl})\)** is already outside both AMC terms. Dropping Term2 and relying on \(P\) does **not** regenerate it (Term1-only ≠ TTS).

2. **Hermiticity of \(\Gamma\) / antihermiticity of \(\Omega\)** relate matrix elements, but AMC’s Term1 and Term2 use **different sixj ladders** (and different \(\chi\) index routing: \(\chi_k(P,d,c,H)\) vs \(\chi_j(P,c,d,H)\)). Inserting \(\chi_k+\chi_j\) into Term1’s sixjs alone ≠ Term1+Term2.

3. **Scalar Factorized** builds one \(\mathrm{CHI}=\Omega W_c\Omega+(\Omega W_d\Omega)^T\) with *outer-leg* occupations, then Pandya×Γ̄. That is a different packaging; dual-oracle \(\Gamma^{\mathrm{III}_c}\) TTS ≠ Factorized IIe even at \(\lambda=0\). Feeding AMC \(\chi_k+\chi_j\) into that Pandya path does not close the gap.

**Conclusion:** In J-scheme / AMC, \(\Gamma^{\mathrm{III}_c}\) stays a **two-strip** object (Path A). A single m-scheme \(\chi\) exists, but it does not collapse the two AMC sixj topologies into one Factorized Pandya contraction that matches DIRECT.
