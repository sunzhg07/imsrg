# Ω₁_ibal × Ω₂_ajkb (no occ) — Path A vs Pandya

## Setup

\[
\chi_{ijkl}=\sum_{ab}\Omega^{(1)}_{ibal}\,\Omega^{(2)}_{ajkb}
\]

| Case | Ω₁ | Ω₂ | χ |
|---|---|---|---|
| **SS** | scalar, `reduce=true` | scalar, `reduce=true` | scalar, `reduce=true` |
| **ST** | scalar, `reduce=true` | tensor | **tensor** (rank λ₂) |
| **TS** | tensor | scalar, `reduce=true` | **tensor** (rank λ₁) |

## Ground truth: Neithan Path B

Parzuchowski / `neithan.tex` ST ph core (no occ, drop \((1-P)\)):

\[
\chi^{(J_0 J_1)\lambda}_{ijkl}
=\sum_{J_3 J_4}\hat J_0\hat J_1\hat J_3\hat J_4
(-1)^{j_j+j_l+J_1+J_4}
\begin{Bmatrix}j_i&j_l&J_3\\ j_j&j_k&J_4\\ J_0&J_1&\lambda\end{Bmatrix}
\sum_{ab}
\bar S^{J_3}_{i\bar l a\bar b}\,
\bar T^{(J_3 J_4)\lambda}_{a\bar b k\bar j}
\]

with Neithan Pandya (scalar uses \(\hat J'\), **not** \(\hat J'^2\)):

\[
\bar S^{J}_{p\bar q r\bar s}=-\sum_{J'}\hat J'\,
\begin{Bmatrix}j_p&j_q&J\\ j_r&j_s&J'\end{Bmatrix}S_{psrq}
\]

\[
\bar T^{(J_1 J_2)\lambda}_{p\bar q r\bar s}
=-\sum\hat J_1\hat J_2\hat J_3\hat J_4
(-1)^{j_q+j_s+J_2+J_4}
\begin{Bmatrix}j_p&j_s&J_3\\ j_q&j_r&J_4\\ J_1&J_2&\lambda\end{Bmatrix}
T^{(J_3 J_4)\lambda}_{psrq}.
\]

**AMC reduced Path B ≡ Neithan** (machine precision): reduced scalar
\(\bar S=-\hat J\sum\hat J'\,6j\,S\) plus product \(\hat J^{-1}\bar S\bar T\) cancels the extra \(\hat J\), recovering Neithan’s bare \(\hat J'\) packaging. Tensor Pandya already matches.

## What was wrong

| Piece | Status |
|---|---|
| Path B (AMC Pandya product + inv) | **correct** for SS, ST, **and TS** |
| AMC printed **ST** direct | correct (matches Path B) |
| AMC printed **TS** direct | **wrong** (rel ~1.7 vs Path B / Neithan) |

TS Path A in `run/test_o1o2_ph_cases.py` is now the **analytic expansion** of Neithan Path B (substitute Pandya defs), not AMC’s `o1o2_ts_direct_*.tex`.

Expanded TS phase keeps \(2j_l\) (for fermions \((-1)^{2j_l}=-1\)); dropping it flips the global sign.

## Numeric proof (`run/test_o1o2_ph_cases.py`, emax=1)

| Case | Result |
|---|---|
| SS (λ=0,0) | **PASS** (~1e-15) |
| ST (λ=0,2) | **PASS** (~1e-15) |
| TS (λ=2,0) | **PASS** (~1e-15) after replacing AMC TS direct |


Cross-check: `run/test_neithan_ph_match.py` — Neithan B ≡ expanded D; AMC_B ≡ Neithan_B(Ĵ).

```bash
PYTHONPATH=... python3 run/test_o1o2_ph_cases.py
PYTHONPATH=... python3 run/test_neithan_ph_match.py
```
