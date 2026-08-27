# Γ^I / χ^ε and Γ^II / χ^ζ — AMC + verification

**Chain:** m-scheme → AMC DIRECT → AMC factorized (+ DGEMM) → ethS  
**Code:** `comm223_232_chi1b_tensor` in `FactorizedDoubleCommutator_eths.cc`

---

## Γ^I / χ^ε

### AMC inputs (`factored_GI/input/`)

| File | Role |
|---|---|
| `chi_epsilon.txt` | χ^ε alone (scalar 1b, **unreduced**) |
| `chi_epsilon_reduced.txt` | same with `reduce=true` → ĵ^{-1} not ĵ^{-2} |
| `G1_from_chi.txt` | χ then 4-term fold with `P(i/j)`, `P(k/l)` |
| `G1_direct.txt` | unfactored ΩΩΓ strings (need restore `(1−P)` for full) |

### Locked equations (re-run AMC)

**χ^ε unreduced** (`output/chi_epsilon.tex`):

\[
\chi^\varepsilon_{ij}{}^{0}
=
\tfrac12\,\delta_{j_j j_i}\,\hat\jmath_i^{-2}
\sum_{abc\,J_0 J_1\lambda}
w\,(-1)^{J_0+J_1+\lambda}\,\hat\lambda^{-1}\,
\Omega^{J_0 J_1\lambda}_{ciab}\,\Omega^{J_1 J_0\lambda}_{abcj}
\]

\(w=\bar n_a\bar n_b n_c+n_a n_b\bar n_c\).

**Fold** (`G1_from_chi.tex`): four terms = \((1-P_{ij})\chi_{ia}\Gamma_{ajkl}+(1-P_{kl})\chi_{ak}\Gamma_{ijal}\).

**Packaging:** ethS / production uses **unreduced** χ (`ĵ^{-2}`).  
`reduce=true` on χ alone changes the print to `ĵ^{-1}` — do not mix.

### Status

| Step | Result |
|---|---|
| m ≡ AMC χ^ε | **PASS** `run/test_chi_epsilon_mscheme.py` |
| AMC PathB (χ×Γ) ≡ ethS slow ≡ DGEMM | **PASS** `run/test_tts_GI_mscheme.py` |
| ethS ≡ tts_GI | **PASS** `run/test_tts_GI.py` (emax=2) |
| DGEMM ≡ slow | **PASS** `run/test_tts_GI_dgemm.py` |

DGEMM: ordinary-channel normalized `GetMatrix`, leading `2×ang`, outer √(1+δ) on extract (LESSONS §Γ^I).

---

## Γ^II / χ^ζ

### AMC inputs (`factored_GII/input/`)

| File | Role |
|---|---|
| `chi_zeta.txt` | χ^ζ tensor 1b (Γ×Ω) |
| `chi_zeta_reduced.txt` | Γ `reduce=true` — drops one Ĵ̂ vs unreduced |
| `G2_direct.txt` | expand of Path B with `P(i/j)`/`P(k/l)`; unreduced ≡ m |
| `G2_from_chi.txt` | **fails** current AMC (Yutsis cannot orient χ^λ×Ω^λ→Γ) |

### Locked χ^ζ (`output/chi_zeta.tex`)

\[
\chi^\zeta_{ij}{}^{\lambda}
=
\tfrac12\,(-1)^{j_j+\lambda}
\sum_{abc\,J_0 J_1}
w\,(-1)^{J_0+j_a}\,
\hat J_0\hat J_1
\begin{Bmatrix}\lambda & J_1 & J_0 \\ j_a & j_i & j_j\end{Bmatrix}
\Gamma^{J_0}_{aibc}\,\Omega^{J_0 J_1\lambda}_{bcaj}
\]

\(w=n_a n_b\bar n_c+\bar n_a\bar n_b n_c\).  
Index order ≡ analyze \(\Gamma_{ciab}\Omega_{abcj}\) by rename.

Phase in code: use integer form \((-1)^{(j_{2j}+j_{2a})/2+\lambda+J_0}\) (same lesson as χ^β).

### Fold (analyze / sample — AMC print unavailable)

\[
\Gamma^{\mathrm{II}}
=
\sum_a\Bigl[
(1-P_{ij})\,\chi^\zeta_{aj}\,\Omega_{iakl}
-
(1-P_{kl})\,\chi^\zeta_{ak}\,\Omega_{ijal}
\Bigr]
\]

with \([\chi^\lambda\times\Omega^\lambda]^{(0)}\) when the LHS is scalar.

### Status

| Step | Result |
|---|---|
| m ≡ AMC χ^ζ | **PASS** `run/test_chi_zeta_mscheme.py` |
| AMC χ×Ω→Γ print | **done via RME** — use `χ_ja Ω_iakl` (not analyze `χ_aj`); [factored_GII/NOTES.md](factored_GII/NOTES.md) |
| Path B RME ≡ m | **PASS** `run/test_tts_GII_pathB_mscheme.py` (m: \(+R_1-R_2-R_3+R_4\)) |
| AMC DIRECT ≡ Path B ≡ m | **PASS** `run/test_tts_GII_direct_mscheme.py` |
| ethS GII λ=0 slow ≡ DGEMM | **PASS** `run/test_tts_GII_mscheme.py` |
| ethS GII λ=0 vs Factorized | **PASS** `run/test_tts_GII.py` (not vs tts_GII) |
| ethS GII λ≠0 | **not implemented** (`do_GII` gated) |

---

## Re-run AMC

```bash
amc -o factored_GI/output/chi_epsilon.tex factored_GI/input/chi_epsilon.txt
amc -o factored_GI/output/G1_from_chi.tex factored_GI/input/G1_from_chi.txt
amc -o factored_GII/output/chi_zeta.tex factored_GII/input/chi_zeta.txt
amc -o factored_GII/output/G2_direct.tex factored_GII/input/G2_direct.txt
```
