# Γ^{IV_b} AMC equation chain (Pandya→normal→RC→GEMM→Inv)

**WRAP-UP (full \(Z=W_1-W_2\)):** \(\mathrm{m}\equiv\mathrm{AMC\,direct}\equiv\mathrm{Path\,B}\) any λ.
Bench: `run/test_G4b_pathB_fold_mscheme.py`. χ links: `test_chi_iota_{m_vs_amc,pathB_vs_direct}.py`.

**Discipline:** AMC = angular only. Pack / \((1-P)^2\) / read-\(W_{jlik}\) restored by hand.
**Gold:** m-scheme / Path B fold. **No TTS.**

Official path: Pandya \(\bar\chi\) → normal \(\chi\) → RC \(\overline{\overline\chi}\) → \(W\) in Pandya → InvPandya.

Inputs: `input/steps/` · Outputs: `output/steps/`.

| Step | AMC input | Hand restore | Status |
|---|---|---|---|
| 0 Pandya \(\bar\chi\) / inv → normal \(\chi\) | (locked Path B) | — | **DONE** ≡ m ≡ AMC |
| 1 Pandya \(\Omega\) | `01_Omega_pandya_scalar.txt` | — | **AMC OK** |
| 3b normal→CC \(\chi\) | `03b_normal_to_CC_{scalar,tensor}.txt` | pack \((\mathrm{abik}-h_Z\,\mathrm{ikab})\) | **AMC OK** |
| 3b′ normal→CC \(\Omega\) (alt) | `03b_Omega_normal_to_CC_scalar.txt` | — | **AMC OK** |
| 4 GEMM | `04_gemm_IVb_pandya_x_RC.txt` | pack on \(\overline{\overline\chi}\) | **AMC OK** |
| 4′ GEMM both-CC | `04_gemm_IVb_both_CC.txt` | pack | **AMC OK** (alt) |
| 4A Path A term1 RC×RC | `04_gemm_IVb_both_CC_tensor*.txt` | replace \(\hat\lambda_0\to c(\lambda,\lambda_0)\) | **LOCKED** bare \(W_1\) |
| 4B Factorized term2 | pack \(\mathrm{RC}[-h_Z\bar\chi^{T}]\) | InvPandya \((1-P)^2\) | **LOCKED** λ=0 \(Z_2\) |
| 4C Path A term2 klij | reuse \(D_1\) @ \((klij,J_1\leftrightarrow J_0)\) | \(h_\Omega(-1)^{J_0-J_1}f\) | **LOCKED** all λ |
| 5 InvPandya | `05_inv_noperm.txt` | read \(W\) at \(jlik\) + \((1-P)^2\) | **AMC OK** |
| 5A InvRC | `05_invRC_*.txt` | — | **LOCKED** with 4A |

### Path A term1 (alt dual) — bare \(W_1\)

AMC prints \(\hat\lambda_0\) in the RC×RC product; **replace by** \(c(\lambda,\lambda_0)\)
(from RC–WE conjugation; odd \(\lambda_0\Rightarrow c=0\); \(c(\lambda,0)=1/\hat\lambda\)).
Then \(\mathrm{InvRC}[D]\equiv\mathrm{extract}(W_1)\) and \(\mathrm{WE}\equiv P(W_1)\).
Bench: `run/test_G4b_pathA_term1_correct.py`.

### Factorized term2 — \(Z_2\) from Pandya Ω × RC pack transpose (λ=0)

\[
\mathrm{RC}_2=\mathrm{RC}[-h_Z\bar\chi^{T}],\quad
\bar W_2=\bar\Omega\cdot\mathrm{RC}_2,\quad
Z_2=\mathrm{InvPandya}_{(1-P)^2}[\bar W_2]\equiv(1-P)^2(-W_2).
\]

Transpose lives on Pandya \(\bar\chi\) inside the pack, **not** on finished \(\mathrm{RC}_1\).
Bench: `run/test_G4b_term2_pandya_RC.py`.

### Path A term2 — klij remap of term1 (**LOCKED** all λ)

\[
W_2^{J_0J_1\lambda_0}_{ijkl}
=h_\Omega(-1)^{J_0-J_1}f(\lambda,\lambda_0)\,
\mathrm{InvRC}[D_1]^{J_1J_0\lambda_0}_{klij},
\]

\[
f=\frac{\sum_{M_0M_1}\mathrm{CG}(\lambda M_1,\lambda_0;\lambda M_0)}
{\sum_{M_0M_1}(-1)^{M_1}\mathrm{CG}(\lambda M_1,\lambda_0;\lambda M_0)}.
\]

No second GEMM. Origin: Ω† in the same RC–WE algebra as \(c(\lambda,\lambda_0)\).
Bench: `run/test_G4b_pathA_term2_klij.py`.
Rectangular Factorized \(\mathrm{RC}_2\) dual still open; fold remains alternate gold.

---

## Step 3b — normal \(\chi^\iota\) → CC (scalar)

\[
\overline{\overline\chi}{}^{\iota\,J}_{ijkl}
=
(-1)^{j_k+j_l}
\sum_{J'}(-1)^{J'}\,\hat J'^2
\begin{Bmatrix} j_j & j_l & J \\ j_k & j_i & J' \end{Bmatrix}
\chi^{\iota\,J'}_{ijkl}.
\]

Scheme `((1,-3),(4,-2))`. **Pack by hand** (not in AMC):
\[
\mathrm{RC}_{ab,ik}
=
\overline{\overline\chi}{}_{abik}
-h_Z\,\overline{\overline\chi}{}_{ikab}.
\]

## Step 4 — GEMM (Pandya \(\bar\Omega\) × RC \(\chi\))

AMC (`04_gemm_IVb_pandya_x_RC.tex`), pack stripped:

\[
W^{J}_{jlik}
=
(-1)^{J+j_j+j_l}
\sum_{ab\,J_2 J_3}
(-1)^{J_2+J_3+j_a+j_b}\,
\hat J_2^2\hat J_3^2
\begin{Bmatrix} J & J_2 & J_3 \\ j_b & j_k & j_j \end{Bmatrix}
\begin{Bmatrix} J & J_3 & J_2 \\ j_a & j_l & j_i \end{Bmatrix}
\bar\Omega^{J_2}_{jlab}\,
\overline{\overline\chi}{}^{\iota\,J_3}_{abik}.
\]

**Restore pack:** replace \(\overline{\overline\chi}_{abik}\) by \(\overline{\overline\chi}_{abik}-h_Z\overline{\overline\chi}_{ikab}\).

## Step 5 — InvPandya (noperm)

AMC:
\[
\Gamma^{J}_{ijkl}
=
-\sum_{J'}\hat J'^2
\begin{Bmatrix} j_l & j_k & J \\ j_j & j_i & J' \end{Bmatrix}
W^{J'}_{ijkl}.
\]

**Restore by hand** (Factorized L1884–2013): evaluate \(W\) at \(jl,ik\) (and three exchange slots) with full \((1-P_{ij})(1-P_{kl})\), plus code phases / \(\sqrt2\).

---

## Not the lock path

Factorized Pandya→RC one-shot (`03a_*`, `test_chi_iota_rc.py`) remains a **code dual** only.
Full dual analysis (Path 1 normal→RC vs Path 2 Pandya→RC, tensor continuous): [`RC_DUAL.md`](RC_DUAL.md).
