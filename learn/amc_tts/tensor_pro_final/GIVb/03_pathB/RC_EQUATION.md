# χ^ι Pandya → RC — follow scalar **code**, then AMC, then tensor

Discipline ([factored_GIIIb/NOTES.md](../../factored_GIIIb/NOTES.md)):
AMC = angular coupling only. Strip pack → AMC → restore pack by hand.

## 1. Code identity (`bar_CHI_V_RC`)

`FactorizedDoubleCommutator.cc` L1833–1841 (ignore stale comment L1752–1754):

| | |
|---|---|
| OUT | `bar_CHI_V_RC`: bra \((a,b)\), ket \((c,d)\), channel \(J\) |
| IN | `bar_CHI_V`: \((a,d;\,b,c)\) and \((b,c;\,a,d)\), channel \(J'\) |
| 6j | \(\mathrm{GetSixJ}(j_a,j_b,J,j_c,j_d,J')\) |
| phase | `phase((jb+jc)/2 + J')` |
| overall | `+=` (ι); η uses `-=` and \(+\) pack |

**Full code equation (implement this):**

\[
\mathrm{RC}[\bar\chi^\iota]^{J}_{ab,cd}
=
\sum_{J'}(2J'+1)\,(-1)^{j_b+j_c+J'}
\begin{Bmatrix} j_a & j_b & J \\ j_c & j_d & J' \end{Bmatrix}
\bigl(
  \bar\chi^{\iota\,J'}_{ad,bc}
  -h_Z\,\bar\chi^{\iota\,J'}_{bc,ad}
\bigr).
\]

## 2. AMC skeleton (pack stripped)

Same-label `dbarChi_ijkl = barChi_ijkl`, schemes
`((1,-4),(3,-2)) → ((1,-3),(4,-2))`.

| | input | output |
|---|---|---|
| scalar | [`03a_RC_scheme_change_scalar.txt`](../input/steps/03a_RC_scheme_change_scalar.txt) | [`03a_RC_scheme_change_scalar.tex`](../output/steps/03a_RC_scheme_change_scalar.tex) |
| tensor | [`03a_RC_scheme_change_tensor.txt`](../input/steps/03a_RC_scheme_change_tensor.txt) | [`03a_RC_scheme_change_tensor.tex`](../output/steps/03a_RC_scheme_change_tensor.tex) (`--collect-ninejs`) |

**Scalar AMC:**

\[
\overline{\overline\chi}{}^{\iota\,J}_{ijkl}
=
(-1)^{J+j_k+j_l}
\sum_{J'}(-1)^{J'}\,\hat J'^2
\begin{Bmatrix} j_l & j_j & J \\ j_k & j_i & J' \end{Bmatrix}
\bar\chi^{\iota\,J'}_{ijkl}.
\]

**Tensor AMC:**

\[
\overline{\overline\chi}{}^{\iota\,J_0 J_1\lambda}_{ijkl}
=
(-1)^{J_0+J_1+j_i+j_j+j_k+j_l+\lambda}\,
\hat J_0\hat J_1
\sum_{J_2 J_3}
(-1)^{J_2+J_3}\,
\hat J_2\hat J_3
\begin{Bmatrix}
\lambda & J_0 & J_1 \\
J_3 & j_k & j_j \\
J_2 & j_i & j_l
\end{Bmatrix}
\bar\chi^{\iota\,J_2 J_3\lambda}_{ijkl}.
\]

Restore pack on the **code** equation (§1), not by feeding \((ad,bc)\) into AMC.

## 3. Tensor Path B RC (same procedure)

Use §1 with rectangular Pandya \(\bar\chi^{\iota\,J_2 J_3\lambda}\) (locked Path B).
At \(\lambda=0\) this is bit-identical to scalar Factorized.
For \(\lambda\neq 0\): keep code pack/rewire; angular factor = tensor continuous of
the **code** 6j (AMC §2 ninej is the Oc=Op twin — validate against Factorized λ=0
before trusting as drop-in).

Bench: `run/test_chi_iota_rc.py` — Path B \(\bar\chi\) → code RC.
