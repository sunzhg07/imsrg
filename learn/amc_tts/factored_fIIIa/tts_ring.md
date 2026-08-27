# TTS ring / Path A ≡ Path B

Ring (ph) product of two reduced tensors → reduced scalar:

\[
Z_{pqsr}
=
\sum_{ab}
\Omega_{p\bar{a}rb}\,
\Omega_{aqsb}
\qquad(\lambda_Z=0).
\]

AMC inputs: `input/z_pbar_aqsb_direct.txt`, `input/z_pbar_aqsb_pandya.txt`.

## Verdict (locked by m-scheme)

**Path A (AMC direct) is correct. Path B as printed = −Path A.**

Benchmark: `run/test_z_ring_mscheme_sign.py` (random reduced \(\Omega^{\lambda=2}\)).

| Comparison | Result |
|---|---|
| Kernel \(K_A\) vs \(K_B\) | \(K_A=K_B\) exactly |
| Full AMC Path A vs Path B | \(A=-B\) |
| **m-scheme oracle vs Path A** | **equal** (\(\max\|\Delta\|\sim10^{-15}\)) |
| m-scheme oracle vs Path B | \(B=-Z_m\) |

### m-scheme oracle (truth)

Naive \(\sum_{ab}\Omega\Omega\) is **not** enough: two rank-\(\lambda\) tensors → scalar need

\[
Z(m)=\sum_{a m_a b m_b}
\begin{pmatrix}\text{CG}\end{pmatrix}\!(\lambda\,\mu,\,\lambda\,{-\mu};\,0\,0)\,
\Omega_{pb,ar}(m)\,\Omega_{aq,sb}(m)
=\sum\frac{(-1)^{\lambda-\mu}}{\hat\lambda}\,\Omega\,\Omega,
\]

with \(\mu=(m_p+m_b)-(m_a+m_r)\). Then project to **reduced** scalar (\(Z\) AMC `reduce=true`):

\[
Z^{J}_{\mathrm{red}}
=\frac{1}{\hat J}\sum_{m}\mathrm{CG}_{pq}\,\mathrm{CG}_{sr}\,Z(m),
\qquad
Z_{\mathrm{unred}}=Z_{\mathrm{red}}/\hat J.
\]

### Practical

- Trust **Path A** / fix Path B by dropping one overall minus (inv Pandya or mid).
- \(Z\) in these AMC inputs is **reduced**; do not mix with unreduced storage without \(\hat J\).

## J ↔ m convention (`UnitTest::GetMschemeMatrixElement_2b`)

m-scheme MEs are **always unreduced** (physical \(\langle m|\cdots|m\rangle\)); there is no “reduced m-scheme.”

| J storage | Unpack to m |
|---|---|
| **Unreduced** scalar (`IsReduced()=false`) | 2 CGs: uncouple bra + ket |
| **Reduced** tensor/scalar (`IsReduced()=true`) | 3 CGs: WE \(\mathrm{CG}(J_2 M_2,\lambda\mu;J_1 M_1)/\hat J_1\) + uncouple bra + ket |

Unreduced scalar:
\[
\langle ab\,m|O|cd\,m\rangle
=\sum_J\mathrm{CG}_{ab}\,\mathrm{CG}_{cd}\,O^{J}.
\]

Reduced (WE + uncouple):
\[
\langle ab\,m|O^\lambda|cd\,m\rangle
=\sum_{J_1 J_2}
\frac{\mathrm{CG}(J_2 M_2,\lambda\mu;J_1 M_1)}{\hat J_1}
\,\mathrm{CG}_{ab}\,\mathrm{CG}_{cd}\,
\langle J_1\|O^\lambda\|J_2\rangle.
\]

**Norm of a reduced J-scheme operator** (m-scheme–equivalent): unreduce first, then weight by magnetic degeneracy, e.g. for a scalar
\[
\|O\|_m \sim \sqrt{\sum_J (2J+1)\,|O^{J}_{\mathrm{unred}}|^2}
=\sqrt{\sum_J |O^{J}_{\mathrm{red}}|^2},
\]
since \(O_{\mathrm{red}}=\hat J\,O_{\mathrm{unred}}\). Plain `TwoBodyNorm()` does **not** apply this weighting by itself.

### Tests

| Script | What it checks |
|---|---|
| `run/test_omega_mscheme_roundtrip.py` | Random \(\Omega^{\lambda=2}\): J→m→J RME round-trip |
| `run/test_z_ring_mscheme_sign.py` | Ring \(Z\): m-oracle vs Path A / Path B (sign lock) |

## Angular-kernel identity

Path A (direct ninej form): one intermediate \(j_0\), **two 6j + one 9j**:

\[
K_A
=
\sum_{j_0}
\hat{j}_0^{\,2}\,
\begin{Bmatrix}
J_\beta & \lambda & J_\alpha \\
j_b & j_p & j_0
\end{Bmatrix}
\begin{Bmatrix}
J_\gamma & \lambda & J_\delta \\
j_b & j_s & j_0
\end{Bmatrix}
\begin{Bmatrix}
j_r & j_a & J_\beta \\
j_s & J_\gamma & j_0 \\
J_0 & j_q & j_p
\end{Bmatrix}.
\]

Path B: sum over Pandya intermediates \(J_p,J_2\), **one inv 6j + two Pandya 9js**:

\[
K_B
=
\sum_{J_p,J_2}
\hat{J}_p^{\,2}\,\hat{J}_2^{\,2}\,
\begin{Bmatrix}
j_r & j_s & J_0 \\
j_q & j_p & J_p
\end{Bmatrix}
\begin{Bmatrix}
\lambda & J_p & J_2 \\
J_\beta & j_r & j_a \\
J_\alpha & j_p & j_b
\end{Bmatrix}
\begin{Bmatrix}
\lambda & J_2 & J_p \\
J_\delta & j_b & j_s \\
J_\gamma & j_a & j_q
\end{Bmatrix}.
\]

Checked numerically (\(\lambda=1,2\), dense half-integer scan):

\[
\boxed{K_A = K_B}
\qquad(\max|K_A-K_B|\sim 10^{-16}).
\]

No phase or factor inside this identity. Summing \(J_p,J_2\) **reduces** Path B to Path A’s \(1\times9\mathrm{j}+2\times6\mathrm{j}\) form.

Full AMC \(Z\): Path A matches m-scheme; Path B is overall minus of Path A.

\[
Z^{\mathrm{A}}_{\mathrm{AMC}} = Z_{m}
\qquad\text{and}\qquad
Z^{\mathrm{A}}_{\mathrm{AMC}} = -\,Z^{\mathrm{B}}_{\mathrm{expanded}}.
\]

## Related

- Pandya ≠ ph: Pandya is recoupling; occupations enter only in the RME midstep.
- Production \(\chi^\gamma\) / \(f^{\mathrm{III}_a}\): Path B packaging — see `NOTES.md`.
- Broader Ω×Ω cases: `OMEGA_TT_TO_SCALAR.md`, `OMEGA_CROSS_NOOCC.md`.
