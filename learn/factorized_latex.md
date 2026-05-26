# Factorized double-commutator equations

This note extracts the equations implemented in [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc) and [src/FactorizedDoubleCommutator.cc](../src/FactorizedDoubleCommutator.cc) for

$$
Z = [\eta,[\eta,\Gamma]_{3b}]_{1b,2b},
$$

where all two-body matrix elements are the non-reduced J-coupled matrix elements used by the code.

## Notation

$$
\bar n_a \equiv 1-n_a, \qquad d_a \equiv 2j_a+1, \qquad \hat J^2 \equiv 2J+1.
$$

The C++ code often stores `j2 = 2j`, so factors such as `op.j2 + 1.0` are written here as $d_p = 2j_p+1$.  Matrix elements are written as

$$
\eta^J_{ab,cd} \equiv \langle ab;J | \eta | cd;J\rangle,
\qquad
\Gamma^J_{ab,cd} \equiv \langle ab;J | \Gamma | cd;J\rangle .
$$

The permutation symbols $P_{pg}$ and $P_{qh}$ mean that the code evaluates the direct term plus the exchanged terms with the corresponding stored two-body phase, e.g. `bra.Phase(J0)` or `ket.Phase(J0)`.  The final matrix storage also divides by $\sqrt2$ when either pair contains identical orbits.

## Reference implementation: `comm223_231_BruteForce`

Source: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L6620-L7352).

This routine builds the one-body part $Z_{pq}$ directly by summing the full diagram expressions.

### Diagram I

Source: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L6631-L6725).

$$
\begin{aligned}
Z^{I}_{pq}
&= \frac{1}{2}\frac{\delta_{j_p j_q}}{d_p}
\sum_{abcdeJ_0J_1}\delta_{j_dj_e}
\frac{\hat J_0^2\hat J_1^2}{d_d}
\Big(
\bar n_a\bar n_c n_b n_d
-\bar n_b\bar n_d n_a n_c
-\bar n_b\bar n_e n_a n_c
+\bar n_a\bar n_c n_b n_e
\Big) \\
&\qquad\times
\eta^{J_0}_{bd,ac}\,
\eta^{J_0}_{ac,be}\,
\Gamma^{J_1}_{ep,dq} .
\end{aligned}
$$

### Diagram IIa

Source: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L6731-L6847).

$$
\begin{aligned}
Z^{IIa}_{pq}
&=\frac{\delta_{j_pj_q}}{d_p}
\sum_{abcdeJ_0J_1J_2J_3}
(-1)^{J_0+J_1+J_2+j_c+j_d}
\hat J_0^2\hat J_1^2\hat J_2^2\hat J_3^2
\Big(\bar n_a\bar n_c n_b n_d-\bar n_b\bar n_d n_a n_c\Big) \\
&\qquad\times
\begin{Bmatrix} j_a&j_b&J_3\\ j_d&j_c&J_0\end{Bmatrix}
\begin{Bmatrix} j_p&j_e&J_3\\ j_d&j_c&J_1\end{Bmatrix}
\begin{Bmatrix} j_b&j_p&J_2\\ j_e&j_a&J_3\end{Bmatrix}
\eta^{J_0}_{bd,ac}\,
\eta^{J_1}_{cp,de}\,
\Gamma^{J_2}_{ae,bq} .
\end{aligned}
$$

### Diagram IIb

Source: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L6853-L6957).

$$
\begin{aligned}
Z^{IIb}_{pq}
&=\frac{1}{4}\frac{\delta_{j_pj_q}}{d_p}
\sum_{abcdeJ_0}
\hat J_0^2
\Big(\bar n_a\bar n_d n_b n_e-\bar n_b\bar n_e n_a n_d\Big)
\eta^{J_0}_{be,ad}\,
\eta^{J_0}_{cp,be}\,
\Gamma^{J_0}_{ad,cq} .
\end{aligned}
$$

### Diagram IIc

Source: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L6963-L7078).

$$
\begin{aligned}
Z^{IIc}_{pq}
&=\frac{\delta_{j_pj_q}}{d_p}
\sum_{abcdeJ_0J_1J_2J_3}
\hat J_0^2\hat J_1^2\hat J_2^2\hat J_3^2
\Big(\bar n_a\bar n_d n_b n_e-\bar n_b\bar n_e n_a n_d\Big) \\
&\qquad\times
\begin{Bmatrix} j_d&j_e&J_3\\ j_b&j_a&J_0\end{Bmatrix}
\begin{Bmatrix} j_c&j_p&J_3\\ j_b&j_a&J_1\end{Bmatrix}
\begin{Bmatrix} j_e&j_c&J_2\\ j_p&j_d&J_3\end{Bmatrix}
\eta^{J_0}_{be,da}\,
\eta^{J_1}_{ca,bq}\,
\Gamma^{J_2}_{dp,ce} .
\end{aligned}
$$

### Diagram IId

Source: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L7084-L7163).

The code keeps this as a commented expression and notes that it is almost the same as IIb, so it is combined in the optimized implementation.

$$
\begin{aligned}
Z^{IId}_{pq}
&=-\frac{1}{4}\frac{\delta_{j_pj_q}}{d_p}
\sum_{abcdeJ_0}
\hat J_0^2
\Big(\bar n_c\bar n_d n_a n_e-\bar n_a\bar n_e n_c n_d\Big)
\eta^{J_0}_{ae,cd}\,
\eta^{J_0}_{cd,bq}\,
\Gamma^{J_0}_{bp,ae} .
\end{aligned}
$$

### Diagram IIIa

Source: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L7171-L7256).

$$
\begin{aligned}
Z^{IIIa}_{pq}
&=\frac{1}{2}\frac{\delta_{j_pj_q}}{d_p}
\sum_{abcdeJ_0J_1}\delta_{j_dj_e}
\frac{\hat J_0^2\hat J_1^2}{d_d}
\Big(\bar n_a\bar n_e n_b n_c-\bar n_b\bar n_c n_a n_e\Big)
\eta^{J_0}_{bc,ae}\,
\eta^{J_1}_{ep,dq}\,
\Gamma^{J_0}_{ad,bc} .
\end{aligned}
$$

### Diagram IIIb

Source: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L7262-L7352).

$$
\begin{aligned}
Z^{IIIb}_{pq}
&=-\frac{1}{2}\frac{\delta_{j_pj_q}}{d_p}
\sum_{abcdeJ_0J_1}\delta_{j_dj_e}
\frac{\hat J_0^2\hat J_1^2}{d_d}
\Big(\bar n_a\bar n_c n_b n_d-\bar n_b\bar n_d n_a n_c\Big)
\eta^{J_0}_{bd,ac}\,
\eta^{J_1}_{ep,dq}\,
\Gamma^{J_0}_{ac,be} .
\end{aligned}
$$

## Factorized implementation: `comm223_231`

Source: [src/FactorizedDoubleCommutator.cc](../src/FactorizedDoubleCommutator.cc#L55-L801).

ArXiv correspondence: this section matches the one-body J-coupled factorized equations in [learn/arxiv.txt](./arxiv.txt), namely
$$
f^{\mathrm{I}},\quad f^{\mathrm{II}},\quad f^{\mathrm{III}_a},\quad f^{\mathrm{III}_b}
$$
in `doubleCommutator_1b_Jcoupled_factorized`, together with the intermediates
$$
\chi^{\alpha},\quad \chi^{\beta},\quad \bar\chi^{\gamma},\quad \chi^{\delta}
$$
in `Jcoupled_chi_1b`.

The optimized one-body implementation splits the same result into a one-body-intermediate path and a two-body-intermediate path.

### One-body intermediate, diagram I

Source: [src/FactorizedDoubleCommutator.cc](../src/FactorizedDoubleCommutator.cc#L99-L249).

ArXiv correspondence tag: this is the code-backed version of
$$
f^{\mathrm{I}}_{ij}
$$
with intermediate
$$
\chi^{\alpha}_{ij}.
$$

For a two-body channel with angular momentum $J_0$, define the occupation-dressed matrix

$$
{\color{red}\text{ArXiv correspondence: contributes to }\chi^{\alpha}\text{ in }\texttt{Jcoupled\_chi\_1b}.}
$$

$$
\left(\eta_{nn\bar n\bar n}^{J_0}\right)_{ij,kl}
=\Big[n_i n_j\bar n_k\bar n_l-\bar n_i\bar n_j n_k n_l\Big] \eta^{J_0}_{ij,kl} .
$$

The code forms

$$
{\color{red}\text{ArXiv correspondence: matrix-product realization of }\chi^{\alpha}.}
$$

$$
T^{J_0}=2\hat J_0^2\,\eta_{nn\bar n\bar n}^{J_0}\eta^{J_0}
+\left(2\hat J_0^2\,\eta_{nn\bar n\bar n}^{J_0}\eta^{J_0}\right)^T,
$$

and contracts it to

$$
{\color{red}\text{ArXiv correspondence: }\chi^{\alpha}_{ij}.}
$$

$$
\chi^{(221a)}_{de}
=\frac{1}{d_d}\sum_{bJ_0}T^{J_0}_{bd,be} .
$$

Then

$$
{\color{red}\text{ArXiv correspondence: }f^{\mathrm{I}}_{ij}\text{ in }\texttt{doubleCommutator\_1b\_Jcoupled\_factorized}.}
$$

$$
Z^{I}_{pq}=\frac{1}{2d_p}\sum_{deJ_1}\hat J_1^2\,
\chi^{(221a)}_{de}\,
\Gamma^{J_1}_{ep,dq} .
$$

### One-body intermediate, diagrams IIIa and IIIb

Source: [src/FactorizedDoubleCommutator.cc](../src/FactorizedDoubleCommutator.cc#L257-L371).

ArXiv correspondence tag: this is the code-backed version of
$$
f^{\mathrm{II}}_{ij}
$$
with intermediate
$$
\chi^{\beta}_{ij}.
$$

The non-Hermitian one-body intermediate is

$$
{\color{red}\text{ArXiv correspondence: }\chi^{\beta}_{ij}\text{ in }\texttt{Jcoupled\_chi\_1b}.}
$$

$$
\chi^{(221b)}_{de}
=\frac{1}{d_d}
\sum_{abcJ_0}\hat J_0^2
\Big(\bar n_a\bar n_e n_b n_c-\bar n_b\bar n_c n_a n_e\Big)
\eta^{J_0}_{bc,ae}\,
\Gamma^{J_0}_{ad,bc} .
$$

The final contraction is

$$
{\color{red}\text{ArXiv correspondence: }f^{\mathrm{II}}_{ij}\text{ in }\texttt{doubleCommutator\_1b\_Jcoupled\_factorized}.}
$$

$$
Z^{IIIa+IIIb}_{pq}
=\frac{1}{2d_p}\sum_{deJ_1}\hat J_1^2
\Big[\chi^{(221b)}_{de}-h_Z\chi^{(221b)}_{ed}\Big]
\eta^{J_1}_{ep,dq} .
$$

### Two-body intermediate, diagrams IIb and IId

Source: [src/FactorizedDoubleCommutator.cc](../src/FactorizedDoubleCommutator.cc#L380-L541).

ArXiv correspondence tag: this is the code-backed version of
$$
f^{\mathrm{III}_b}_{ij}
$$
with intermediate
$$
\chi^{\delta}_{ijkl}.
$$

In each ordinary two-body channel, define

$$
{\color{red}\text{ArXiv correspondence: contributes to }\chi^{\delta J}_{ijkl}\text{ in }\texttt{Jcoupled\_chi\_1b}.}
$$

$$
\eta^{J}_{occ}(ij,kl)
=\Big[n_i n_j\bar n_k\bar n_l-\bar n_i\bar n_j n_k n_l\Big]\eta^J_{ij,kl} .
$$

The optimized code builds the matrix product

$$
{\color{red}\text{ArXiv correspondence: matrix-product realization of }\chi^{\delta J}.}
$$

$$
\chi^{(222b),J}
=4\hat J^2\left(
\eta^J\eta^J_{occ}\Gamma^J
-\Gamma^J\eta^J_{occ}\eta^J
\right),
$$

with the equal-channel case using transpose symmetry.  The one-body contraction is

$$
{\color{red}\text{ArXiv correspondence: }f^{\mathrm{III}_b}_{ij}\text{ in }\texttt{doubleCommutator\_1b\_Jcoupled\_factorized}.}
$$

$$
Z^{IIb+IId}_{pq}
=\frac{1}{4d_p}\sum_{cJ_0}\chi^{(222b),J_0}_{cp,cq} .
$$

### Pandya transform for diagrams IIa and IIc

Source: [src/FactorizedDoubleCommutator.cc](../src/FactorizedDoubleCommutator.cc#L546-L801).

ArXiv correspondence tag: this is the code-backed version of
$$
f^{\mathrm{III}_a}_{ij}
$$
with intermediate
$$
\bar\chi^{\gamma J}_{i\bar j k\bar l},
$$
written in the Pandya/cross-coupled channel used by the code.

The code uses cross-coupled channels and the Pandya transform

$$
{\color{red}\text{ArXiv correspondence: Pandya transform entering }\bar\chi^{\gamma J}.}
$$

$$
\bar X^J_{ab,cd}
=\sum_{J'}(-1)^{j_b+j_c+J'}\hat J'^2
\begin{Bmatrix} j_a&j_b&J\\ j_c&j_d&J'\end{Bmatrix}
X^{J'}_{ad,bc} .
$$

With

$$
{\color{red}\text{ArXiv correspondence: occupation factor entering }\bar\chi^{\gamma J}.}
$$

$$
f_{ab,cd}=\bar n_c\bar n_b n_a n_d-n_c n_b\bar n_a\bar n_d,
$$

the intermediate matrix is

$$
{\color{red}\text{ArXiv correspondence: matrix-product realization of }\bar\chi^{\gamma J}_{i\bar j k\bar l}.}
$$

$$
\bar\chi^{(222a),J}=\hat J^2\,\bar\eta^J\,\bar\eta^J_{occ}\,\bar\Gamma^J .
$$

The one-body contractions are implemented as the two terms

$$
{\color{red}\text{ArXiv correspondence: }f^{\mathrm{III}_a}_{ij}\text{ in }\texttt{doubleCommutator\_1b\_Jcoupled\_factorized}.}
$$

$$
Z^{IIa}_{pq}=\frac{1}{d_p}\sum_{eJ}\bar\chi^{(222a),J}_{pe,qe},
\qquad
Z^{IIc}_{pq}=-\frac{1}{d_p}\sum_{eJ}\bar\chi^{(222a),J}_{eq,ep},
$$

with signs from the stored cross-coupled ordering.

## Reference implementation: `comm223_232_BruteForce`

Source: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L7359-L10303).

This routine builds the two-body part $Z^J_{pg,qh}$ directly.  The direct forms below use the permutation symbols that the code expands by explicit direct and exchange loops.

### Diagrams Ia, Ib, IVa, IVb

Source: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L7386-L7818).

$$
\begin{aligned}
Z^{Ia,J_0}_{pg,qh}
&=\frac{1}{2}P_{pg}\sum_{abcdJ_2}
\frac{\delta_{j_dj_p}}{d_p}\hat J_2^2
\Big(\bar n_a\bar n_c n_b+\bar n_b n_a n_c\Big)
\eta^{J_2}_{bp,ac}\,
\eta^{J_2}_{ac,bd}\,
\Gamma^{J_0}_{dg,qh},\\
Z^{Ib,J_0}_{pg,qh}
&=\frac{1}{2}P_{qh}\sum_{abcdJ_2}
\frac{\delta_{j_dj_q}}{d_q}\hat J_2^2
\Big(\bar n_a n_b n_c+\bar n_b\bar n_c n_a\Big)
\eta^{J_2}_{ad,bc}\,
\eta^{J_2}_{bc,aq}\,
\Gamma^{J_0}_{pg,dh},\\
Z^{IVa,J_0}_{pg,qh}
&=-\frac{1}{2}P_{qh}\sum_{abcdJ_2}
\frac{\delta_{j_dj_q}}{d_q}\hat J_2^2
\Big(\bar n_a n_b n_c+\bar n_b\bar n_c n_a\Big)
\eta^{J_2}_{bc,aq}\,
\eta^{J_0}_{pg,dh}\,
\Gamma^{J_2}_{ad,bc},\\
Z^{IVb,J_0}_{pg,qh}
&=-\frac{1}{2}P_{pg}\sum_{abcdJ_2}
\frac{\delta_{j_dj_p}}{d_p}\hat J_2^2
\Big(\bar n_a\bar n_c n_b+\bar n_b n_a n_c\Big)
\eta^{J_2}_{bp,ac}\,
\eta^{J_0}_{dg,qh}\,
\Gamma^{J_2}_{ac,bd}.
\end{aligned}
$$

### Diagrams IIa and IIc

Source: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L7823-L8135).

$$
\begin{aligned}
Z^{IIa,J_0}_{pg,qh}
&=-P_{pg}\sum_{abcdJ_2J_3J_4}
\hat J_2^2\hat J_3^2\hat J_4^2
\Big(\bar n_b\bar n_d n_c+\bar n_c n_b n_d\Big) \\
&\qquad\times
\begin{Bmatrix} j_d&j_g&J_4\\ j_c&j_b&J_2\end{Bmatrix}
\begin{Bmatrix} j_p&j_a&J_4\\ j_c&j_b&J_3\end{Bmatrix}
\begin{Bmatrix} j_g&j_p&J_0\\ j_a&j_d&J_4\end{Bmatrix}
\eta^{J_2}_{cg,db}\,
\eta^{J_3}_{pb,ca}\,
\Gamma^{J_0}_{da,qh},\\
Z^{IIc,J_0}_{pg,qh}
&=-P_{qh}\sum_{abcdJ_2J_3J_4}
\hat J_2^2\hat J_3^2\hat J_4^2
\Big(\bar n_b n_c n_d+\bar n_c\bar n_d n_b\Big) \\
&\qquad\times
\begin{Bmatrix} j_q&j_d&J_4\\ j_c&j_b&J_2\end{Bmatrix}
\begin{Bmatrix} j_a&j_h&J_4\\ j_c&j_b&J_3\end{Bmatrix}
\begin{Bmatrix} j_q&j_h&J_0\\ j_a&j_d&J_4\end{Bmatrix}
\eta^{J_2}_{cd,qb}\,
\eta^{J_3}_{ab,ch}\,
\Gamma^{J_0}_{pg,ad}.
\end{aligned}
$$

### Diagrams IIb and IId

Source: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L8140-L8628).

$$
\begin{aligned}
Z^{IIb,J_0}_{pg,qh}
&=-P_{pg}P_{qh}\sum_{abcdJ_2J_3J_4J_5}
\hat J_2^2\hat J_3^2\hat J_4^2\hat J_5^2
\Big(\bar n_b n_c n_d+\bar n_c\bar n_d n_b\Big) \\
&\qquad\times
\begin{Bmatrix} j_q&j_d&J_5\\ j_c&j_b&J_2\end{Bmatrix}
\begin{Bmatrix} j_p&j_a&J_5\\ j_c&j_b&J_3\end{Bmatrix}
\begin{Bmatrix} J_0&J_5&J_4\\ j_d&j_h&j_q\end{Bmatrix}
\begin{Bmatrix} J_0&J_4&J_5\\ j_a&j_p&j_g\end{Bmatrix}
\eta^{J_2}_{dc,bq}\,
\eta^{J_3}_{bp,ac}\,
\Gamma^{J_4}_{ga,hd},\\
Z^{IId,J_0}_{pg,qh}
&=-P_{pg}P_{qh}\sum_{abcdJ_2J_3J_4J_5}
\hat J_2^2\hat J_3^2\hat J_4^2\hat J_5^2
\Big(\bar n_c n_b n_d+\bar n_b\bar n_d n_c\Big) \\
&\qquad\times
\begin{Bmatrix} j_d&j_g&J_5\\ j_c&j_b&J_2\end{Bmatrix}
\begin{Bmatrix} j_a&j_h&J_5\\ j_c&j_b&J_3\end{Bmatrix}
\begin{Bmatrix} J_0&J_5&J_4\\ j_d&j_p&j_g\end{Bmatrix}
\begin{Bmatrix} J_5&J_4&J_0\\ j_q&j_h&j_a\end{Bmatrix}
\eta^{J_2}_{gc,bd}\,
\eta^{J_3}_{ba,hc}\,
\Gamma^{J_4}_{dp,aq}.
\end{aligned}
$$

### Diagrams IIe and IIf

Source: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L8633-L9061).

$$
\begin{aligned}
Z^{IIe,J_0}_{pg,qh}
&=-\frac{1}{2}P_{pg}P_{qh}\sum_{abcdJ_2J_3J_4}
\hat J_2^2\hat J_3^2\hat J_4^2
\Big(\bar n_b n_a n_c+\bar n_a\bar n_c n_b\Big) \\
&\qquad\times
\begin{Bmatrix} j_p&j_h&J_4\\ j_b&j_d&J_2\end{Bmatrix}
\begin{Bmatrix} j_q&j_g&J_4\\ j_b&j_d&J_3\end{Bmatrix}
\begin{Bmatrix} j_h&j_q&J_0\\ j_g&j_p&J_4\end{Bmatrix}
\eta^{J_2}_{ac,bh}\,
\eta^{J_2}_{pd,ac}\,
\Gamma^{J_3}_{bg,qd},\\
Z^{IIf,J_0}_{pg,qh}
&=-\frac{1}{2}P_{pg}P_{qh}\sum_{abcdJ_2J_3J_4}
\hat J_2^2\hat J_3^2\hat J_4^2
\Big(\bar n_a\bar n_c n_b+\bar n_b n_a n_c\Big) \\
&\qquad\times
\begin{Bmatrix} j_h&j_p&J_4\\ j_b&j_d&J_2\end{Bmatrix}
\begin{Bmatrix} j_g&j_q&J_4\\ j_b&j_d&J_3\end{Bmatrix}
\begin{Bmatrix} j_h&j_q&J_0\\ j_g&j_p&J_4\end{Bmatrix}
\eta^{J_2}_{pb,ac}\,
\eta^{J_2}_{ac,dh}\,
\Gamma^{J_3}_{dg,qb}.
\end{aligned}
$$

### Diagrams IIIa and IIIb

Source: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L9066-L9551).

$$
\begin{aligned}
Z^{IIIa,J_0}_{pg,qh}
&=P_{pg}P_{qh}\sum_{abcdJ_2J_3J_4J_5}
\hat J_2^2\hat J_3^2\hat J_4^2\hat J_5^2
\Big(\bar n_a\bar n_c n_b+\bar n_b n_a n_c\Big) \\
&\qquad\times
\begin{Bmatrix} j_a&j_b&J_5\\ j_p&j_c&J_2\end{Bmatrix}
\begin{Bmatrix} J_3&J_0&J_5\\ j_p&j_c&j_g\end{Bmatrix}
\begin{Bmatrix} j_q&j_d&J_5\\ j_a&j_b&J_4\end{Bmatrix}
\begin{Bmatrix} J_3&J_0&J_5\\ j_q&j_d&j_h\end{Bmatrix}
\eta^{J_2}_{bp,ca}\,
\eta^{J_3}_{gc,hd}\,
\Gamma^{J_4}_{da,bq},\\
Z^{IIIb,J_0}_{pg,qh}
&=P_{pg}P_{qh}\sum_{abcdJ_2J_3J_4J_5}
\hat J_2^2\hat J_3^2\hat J_4^2\hat J_5^2
\Big(\bar n_a n_b n_c+\bar n_b\bar n_c n_a\Big) \\
&\qquad\times
\begin{Bmatrix} j_a&j_b&J_5\\ j_c&j_q&J_2\end{Bmatrix}
\begin{Bmatrix} J_0&J_3&J_5\\ j_c&j_q&j_h\end{Bmatrix}
\begin{Bmatrix} j_d&j_p&J_5\\ j_a&j_b&J_4\end{Bmatrix}
\begin{Bmatrix} J_0&J_3&J_5\\ j_d&j_p&j_g\end{Bmatrix}
\eta^{J_2}_{cb,aq}\,
\eta^{J_3}_{gd,hc}\,
\Gamma^{J_4}_{ap,db}.
\end{aligned}
$$

### Diagrams IIIc and IIId

Source: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L9556-L9868).

$$
\begin{aligned}
Z^{IIIc,J_0}_{pg,qh}
&=-P_{qh}\sum_{abcdJ_2J_3J_4}
\hat J_2^2\hat J_3^2\hat J_4^2
\Big(\bar n_a n_b n_c+\bar n_b\bar n_c n_a\Big) \\
&\qquad\times
\begin{Bmatrix} j_q&j_b&J_4\\ j_c&j_a&J_2\end{Bmatrix}
\begin{Bmatrix} J_3&J_0&J_4\\ j_c&j_a&j_d\end{Bmatrix}
\begin{Bmatrix} J_4&J_3&J_0\\ j_h&j_q&j_b\end{Bmatrix}
\eta^{J_2}_{bc,aq}\,
\eta^{J_0}_{pg,cd}\,
\Gamma^{J_3}_{da,hb},\\
Z^{IIId,J_0}_{pg,qh}
&=-P_{pg}\sum_{abcdJ_2J_3J_4}
\hat J_2^2\hat J_3^2\hat J_4^2
\Big(\bar n_a\bar n_c n_b+\bar n_b n_a n_c\Big) \\
&\qquad\times
\begin{Bmatrix} j_a&j_p&J_4\\ j_b&j_c&J_2\end{Bmatrix}
\begin{Bmatrix} J_0&J_3&J_4\\ j_b&j_c&j_d\end{Bmatrix}
\begin{Bmatrix} J_4&J_3&J_0\\ j_g&j_p&j_a\end{Bmatrix}
\eta^{J_2}_{bp,ac}\,
\eta^{J_0}_{cd,qh}\,
\Gamma^{J_3}_{ga,db}.
\end{aligned}
$$

### Diagrams IIIe and IIIf

Source: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L9873-L10299).

$$
\begin{aligned}
Z^{IIIe,J_0}_{pg,qh}
&=-\frac{1}{2}P_{pg}P_{qh}\sum_{abcdJ_2J_3J_4}
\hat J_2^2\hat J_3^2\hat J_4^2
\Big(\bar n_d n_a n_b+\bar n_a\bar n_b n_d\Big) \\
&\qquad\times
\begin{Bmatrix} J_2&J_3&J_4\\ j_p&j_h&j_d\end{Bmatrix}
\begin{Bmatrix} J_2&J_3&J_4\\ j_q&j_g&j_c\end{Bmatrix}
\begin{Bmatrix} j_h&j_q&J_0\\ j_g&j_p&J_4\end{Bmatrix}
\eta^{J_2}_{ab,hd}\,
\eta^{J_3}_{dp,cq}\,
\Gamma^{J_2}_{gc,ab},\\
Z^{IIIf,J_0}_{pg,qh}
&=-\frac{1}{2}P_{pg}P_{qh}\sum_{abcdJ_2J_3J_4}
\hat J_2^2\hat J_3^2\hat J_4^2
\Big(\bar n_c n_a n_b+\bar n_a\bar n_b n_c\Big) \\
&\qquad\times
\begin{Bmatrix} J_2&J_3&J_4\\ j_q&j_g&j_c\end{Bmatrix}
\begin{Bmatrix} J_2&J_3&J_4\\ j_p&j_h&j_d\end{Bmatrix}
\begin{Bmatrix} j_h&j_q&J_0\\ j_g&j_p&J_4\end{Bmatrix}
\eta^{J_2}_{gc,ab}\,
\eta^{J_3}_{dp,cq}\,
\Gamma^{J_2}_{ab,hd}.
\end{aligned}
$$

## Factorized implementation: `comm223_232`

Source: [src/FactorizedDoubleCommutator.cc](../src/FactorizedDoubleCommutator.cc#L805-L3049).

ArXiv correspondence: this section matches the two-body J-coupled factorized equations in [learn/arxiv.txt](./arxiv.txt), namely
$$
\Gamma^{\mathrm{I}J},\ \Gamma^{\mathrm{II}J},\ \Gamma^{\mathrm{III}_aJ},\ \overline{\overline\Gamma}^{\mathrm{III}_bJ},\ \bar\Gamma^{\mathrm{III}_cJ},\ \Gamma^{\mathrm{IV}_aJ},\ \overline{\overline\Gamma}^{\mathrm{IV}_bJ},\ \bar\Gamma^{\mathrm{IV}_cJ}
$$
in `Jcoupled_Factorized_DoubleCommutator_twobody`, together with the intermediates
$$
\chi^{\epsilon},\ \chi^{\zeta},\ \bar\chi^{\eta},\ \chi^{\theta},\ \bar\chi^{\iota},\ \overline{\overline\chi}^{\kappa},\ \chi^{\lambda}
$$
in `chi_2b_Jcoupled`.

### One-body intermediates for Ia, Ib, IVa, IVb

Source: [src/FactorizedDoubleCommutator.cc](../src/FactorizedDoubleCommutator.cc#L850-L1144).

ArXiv correspondence tag: this block packages the two one-body-intermediate families
$$
\Gamma^{\mathrm{I}J}_{ijkl}
\quad\text{and}\quad
\Gamma^{\mathrm{II}J}_{ijkl},
$$
with intermediates
$$
\chi^{\epsilon}_{ij}
\quad\text{and}\quad
\chi^{\zeta}_{ij}.
$$
In the code-backed form above, the first four insertions are the `\chi^I` / `\chi^{\epsilon}` family and the second four insertions are the `\chi^{II}` / `\chi^{\zeta}` family.

The optimized one-body-intermediate path forms two matrices.  The code loops over holes $i,j$ for the $\bar n_a n_i n_j$ term and over an unrestricted $b$ for the $\bar n_a\bar n_b n_i$ term:

$$
{\color{red}\text{ArXiv correspondence: }\chi^{\epsilon}_{ij}\text{ and }\chi^{\zeta}_{ij}\text{ in }\texttt{chi\_2b\_Jcoupled}.}
$$

$$
\begin{aligned}
\chi^{I}_{pq}
&=\frac{1}{2d_q}\sum_{aijJ_2}\hat J_2^2\bar n_a n_i n_j\,
\eta^{J_2}_{ap,ij}\eta^{J_2}_{ij,aq}
+\frac{1}{2d_q}\sum_{abiJ_2}\hat J_2^2\bar n_a\bar n_b n_i\,
\eta^{J_2}_{ip,ab}\eta^{J_2}_{ab,iq},\\
\chi^{II}_{pq}
&=\frac{1}{2d_q}\sum_{aijJ_2}\hat J_2^2\bar n_a n_i n_j\,
\Gamma^{J_2}_{ap,ij}\eta^{J_2}_{ij,aq}
+\frac{1}{2d_q}\sum_{abiJ_2}\hat J_2^2\bar n_a\bar n_b n_i\,
\Gamma^{J_2}_{ip,ab}\eta^{J_2}_{ab,iq}.
\end{aligned}
$$

The resulting two-body update can be written schematically as the ordinary one-body/two-body insertion on the bra and ket sides:

$$
{\color{red}\text{ArXiv correspondence: }\Gamma^{\mathrm{I}J}_{ijkl}\text{ and }\Gamma^{\mathrm{II}J}_{ijkl}\text{ in }\texttt{Jcoupled\_Factorized\_DoubleCommutator\_twobody}.}
$$

$$
\begin{aligned}
Z^{J}_{pq,rs} &\mathrel{+}=
\sum_b\Big[
\chi^I_{pb}\Gamma^J_{bq,rs}
+\chi^I_{qb}\Gamma^J_{pb,rs}
+\Gamma^J_{pq,bs}\chi^I_{br}
+\Gamma^J_{pq,rb}\chi^I_{bs}
\Big] \\
&\quad
+h_Z\sum_b\Big[
\chi^{II}_{bp}\eta^J_{bq,rs}
+\chi^{II}_{bq}\eta^J_{pb,rs}
-\eta^J_{pq,bs}\chi^{II}_{br}
-\eta^J_{pq,rb}\chi^{II}_{bs}
\Big],
\end{aligned}
$$

with pair-normalization and exchange phases supplied by `GetTBME_norm`.

### Two-body intermediates for diagrams II and III

Source: [src/FactorizedDoubleCommutator.cc](../src/FactorizedDoubleCommutator.cc#L1151-L3049).

ArXiv correspondence tag: this block packages the remaining J-coupled two-body families.  The correspondence is grouped by coupling scheme as follows.

- `\chi^{III}\Gamma + h_Z\Gamma(\chi^{III})^T` corresponds to the ordinary-channel `\chi^{\eta}` family, i.e.
$$
\Gamma^{\mathrm{III}_aJ}_{ijkl}.
$$
- `\bar\chi^{V}_{RC}` and `\chi^V_{final}=\bar\eta\,\bar\chi^V_{RC}` correspond to the double-bar / cross-coupled `\chi^{\eta}` family, i.e.
$$
\overline{\overline\Gamma}^{\mathrm{III}_bJ}_{j\bar l k\bar i}.
$$
- `\bar\chi_\Gamma=\bar\chi^{IV}\bar\Gamma` corresponds to the barred `\chi^{\theta}` family, i.e.
$$
\bar\Gamma^{\mathrm{III}_cJ}_{i\bar l k\bar j}.
$$
- `-\eta\chi^{VI}-\chi^{VI,II}\eta` is the code-packed ordinary-channel realization of the `\chi^{\kappa}` family, i.e.
$$
\Gamma^{\mathrm{IV}_aJ}_{ijkl}.
$$
- `\bar\chi^{VII}` and its inverse-Pandya insertion correspond to the barred/double-barred `\chi^{\iota}` and `\chi^{\lambda}` families, i.e.
$$
\overline{\overline\Gamma}^{\mathrm{IV}_bJ}_{j\bar l k\bar i}
\quad\text{and}\quad
\bar\Gamma^{\mathrm{IV}_cJ}_{i\bar l k\bar j}.
$$

Manual-check note: the code bundles several arXiv subterms into matrix products, so these tags indicate equation families rather than one printed line per C++ line.

The optimized two-body-intermediate path constructs full cross-coupled matrices $\bar\eta$ and $\bar\Gamma$ using the Pandya transform

$$
{\color{red}\text{ArXiv correspondence: common Pandya transform for }\bar\chi^{\eta},\ \bar\chi^{\iota},\ \bar\chi^{\lambda}\text{ families.}}
$$

$$
\bar X^J_{ab,cd}
=\sum_{J'}(-1)^{j_b+j_c+J'}\hat J'^2
\begin{Bmatrix} j_a&j_b&J\\ j_c&j_d&J'\end{Bmatrix}
X^{J'}_{ad,bc} .
$$

The occupation-dressed matrices in the first block are

$$
{\color{red}\text{ArXiv correspondence: occupation factors entering }\chi^{\eta},\ \chi^{\kappa},\ \chi^{\lambda}\text{ families.}}
$$

$$
\begin{aligned}
f_{abc} &= \bar n_a n_b n_c+n_a\bar n_b\bar n_c,\\
f_{abd} &= n_a\bar n_b n_d+\bar n_a n_b\bar n_d,\\
f_{bcd} &= n_b n_c\bar n_d+\bar n_b\bar n_c n_d,\\
f_{acd} &= n_a\bar n_c n_d+\bar n_a n_c\bar n_d .
\end{aligned}
$$

The code then builds the matrix products

$$
{\color{red}\text{ArXiv correspondence: packed intermediates for }\chi^{\eta}\text{ and }\chi^{\kappa}\text{ families.}}
$$

$$
\begin{aligned}
\bar\chi^{III} &= \bar\eta\,\bar\eta_{occ},\\
\bar\chi^{V} &= \bar\Gamma\,\bar\eta_{occ},\\
\bar\chi^{VI} &= \bar\Gamma\,\bar\eta_{occ,d},\\
\bar\chi^{VI,II} &= h_\eta\,\bar\eta_{occ,d}^{T}\bar\Gamma .
\end{aligned}
$$

After inverse/ordinary channel recoupling, the first direct matrix update is

$$
{\color{red}\text{ArXiv correspondence: packed realization of }\Gamma^{\mathrm{III}_aJ}_{ijkl}\text{ and }\Gamma^{\mathrm{IV}_aJ}_{ijkl}.}
$$

$$
Z^{IIa+IIc+IIIc+IIId}
\mathrel{+}=
\chi^{III}\Gamma+h_Z\Gamma(\chi^{III})^T
-\eta\chi^{VI}-\chi^{VI,II}\eta .
$$

The IIb/IId and IIIa/IIIb terms are recoupled through

$$
{\color{red}\text{ArXiv correspondence: recoupling step for }\overline{\overline\Gamma}^{\mathrm{III}_bJ}.}
$$

$$
\bar X^{J}_{ab,cd}
=\sum_{J'}(-1)^{j_b+j_c+J'}\hat J'^2
\begin{Bmatrix} j_a&j_b&J\\ j_c&j_d&J'\end{Bmatrix}
\bar X^{J'}_{ad,bc},
$$

followed by

$$
{\color{red}\text{ArXiv correspondence: matrix-product realization of the }\chi^{\eta}\text{ cross-coupled family.}}
$$

$$
\chi^V_{final}=\bar\eta\,\bar\chi^V_{RC}
$$

and the inverse Pandya transform

$$
{\color{red}\text{ArXiv correspondence: inverse transform returning }\overline{\overline\Gamma}^{\mathrm{III}_bJ}\text{ and related packed families to standard coupling.}}
$$

$$
X^J_{ij,kl}=-(1-P_{ij})(1-P_{kl})(-1)^{J+j_i+j_j}
\sum_{J'}(-1)^{J'+j_i+j_k}\hat J'^2
\begin{Bmatrix} j_j&j_i&J\\ j_k&j_l&J'\end{Bmatrix}
\bar X^{J'}_{jl,ki} .
$$

For IIe/IIf the code builds

$$
{\color{red}\text{ArXiv correspondence: packed intermediate for }\bar\chi^{\theta J}.}
$$

$$
\chi^{IV}=\eta\eta_{occ,c}+(\eta\eta_{occ,d})^T,
$$

and for IIIe/IIIf it builds

$$
{\color{red}\text{ArXiv correspondence: packed intermediate for }\bar\chi^{\iota J}\text{ and }\chi^{\lambda J}\text{ families.}}
$$

$$
\chi^{VII}=\Gamma\eta_{occ,d}+h_\eta\eta_{occ,d}^T\Gamma .
$$

These are recoupled to $\bar\chi^{IV}$ and $\bar\chi^{VII}$, then IIe/IIf use

$$
{\color{red}\text{ArXiv correspondence: matrix-product realization of }\bar\Gamma^{\mathrm{III}_cJ}_{i\bar l k\bar j}.}
$$

$$
\bar\chi_\Gamma=\bar\chi^{IV}\bar\Gamma,
$$

with the same inverse Pandya transform.  The comments in the code write these final two as

$$
{\color{red}\text{ArXiv correspondence: packed inverse-Pandya forms for }\bar\Gamma^{\mathrm{III}_cJ}\text{ and }\bar\Gamma^{\mathrm{IV}_cJ}.}
$$

$$
\begin{aligned}
Z^{IIe,J}_{ij,kl}
&=-\frac12(1-P_{ij})(1-P_{kl})
\sum_{J'}\hat J'^2
\begin{Bmatrix} j_i&j_j&J\\ j_k&j_l&J'\end{Bmatrix}
\bar\chi_\Gamma^{J'}_{il,kj},\\
Z^{IIf,J}_{ij,kl}
&=-\frac12(1-P_{ij})(1-P_{kl})
\sum_{J'}\hat J'^2
\begin{Bmatrix} j_i&j_j&J\\ j_k&j_l&J'\end{Bmatrix}
\bar\chi_{\Gamma,II}^{J'}_{il,kj}.
\end{aligned}
$$

## Reference and factorized implementation: `comm223_132`

Sources: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L13338-L13799), [src/FactorizedDoubleCommutator.cc](../src/FactorizedDoubleCommutator.cc#L3051-L3565).

The `223_132` contribution is the remaining two-body part of

$$
[\eta_2,[\eta_2,\Gamma_2]_{3b}]_{2b}
$$

where one of the contractions is represented by the one-body part of `Eta`.  In the reference code it is implemented by `comm223_132_impl` with three independently callable pieces:

$$
Z^{223\_132}=Z^{\rm ladder}+Z^{\rm cross}+Z^{\rm onebody}.
$$

The wrappers `comm223_132_ladder`, `comm223_132_cross`, and `comm223_132_onebody` select these pieces separately.

ArXiv correspondence tag: `223_132` is not covered by the scalar factorized equations in [learn/arxiv.txt](./arxiv.txt).  Manual checking for this section should instead use the dedicated AMC-derived note and the code-backed equations in this file.

### Common one-body contraction

Source: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L13377-L13403), [src/FactorizedDoubleCommutator.cc](../src/FactorizedDoubleCommutator.cc#L3114-L3131).

All three pieces use the occupation-weighted one-body contraction

$$
f_{ab}\equiv \big(\bar n_a n_b-n_a\bar n_b\big)\eta_{ba},
\qquad j_a=j_b,
$$

or the transposed matrix element $\eta_{ab}$ depending on the orientation of the diagram.  The code enforces the angular-momentum delta by looping over `GetOneBodyChannel(oa.l, oa.j2, oa.tz2)`.

### Ladder piece

Source: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L13357-L13424), [src/FactorizedDoubleCommutator.cc](../src/FactorizedDoubleCommutator.cc#L3058-L3142).

For a two-body channel with $J_0=J_1=J$, the direct ladder contribution is

$$
\begin{aligned}
Z^{\rm ladder,J}_{ij,kl}
&=\sum_{abc} f_{ab}
\Big(
\eta^J_{ca,kl}\Gamma^J_{ij,cb}
-\Gamma^J_{ca,kl}\eta^J_{ij,cb}
\Big) .
\end{aligned}
$$

The C++ loops also include the exchanged storage of the intermediate ket pair $(c,a)$ and apply the usual $1/\sqrt2$ pair normalization for $i=j$ or $k=l$ before storing the TBME.

### Direct one-body-contraction piece

Source: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L13553-L13775), [src/FactorizedDoubleCommutator.cc](../src/FactorizedDoubleCommutator.cc#L3150-L3278).

The reference implementation forms two one-body intermediates

$$
\begin{aligned}
\chi^\eta_{pq}
&=\frac{1}{d_q}\sum_{ijJ_2}\hat J_2^2
\big(\bar n_i n_j-n_i\bar n_j\big)
\eta^{J_2}_{pi,qj}\eta_{ji},\\
\chi^\Gamma_{pq}
&=\frac{1}{d_q}\sum_{ijJ_2}\hat J_2^2
\big(\bar n_i n_j-n_i\bar n_j\big)
\Gamma^{J_2}_{pi,qj}\eta_{ji}.
\end{aligned}
$$

These are then inserted into the external two-body line as

$$
\begin{aligned}
Z^{\rm onebody,J}_{ij,kl}
&=\sum_a\Big(
\chi^\eta_{ia}\Gamma^J_{aj,kl}
+\chi^\eta_{ja}\Gamma^J_{ia,kl}
-\Gamma^J_{ij,al}\chi^\eta_{ak}
-\Gamma^J_{ij,ka}\chi^\eta_{al}
\Big)\\
&\quad-\sum_a\Big(
\chi^\Gamma_{ia}\eta^J_{aj,kl}
+\chi^\Gamma_{ja}\eta^J_{ia,kl}
-\eta^J_{ij,al}\chi^\Gamma_{ak}
-\eta^J_{ij,ka}\chi^\Gamma_{al}
\Big),
\end{aligned}
$$

with the exchange phases supplied explicitly in the reference routine as

$$
(-1)^{J+j_i+j_j},\qquad (-1)^{J+j_k+j_l}
$$

when the insertion is on the exchanged bra or ket leg.  Written with those phases displayed, the eight explicit reference-code terms are

$$
\begin{aligned}
Z^{\rm onebody,J}_{ij,kl}
&=\sum_{abcJ_2}\hat J_2^2 f_{ab}
\Bigg[
\frac{1}{d_i}\eta^{J_2}_{ia,cb}\Gamma^J_{cj,kl}
-\frac{(-1)^{J+j_i+j_j}}{d_j}\eta^{J_2}_{ja,cb}\Gamma^J_{ci,kl}\\
&\qquad
-\frac{1}{d_k}\Gamma^J_{ij,cl}\eta^{J_2}_{ac,bk}
+\frac{(-1)^{J+j_k+j_l}}{d_l}\Gamma^J_{ij,ck}\eta^{J_2}_{ac,bl}\\
&\qquad
+\frac{1}{d_l}\eta^J_{ij,kc}\Gamma^{J_2}_{ca,lb}
-\frac{(-1)^{J+j_k+j_l}}{d_k}\eta^J_{ij,lc}\Gamma^{J_2}_{ca,kb}\\
&\qquad
-\frac{1}{d_j}\Gamma^{J_2}_{ja,cb}\eta^J_{ic,kl}
+\frac{(-1)^{J+j_i+j_j}}{d_i}\Gamma^{J_2}_{ia,cb}\eta^J_{jc,kl}
\Bigg].
\end{aligned}
$$

### Cross-coupled reference piece

Source: [src/ReferenceImplementations.cc](../src/ReferenceImplementations.cc#L13428-L13549).

The direct reference implementation writes the cross term in the normal-coupled basis with a nine-j recoupling.  A representative direct term is

$$
\begin{aligned}
Z^{{\rm cross},J_0}_{ij,kl}
&\supset \sum_{abcJ_2J_3}
(-1)^{J_2+J_3+j_a+j_c+J_0}
\hat J_2^2\hat J_3^2
\begin{Bmatrix}
j_i&j_c&J_2\\
j_j&J_3&j_a\\
J_0&j_l&j_k
\end{Bmatrix}
\eta^{J_2}_{ic,ka}\,\eta_{ab}\,\Gamma^{J_3}_{bj,cl}
\big(\bar n_a n_b-n_a\bar n_b\big).
\end{aligned}
$$

The code expands this seed into the eight antisymmetrized variants generated by exchanging $i\leftrightarrow j$ and $k\leftrightarrow l$, with signs and extra phases

$$
1,
\quad (-1)^{j_i+j_j},
\quad (-1)^{j_k+j_l},
\quad (-1)^{J_0+j_i+j_j+j_k+j_l}.
$$

The alternate direct term in the lambda is the same recoupling with the one-body contraction placed on the other side of the two-body matrix element,

$$
-\eta^{J_2}_{ib,kc}\,\eta_{ab}\,\Gamma^{J_3}_{cj,al},
$$

and its exchanged partners.

### Factorized cross implementation

Source: [src/FactorizedDoubleCommutator.cc](../src/FactorizedDoubleCommutator.cc#L3287-L3565).

The optimized `FactorizedDoubleCommutator::comm223_132_cross` evaluates the same cross piece through cross-coupled matrix products.  It first forms

$$
\bar X^J_{ab,cd}
=\sum_{J'}(-1)^{j_b+j_c+J'}\hat J'^2
\begin{Bmatrix} j_a&j_b&J\\ j_c&j_d&J'\end{Bmatrix}
X^{J'}_{ad,bc},
\qquad X\in\{\eta,\Gamma\}.
$$

The one-body-contraction matrix in the same cross-coupled space is

$$
\bar\eta^{(1),J}_{ab,cd}
=\delta_{ac}\,\eta_{db}\big(\bar n_b n_d-n_b\bar n_d\big)
+\delta_{bd}\,\eta_{ac}\big(\bar n_a n_c-n_a\bar n_c\big).
$$

The cross intermediate is the simple matrix product

$$
\bar\chi^{132,J}=\bar\eta^J\,\bar\eta^{(1),J}\,\bar\Gamma^J .
$$

It is transformed back to the normal-coupled TBME through the same inverse Pandya antisymmetrizer used elsewhere in the factorized double commutator:

$$
X^J_{ij,kl}=-(1-P_{ij})(1-P_{kl})(-1)^{J+j_i+j_j}
\sum_{J'}(-1)^{J'+j_i+j_k}\hat J'^2
\begin{Bmatrix} j_j&j_i&J\\ j_k&j_l&J'\end{Bmatrix}
\bar X^{J'}_{jl,ki} .
$$

In component form the code accumulates four inverse-Pandya pieces,

$$
Z^{\rm cross,J}_{ij,kl}
=C_{ij,kl}
-(-1)^{j_i+j_j-J}C_{ji,kl}
-(-1)^{j_k+j_l-J}C_{ij,lk}
+(-1)^{j_i+j_j+j_k+j_l}C_{ji,lk},
$$

where each $C$ is a $J'$ sum over the appropriate element of $\bar\chi^{132,J'}$.