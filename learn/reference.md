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
