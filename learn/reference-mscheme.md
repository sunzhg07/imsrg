# M-scheme form of `reference_eq.md`

This file transcribes the documented J-coupled equations in `learn/reference_eq.md` into m-scheme form.
All explicit angular-momentum recoupling factors, $(2J+1)$ weights, and sums over $J$ are removed.
The two-body matrix elements are written in the uncoupled m-scheme basis as $\eta_{ab,cd}$, $\Gamma_{ab,cd}$, and $\gamma_{ab,cd}$.

## comm223_231

### Diagram I

$$
Z^{I}_{pq}
= \frac{1}{2} \sum_{abcde}
\Bigl(
\bar n_a \bar n_c n_b n_d
- \bar n_b \bar n_d n_a n_c
- \bar n_b \bar n_e n_a n_c
+ \bar n_a \bar n_c n_b n_e
\Bigr) \\
\times
\left(
\eta_{bd,ac} \, \eta_{ac,be}
+ \eta_{bd,ca} \, \eta_{ca,be}
\right)
\Gamma_{ep,dq}
$$

### Diagram IIa

$$
Z^{IIa}_{pq}
= \sum_{abcde}
\Bigl(
\bar n_a \bar n_c n_b n_d
- \bar n_b \bar n_d n_a n_c
\Bigr)
\eta_{bd,ac} \, \eta_{cp,de} \, \Gamma_{ae,bq}
$$

### Diagram IIb

$$
Z^{IIb}_{pq}
= \frac{1}{4} \sum_{abcde}
\Bigl(
\bar n_a \bar n_d n_b n_e
- \bar n_b \bar n_e n_a n_d
\Bigr)
\eta_{be,ad} \, \eta_{cp,be} \, \Gamma_{ad,cq}
$$

### Diagram IIc

$$
Z^{IIc}_{pq}
= \sum_{abcde}
\Bigl(
\bar n_a \bar n_d n_b n_e
- \bar n_b \bar n_e n_a n_d
\Bigr)
\eta_{be,da} \, \eta_{ca,bq} \, \Gamma_{dp,ce}
$$

### Diagram IId

$$
Z^{IId}_{pq}
= -\frac{1}{4} \sum_{abcde}
\Bigl(
\bar n_c \bar n_d n_a n_e
- \bar n_a \bar n_e n_c n_d
\Bigr)
\eta_{ae,cd} \, \eta_{cd,bq} \, \Gamma_{bp,ae}
$$

### Diagram IIIa

$$
Z^{IIIa}_{pq}
= \frac{1}{2} \sum_{abcde}
\Bigl(
\bar n_a \bar n_e n_b n_c
- \bar n_b \bar n_c n_a n_e
\Bigr)
\eta_{bc,ae} \, \eta_{ep,dq} \, \Gamma_{ad,bc}
$$

### Diagram IIIb

$$
Z^{IIIb}_{pq}
= -\frac{1}{2} \sum_{abcde}
\Bigl(
\bar n_a \bar n_c n_b n_d
- \bar n_b \bar n_d n_a n_c
\Bigr)
\eta_{bd,ac} \, \eta_{ep,dq} \, \Gamma_{ac,be}
$$

## comm223_232

### CHI_I and CHI_II

$$
\Chi_{I,pq}
= \frac{1}{2} \sum_{a i j}
\Bigl(
\bar n_a \bar n_c n_b
- \bar n_b n_a n_c
\Bigr)
\eta_{b p, a c} \, \eta_{a c, b q}
$$

$$
\Chi_{II,pq}
= \frac{1}{2} \sum_{a i j}
\Bigl(
\bar n_b \bar n_c n_a
- \bar n_a n_b n_c
\Bigr)
\eta_{b c, a q} \, \gamma_{a p, b c}
$$

### Diagram IIa

$$
Z^{IIa}_{p g, q h}
= \sum_{d a}
\Chi_{III,p g, d a} \, \Gamma_{d a, q h}
$$

### Diagram IIc

$$
Z^{IIc}_{p g, q h}
= \sum_{d a}
\Gamma_{p g, a d} \, \Chi_{III,a d, q h}
$$

### Diagram IIb

$$
Z^{IIb}_{p g, q h}
= -P_{p g} P_{q h}
\sum_{a d}
(\bar\Chi_{III})^T_{d q a p} \, \Gamma_{g a, h d}
$$

### Diagram IId

$$
Z^{IId}_{p g, q h}
= -P_{p g} P_{q h}
\sum_{a b c d}
(-1)^{\phi}
(\bar\Chi_{III})^T_{d g h a} \, \Gamma_{d p, a q}
$$

In these m-scheme forms, $P_{p g}$ and $P_{q h}$ denote label exchange operators, and $(-1)^{\phi}$ is the phase factor inherited from the J-coupled IId recoupling.
