# Reference equations from `ReferenceImplementations.cc`

This file contains only the equations mapped to the code in `src/ReferenceImplementations.cc` for `comm223_231` and `comm223_232`.

- `comm223_231` begins at line **10309** in `src/ReferenceImplementations.cc`
- `comm223_232` begins at line **11067** in `src/ReferenceImplementations.cc`

## comm223_231

### Diagram I (around line 10322)

- Intermediate one-body operator:

$$
\Chi_{221_a} = \sum_{J_0} \hat J_0 \, (\bar n_a \bar n_c n_b n_d - \bar n_b \bar n_d n_a n_c - \bar n_b \bar n_e n_a n_c + \bar n_a \bar n_c n_b n_e) \, \eta^{J_0}_{bd,ac} \, \eta^{J_0}_{ac,be}
$$

- Final diagram I contribution in code:

$$
I_{pq} = \frac{1}{2 (2 j_p + 1)} \sum_{de J_1} \Chi_{221_a}(d,e) \, (2 J_1 + 1) \, \Gamma^{J_1}_{e p, d q}
$$

### Diagram IIa (Pandya transform)

- Pandya transform used in code:

$$
\bar X^J_{i j' k l'} = - \sum_{J'} (2 J' + 1) 
\begin{Bmatrix} i & j & J \\ k & l & J' \end{Bmatrix} X^{J'}_{i l k j}
$$

### Chi_222_a (around line 10800)

- Intermediate two-body operator in code:

$$
\Chi_{222_a} = \sum_{c d e f} (\bar n_e \bar n_d n_f n_c - \bar n_f \bar n_c n_e n_d) \, \bar\eta_{p e d c} \, \bar\eta_{c d a b}
$$

### Diagram IIa and IIc

- Diagram IIa contribution in code:

$$
IIa_{pq} = \frac{1}{2 j_p + 1} \sum_{a b e J_3} \Chi_{222_a,p e a b} \, \Gamma^{J_3}_{a b q e}
$$

- Diagram IIc contribution in code:

$$
IIc_{pq} = -\frac{1}{2 j_p + 1} \sum_{a b e J_3} \Chi_{222_a,e q a b} \, \Gamma^{J_3}_{a b e p}
$$

### Diagram IIb

- Intermediate two-body operator in code:

$$
\Chi_{222_b} = \sum_{b e} (\bar n_a \bar n_d n_b n_e - \bar n_b \bar n_e n_a n_d) \, \eta^{J_0}_{c p, b e} \, \eta^{J_0}_{b e, a d}
$$

- Diagram IIb contribution in code:

$$
IIb_{pq} = \frac{1}{4 (2 j_p + 1)} \sum_{a c d J_0} \Chi_{222_b,c p a d} \, \Gamma^{J_0}_{a d, c q}
$$

### Diagram IId

- Diagram IId contribution in code:

$$
IId_{pq} = -\frac{1}{4 (2 j_p + 1)} \sum_{a b e J_0} \Chi_{222_b,b q a e} \, \Gamma^{J_0}_{b p, a e}
$$

### Diagram IIIa / IIIb

- Intermediate one-body operator in code:

$$
\Chi_{221_b} = \sum_{J_0} \hat J_0 \, (\bar n_a \bar n_e n_b n_c - \bar n_b \bar n_c n_a n_e) \, \eta^{J_0}_{b c, a e} \, \Gamma^{J_0}_{a d, b c}
$$

- Final III contributions in code:

$$
Z_{pq} = \frac{1}{2 (2 j_p + 1)} \sum_{d e J_1} (2J_1+1) \, \Chi_{221_b}(d,e) \, \eta^{J_1}_{e p, d q} 
- \frac{1}{2 (2 j_p + 1)} h_Z \sum_{d e J_1} (2J_1+1) \, \Chi_{221_b}(e,d) \, \eta^{J_1}_{e p, d q}
$$

## comm223_232

### CHI_I and CHI_II (around line 11090)

- Code comment for the intermediate one-body operators:

$$
\Chi_{I,pq} = \frac{1}{2} \sum_{a i j J_2} \hat J_2 \, (\bar n_a \bar n_c n_b - \bar n_b n_a n_c) \, \eta^{J_2}_{b p a c} \, \eta^{J_2}_{a c b q}
$$

$$
\Chi_{II,pq} = \frac{1}{2} \sum_{a i j J_2} \hat J_2 \, (\bar n_b \bar n_c n_a - \bar n_a n_b n_c) \, \eta^{J_2}_{b c a q} \, \gamma^{J_2}_{a p b c}
$$

### Diagram I / IV structure

- The code computes scalar or non-scalar two-body contributions using these intermediates and then adds to `Z.TwoBody`.

### Factorization of IIa and IIc (around line 11480)

- Intermediate cross-coupled operator in code:

$$
\bar\Chi_{III} = \sum_{b c J_2 J_3} (\bar n_b \bar n_d n_c - \bar n_c n_b n_d) \, \bar\eta_{p a c b} \, \bar\eta_{c b d g}
$$

- Inverse Pandya transformation used in code:

$$
X^J_{i j k l} = - (1 - P_{ij}) \sum_{J'} (2J' + 1) 
\begin{Bmatrix} i & j & J \\ k & l & J' \end{Bmatrix} \bar X^{J'}_{i l' k j'}
$$

- IIa and IIc in code:

$$
IIa_{p g, q h} = \sum_{d a} \Chi_{III}^{J_0}{}_{p g, d a} \Gamma^{J_0}_{d a, q h}
$$

$$
IIc_{p g, q h} = \sum_{d a} \Gamma^{J_0}_{p g, a d} \Chi_{III}^{J_0}{}_{a d, q h}
$$

### Factorization of IIb and IId (later in comm223_232)

- Code comment for the two-body parts:

$$
IIb^{J_0}_{p g q h} = - P_{pg} P_{qh} \sum_{a d J_4 J_5} (2 J_4 + 1) (2 J_5 + 1) 
\begin{Bmatrix} J_0 & J_4 & J_5 \\ j_d & j_q & j_h \end{Bmatrix} 
\begin{Bmatrix} J_0 & J_4 & J_5 \\ j_a & j_p & j_g \end{Bmatrix} \, (\bar{\Chi}_{III}^{J_5})^T_{d q a p} \Gamma^{J_4}_{g a h d}
$$

$$
IId^{J_0}_{p g q h} = - P_{pg} P_{qh} \sum_{a b c d J_4 J_5} (2 J_4 + 1) (2 J_5 + 1) 
\begin{Bmatrix} J_0 & J_5 & J_4 \\ j_d & j_p & j_g \end{Bmatrix} 
\begin{Bmatrix} J_5 & J_4 & J_0 \\ j_q & j_h & j_a \end{Bmatrix} \, (-1)^{j_d+j_g+j_h+j_a} (\bar{\Chi}_{III}^{J_5})^T_{d g h a} \Gamma^{J_4}_{d p a q}
$$
