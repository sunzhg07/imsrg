# Analysis of `comm223ss` in `ReferenceImplementations.cc`

## Function signature and expression

```
ZJ1J2J3_ijklmn = sum_a  PJ1J3(ij/k) PJ1J2(lm/n) sqrt((2J1+1)(2J2+1))
                         { n  a  J1 }
                         { k  J3 J2 }
                         * ( XJ1_ijna YJ2_kalm - YJ1_ijna XJ2_kalm )
```

Located at `src/ReferenceImplementations.cc` around line 1447.

---

## Code structure

- Outer loops: `ch3bra`, `ch3ket`, `ibra`, `iket` over three-body channels and ket pairs.
- `bra = (p=i, q=j, r=k, Jpq=J1)`, `ket = (p=l, q=m, r=n, Jpq=J2)`, both with total `twoJ`.
- Three permutation cases: `{ABC, CBA, ACB}` (the P(ij/k) antisymmetrizer).
  - `ABC`: identity, no J1' sum, rec = +1
  - `CBA`: odd permutation, J1' summed via 6j
  - `ACB`: odd permutation, J1' summed via 6j
- The commutator 6j symbol (`sixj` in code) is:
  ```
  GetCachedSixJ(o3.j2, twoJ, J1p, o6.j2, j2a, J2p)
  = { j_k   J    J1' }
    { j_I6  j_a  J2' }
  ```
  where `I6` comes from the **ket** permutation (the third index after permuting lmn) and `I3` from the **bra** permutation (third index after permuting ijk).

- Accumulation line:
  ```cpp
  zijklmn += rec_ijk * rec_lmn * sixj
           * sqrt((2*J1p+1)*(2*J2p+1))
           * (x_126a * y_3a45 - y_126a * x_3a45);
  ```
  where `x_126a = X2.GetTBME_J(J1p, J1p, I1, I2, I6, a)` and
        `x_3a45 = X2.GetTBME_J(J2p, J2p, I3, a, I4, I5)`.

---

## RecouplingCoefficient and PermutationPhase

From `src/ThreeBodyStorage.cc`:

```cpp
RecouplingCoefficient(perm, ja, jb, jc, Jab_in, Jab, twoJ):
  ABC -> delta(Jab_in, Jab)
  CBA -> -sqrt((2Jab_in+1)(2Jab+1)) * SixJ(ja, jb, Jab, jc, J, Jab_in)
  ACB -> (-1)^(jb+jc+Jab_in-Jab) * sqrt((2Jab_in+1)(2Jab+1)) * SixJ(jb, ja, Jab, jc, J, Jab_in)

PermutationPhase:
  ABC, BCA, CAB -> +1
  BAC, CBA, ACB -> -1
```

`P(ij/k) = PermutationPhase * RecouplingCoefficient`.

---

## Case study: bra = |11(J1)2; J⟩  (i=1, j=1, k=2)

Here j_i = j_j (same orbit), j_k is orbit 2's j.

### Term 1 — ABC: (I1,I2,I3) = (1,1,2)

- J1' = J1 (fixed, no sum)
- `rec_ABC = +1`
- 2B bra ME: `X^{J1}(1,1 | I6, a)`

### Term 2 — CBA: (I1,I2,I3) = (2,1,1)

- J1' summed over allowed range
- `RecouplingCoefficient(CBA, j_i, j_i, j_k, J1', J1, 2J)`:

$$-\sqrt{(2J_1'+1)(2J_1+1)}\,\begin{Bmatrix}j_i & j_i & J_1\\j_k & J & J_1'\end{Bmatrix}$$

- `PermutationPhase(CBA) = -1`, so:

$$\text{rec}_{CBA} = +\sqrt{(2J_1'+1)(2J_1+1)}\,\begin{Bmatrix}j_i & j_i & J_1\\j_k & J & J_1'\end{Bmatrix}$$

- 2B bra ME: `X^{J1'}(2,1 | I6, a)`

### Term 3 — ACB: (I1,I2,I3) = (1,2,1)

- J1' summed over allowed range
- `RecouplingCoefficient(ACB, j_i, j_i, j_k, J1', J1, 2J)`:

$$(-1)^{j_i+j_k+J_1'-J_1}\,\sqrt{(2J_1'+1)(2J_1+1)}\,\begin{Bmatrix}j_i & j_i & J_1\\j_k & J & J_1'\end{Bmatrix}$$

- `PermutationPhase(ACB) = -1`, so:

$$\text{rec}_{ACB} = (-1)^{j_i+j_k+J_1'-J_1+1}\,\sqrt{(2J_1'+1)(2J_1+1)}\,\begin{Bmatrix}j_i & j_i & J_1\\j_k & J & J_1'\end{Bmatrix}$$

- 2B bra ME: `X^{J1'}(1,2 | I6, a)`

### Full bra permutation operator

$$P^{J_1 J}(ij/k)\Big|_{i=j=1,\,k=2}
= X^{J_1}_{1\,1,\,I_6 a}
+ \sum_{J_1'}\sqrt{(2J_1'+1)(2J_1+1)}\,\begin{Bmatrix}j_i & j_i & J_1\\j_k & J & J_1'\end{Bmatrix}
\left[X^{J_1'}_{2\,1,\,I_6 a}
+ (-1)^{j_i+j_k+J_1'-J_1+1}\,X^{J_1'}_{1\,2,\,I_6 a}\right]$$

Note: CBA and ACB yield the **same 6j** (since j_i = j_j), differing only by the phase from ACB.

---

## Case study: ket = |11(J2)2; J⟩  (l=1, m=1, n=2)

Same orbital configuration, same structure, but with J2, and the 2B ME is `Y^{J2'}(I3, a | I4, I5)`.

### Term 1 — ABC: (I4,I5,I6) = (1,1,2)

- J2' = J2 (fixed)
- `rec_ABC = +1`
- Commutator 6j: `{ j_k  J  J1' ; j_k  j_a  J2 }`
- 2B ket ME: `Y^{J2}(I3, a | 1, 1)`

### Term 2 — CBA: (I4,I5,I6) = (2,1,1)

- J2' summed

$$\text{rec}_{CBA} = +\sqrt{(2J_2'+1)(2J_2+1)}\,\begin{Bmatrix}j_i & j_i & J_2\\j_k & J & J_2'\end{Bmatrix}$$

- Commutator 6j: `{ j_k  J  J1' ; j_i  j_a  J2' }`
- 2B ket ME: `Y^{J2'}(I3, a | 2, 1)`

### Term 3 — ACB: (I4,I5,I6) = (1,2,1)

- J2' summed

$$\text{rec}_{ACB} = (-1)^{j_i+j_k+J_2'-J_2+1}\,\sqrt{(2J_2'+1)(2J_2+1)}\,\begin{Bmatrix}j_i & j_i & J_2\\j_k & J & J_2'\end{Bmatrix}$$

- Commutator 6j: `{ j_k  J  J1' ; j_i  j_a  J2' }`
- 2B ket ME: `Y^{J2'}(I3, a | 1, 2)`

### Full ket permutation operator

$$P^{J_2 J}(lm/n)\Big|_{l=m=1,\,n=2}
= Y^{J_2}_{I_3 a,\,1\,1}
+ \sum_{J_2'}\sqrt{(2J_2'+1)(2J_2+1)}\,\begin{Bmatrix}j_i & j_i & J_2\\j_k & J & J_2'\end{Bmatrix}
\Bigl[Y^{J_2'}_{I_3 a,\,2\,1}
+(-1)^{j_i+j_k+J_2'-J_2+1}\,Y^{J_2'}_{I_3 a,\,1\,2}\Bigr]$$

---

## Symmetry observation

The bra and ket permutation operators are **structurally identical**:
- Same 6j symbol form (because i=j=l=m with j_i=j_j=j_l=j_m)
- Same phase pattern for ABC / CBA / ACB terms
- Differ only in: which coupling quantum number (J1 vs J2), and which two-body ME (bra indices 1,2 vs ket indices I3,a)

---

## Full matrix element (combining bra+ket+commutator 6j)

$$Z^{J_1 J_2 J}_{11\,2,\,11\,2} = \sum_{J_1' J_2' a}
\;\text{rec}_{ijk}(J_1')\;\text{rec}_{lmn}(J_2')
\;\sqrt{(2J_1'+1)(2J_2'+1)}
\begin{Bmatrix}j_k & J & J_1'\\j_{I_6} & j_a & J_2'\end{Bmatrix}
\bigl(X^{J_1'}_{I_1 I_2,\,I_6 a}\;Y^{J_2'}_{I_3 a,\,I_4 I_5} - X\leftrightarrow Y\bigr)$$

where the three (I1,I2,I3) and (I4,I5,I6) index sets come from the three permutation terms above,
and the commutator 6j has argument `(o3.j2, twoJ, J1p, o6.j2, j2a, J2p)` in the code,
i.e. `{ j_{I3}  J  J1' ; j_{I6}  j_a  J2' }`.
