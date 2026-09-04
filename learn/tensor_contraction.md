# T×T → scalar: gold rule for future coding

**Lock \(Z\) to m-scheme, not to AMC’s printed hats.**
AMC is the **coupled product** \([A^{(\lambda)}\times B^{(\lambda)}]^{(0)}\) (\(\hat\lambda^{-1}\)).
Code \(Z\) is the **m-trace / leftover \(m\)-average** (\(\hat\lambda^{-2}\)). Convert hats; do not retune AMC input.

Also read `learn/amc_tts/REDUCED_UNREDUCED.md` before any m ↔ J numeric compare.

**Where the code lives** (same split as ss / st):

| | file |
|---|---|
| m-scheme gold | `UnitTest.cc` (`Mscheme_Test_commXXXtts`) — calls **production** `Commutator::` |
| production tts | `src/TensorCommutators.cc` |
| naive J-scheme / tested-not-production | `src/ReferenceImplementations.cc` (`Test_tts_against_ref`; `comm232tts_bare`; factorized 223_231 / 223_232) |

---

## Gold rule (do this every time)

Three different angular objects. Name which one you are coding.

| object | definition | hat | where |
|---|---|---|---|
| **coupled product** | \([A^{(\lambda)}\times B^{(\lambda)}]^{(0)}_0=\hat\lambda^{-1}\sum_M(-1)^{\lambda-M}A_M B_{-M}\) | \(\hat\lambda^{-1}\) | AMC `.tex` (`reduce=false` on scalar \(Z\)) |
| **dot product** (Edmonds) | \(A\cdot B=\sum_M(-1)^M A_M B_{-M}=\hat\lambda\,[A\times B]^{(0)}\) | — | textbooks / EOM overlap (not the commutator) |
| **m-trace / leftover average (gold)** | same \(m\)-labels, **no** extra \((-1)^{\lambda-M}\) | \(\hat\lambda^{-2}\) after WE | `GetMscheme` + `Mscheme_Test_commXXXtts` |

Leftover of the scalar \(Z\) (this is global, not per 6j):

| leftover \(Z\) | gold |
|---|---|
| 0-body | sum all \(m\) of every contracted index |
| scalar 1b | \(Z_{ij}=\dfrac{1}{2j_i+1}\sum_m\) Wick (AMC \(\hat\jmath_i^{-2}\)) |
| scalar 2b | \(\tilde Z^{J_0}=\dfrac{1}{2J_0+1}\sum_m\mathrm{CG}_{ij}\,\mathrm{CG}_{kl}\,Z(m)\) vs tilde `GetTBME_J` (AMC \(\hat J_0^{-2}\)) |
| scalar 3b | \(Z^{J_{ab} j,\,J_{de} j}=\dfrac{1}{2j+1}\sum_m\mathrm{CG}_{ij}\,\mathrm{CG}_{ijk}\,\mathrm{CG}_{lm}\,\mathrm{CG}_{lmn}\,Z(m)\) vs `GetME_pn` |

`GetMschemeMatrixElement_*` is always the **physical** ME of one component \(T^\lambda_M\). Same \(m\)-labels already force \(M_Y=-M_X\). Do **not** put the coupling CG or \((-1)^{\lambda-M}\) on that loop.

**Code \(Z\) = leftover average = \(\hat\lambda^{-2}\).** Keep AMC’s leftover \(\hat\jmath_i^{-2}\) / \(\hat J_0^{-2}\). Convert only \(\hat\lambda^{-1}\to\hat\lambda^{-2}\).

If a ratio clusters at \(\hat\lambda^{\pm 1}\) or \(1/\sqrt{2J+1}\), **stop** — packaging, not a 6j.

**Dummy in/out.** A contracted index must leave one operator and enter the other. At \(\lambda=0\), \(\Omega_{ijkl}=h_\Omega\Omega_{klij}\) hides illegal slots; at \(\lambda\neq 0\) that identity is false for `GetMscheme`. Rewrite slots + leftover CG; do not substitute swapped MEs. Occupation-weighted \(\chi\) is often NH. Full note: `learn/tensor_contraction` (closed vs open leftover, and the hermiticity section).

---

## Wick string vs m-scheme **code** (ss / st / tts)

The Wick algebra is the same string as `ss`. The m-scheme **implementation** is not. Write a separate `Mscheme_Test_commXXXtts`. Do **not** reuse the `ss`/`st` gold loop.

| | \(X,Y\) | \(Z\) | gold on leftover |
|---|---|---|---|
| **ss** | unreduced scalars | unreduced scalar | stretched leftover \(m\) is enough |
| **st** | scalar + rank \(\lambda\) | **tensor** | one \(\mu\); `GetMscheme(Z)` |
| **tts** | two rank-\(\lambda\) tensors | **scalar** unreduced | project leftover (table above). Stretched \(m\) is **not** a scalar |

`comm111tts` first FAIL: same \(\sum_a(X_{ia}Y_{aj}-Y_{ia}X_{aj})\) as 111ss, but at stretched \(m\) the product is ranks \(0\ldots 2\lambda\). \(j=1/2\) fake-passed. The scalar is the \(m\)-average.

`comm222_pp_hhtts` first FAIL: same ss Wick compared at one leftover \(m\)-tuple to `GetMscheme(Z)`. A \(J=0\) ket fake-passed as a global minus; mixed \(J\) ratios were not a constant. Fix = 111 leftover with \(j_i\to J_0\), 220 tilde store.

Contracted indices: sum **all** their \(m\) (111: all \(m_a\); ladder: all \(m_a,m_b\)). Do not force \(M_{\mathrm{ab}}=M_{ij}\) (that is ss \(\mu=0\)).

Do **not** use 132tts “Wick at each leftover \(m\)” as the 2b\(\times\)2b template. 132 is 1b\(\times\)3b after the 1b \(\mu\) is already contracted.

Pair \(J_0,J_2,\lambda\) are **integers**. Orbit \(j=1/2\) is not a two-body \(J\).

---

## Try → FAIL → what actually worked

Name leftover gold first. Convert AMC hats to match it. Almost every FAIL was packaging, a missing Wick factor, or copying the wrong ss coupling — not a wrong 6j.

Print `Zm/ZJ` on the first ~20 fails, plus \(\hat\lambda\), \(1/\hat\lambda\), \(\hat J\), \(2J+1\). Do not retune 6js until the cluster is gone.

| cluster | meaning |
|---|---|
| all \(\hat\lambda\) or \(1/\hat\lambda\) | leftover vs coupled-product hats |
| all \(1/\sqrt{2J+1}\) or \(2J+1\) | leftover \(J\) packaging / reduced vs unreduced |
| all \(-1\) | overall sign or Eq3/4 commutator minus |
| mixed, not a hat | wrong perm, wrong ME slots, or ss contraction on a tts 6j |

**1. Reusing the ss m-scheme test.** FAIL: 111tts compared ss Wick at one leftover \(m\) to `GetMscheme(Z)`. \(j=1/2\) passed; larger \(j\) did not. Worked: leftover average \(Z_{ij}=\frac{1}{2j_i+1}\sum_m\) Wick. Same for leftover 2b/3b with the CG project vs tilde `GetTBME_J` / `GetME_pn`.

**2. Copying one leftover \(m\) onto a 2b ladder.** FAIL: 222 pp/hh, one \(m\)-tuple vs `GetMscheme`. A \(J=0\) ket looked like a global minus; mixed \(J\) ratios were not constant. Worked: leftover is 111 with \(j_i\to J_0\), store is 220 tilde.

**3. Forcing contracted \(M\) like ss \(\mu=0\).** FAIL: requiring \(M_{ab}=M_{ij}\). Worked: sum **all** \(m\) of contracted legs. Leftover \(m\) is averaged separately.

**4. AMC omits combinatorial factors.** FAIL: 220/221/222/231/232 missing \(\tfrac14\) or \(\tfrac12\), or missing \((1-P)\). Worked: restore from the ss Wick, not from the tex. AMC is the primitive product.

**5. Phase: sometimes drop \(\lambda\), usually keep it.** FAIL: 221 needed drop \(\lambda\) from \((-1)^{J_0+J_1+\lambda}\). Copying that drop onto 111/220/222 broke even \(\lambda\) (the extra minus is invisible there). Worked: drop \(\lambda\) **only** on 221.

**6. Leftover \(\hat J_0\) is topology-dependent.** FAIL: assuming leftover 2b is always \(\hat J_0^{-2}\) like 222. Worked: 222/132 keep \(\hat J_0^{-2}\). **232 and 122 AMC leftover is \(\hat J_0^{-1}\)**. Convert only \(\hat\lambda\). If you force \(\hat J_0^{-2}\) on those, ratios sit at \(\hat J_0\).

**7. \((1-P)\) is coupled, not a bare minus.** FAIL: 232 used `psign=-1`. That zeros \(i=j\) and misses even/odd \(J_0\). Worked: \((1-P_{ij}^{J_0})\) with \(P=(-1)^{j_i+j_j-J_0}\). Eq3/4 = \(-(X\leftrightarrow Y)\) of Eq1/2, not the same overall minus on all four AMC products.

**8. Native antiherm 3b vs leftover 2b gold.** FAIL: 132/232 leftover 2b vs `GetMscheme` with native antiherm 3b (holes like \(|abb\rangle\)). Worked: **native herm 3b** on the 3b tensor. Native antiherm **2b** is fine when leftover is 3b (223).

**9. Do not steal the ss J-scheme contraction.** FAIL: 223tts tried ss slots \(X_{I_1 I_2 I_6 a} Y_{I_3 a I_4 I_5}\) plus a 9j “generalization.” \(\lambda=0\) did not reproduce leftover gold. ss `{ABC,CBA,ACB}` on the ket gave mixed ratios. Worked: **strict AMC primitive** \(X_{ijla}Y_{akmn}\) with AMC’s 6j×9j. Wick \((1-P_{ik}-P_{jk})(1-P_{lm}-P_{ln})\): bra ABC/CBA/ACB; ket ABC/CBA/**BAC**. ss uses ACB on the ket because its primitive is \(X_{ijna}\), not \(X_{ijla}\). Overall minus vs leftover (AMC Eq1 is \(+XY\)). ABC-only MEs were exact \(-1\) (primitive sign); mixed ratios were the wrong ket perm.

**10. \(\lambda=0\): same kernel, not a ss call.** FAIL: `if (not reduced) commXXss` is two codes. Worked (223 prototype): one kernel for all \(\lambda\). If \(X\) or \(Y\) is an unreduced scalar, `MakeReduced` then the same kernel. \(Z\) stays unreduced leftover. Older tts still fall through to `ss`.

**Do not copy.** 132 leftover-\(m\) Wick is not the 2b×2b template. 221’s drop-\(\lambda\) phase is not global. 232’s \(\hat J_0^{-1}\) is not 222’s \(\hat J_0^{-2}\). ss `{ABC,CBA,ACB}` is not a universal 3-body perm set. `GetTBME_J(J_{\mathrm{bra}},J_{\mathrm{ket}},\ldots)\): hats sit on the **pairs**.

Lock: He4, seed 17, \(e_{\max}=1\), \(\lambda=0,1,2\), `OMP_NUM_THREADS=1`. Odd parity when you need \(n_a-n_b\) hp (132).

---

## Locked J-scheme formulas

MEs: 1b = `OneBody` (Edmonds reduced if `IsReduced()`). 2b = **tilde** `GetTBME_J`. Occupations as written. Free sums unless noted. Unreduced \(\lambda=0\): scalar `ss` (no extra hats).

`GetTBME_J(j_{\mathrm{bra}},j_{\mathrm{ket}},a,b,c,d)`: \((ab)\) couples to \(j_{\mathrm{bra}}\). Hats sit on the **pairs**, not on loop names in the wrong order.

### `comm110tts` (locked \(\lambda=0,1,2\))

m-gold: \(Z_0=\sum_{ia,m}n_i\bar n_a(X_{ia}Y_{ai}-Y_{ia}X_{ai})\). \(X_{ia}\) and \(X_{ai}\) are different MEs.

Reduced: \(Z_0=\sum_{ia}n_i\bar n_a(-1)^{j_i+j_a+\lambda}\hat\lambda^{-2}(Y_{ia}X_{ai}-X_{ia}Y_{ai})\), \(\triangle(j_i,j_a,\lambda)\).

Unreduced: extra \((2j_i+1)\), no \(\hat\lambda\).

### `comm220tts` (locked \(\lambda=0,1,2\))

m-gold (free \(ijab\), all four \(m\)):

\[
Z_0=\frac14\sum_{ijab,m}n_i n_j\bar n_a\bar n_b\bigl(X_{ijab}Y_{abij}-X_{abij}Y_{ijab}\bigr).
\]

Reduced. \(J_0\) on \((ij)\), \(J_1\) on \((ab)\):

\[
Z_0
=\frac14\sum_{ijab,J_0 J_1}
n_i n_j\bar n_a\bar n_b
\frac{(-1)^{J_0+J_1+\lambda}}{\hat\lambda^2}
\bigl(
\tilde X^{J_0 J_1}_{ijab}\tilde Y^{J_1 J_0}_{abij}
-
\tilde X^{J_1 J_0}_{abij}\tilde Y^{J_0 J_1}_{ijab}
\bigr).
\]

```
GetTBME_J(J0, J1, i, j, a, b)   // X_ijab
GetTBME_J(J1, J0, a, b, i, j)   // Y_abij and X_abij
GetTBME_J(J0, J1, i, j, a, b)   // Y_ijab
```

Unreduced: \(\tfrac14\sum_J(2J+1)\) same \(J\) on bra/ket. NAS + \(i\le j,\,a\le b\): **no** \(\tfrac14\). Do not mix tilde+free+\(\tfrac14\) with NAS+\(i\le j\).

### `comm111tts` (locked \(\lambda=0,1,2\))

m-gold: \(Z_{ij}=\dfrac{1}{2j_i+1}\sum_m\sum_{a,m_a}(X_{ia}Y_{aj}-Y_{ia}X_{aj})\).

AMC: \(\hat\jmath_i^{-2}\hat\lambda^{-1}\) and \((-1)^{j_i+j_a+\lambda}\). Code \(\hat\lambda^{-2}\), **keep** \(\lambda\) in the phase:

\[
Z_{ij}
=\delta_{j_i j_j}\sum_a
\frac{(-1)^{j_i+j_a+\lambda}}{(2j_i+1)\,\hat\lambda^{2}}
(Y_{ia}X_{aj}-X_{ia}Y_{aj}).
\]

### `comm121tts` / `comm221tts` (locked \(\lambda=0,1,2\))

Same leftover rule as 111 (average \(m_i=m_j\)). Code \(\hat\lambda^{-2}\). 221: restore AMC’s omitted \(\tfrac12\). 2b hats: \(J_0\) on \((ci)/(cj)\), \(J_1\) on \((ab)\); `GetTBME_J(J1,J0,a,b,c,j)` for \(Y_{abcj}\).

221 AMC phase is \((-1)^{J_0+J_1+\lambda}\). The m-average wants \((-1)^{J_0+J_1}\) (drop \(\lambda\)). Even \(\lambda\) hides it. **Do not copy this drop onto 111/220/222** — those keep \(\lambda\) in the phase.

### `comm122tts` (locked \(\lambda=0,1,2\); leftover 2b; 1b\(\times\)2b)

Wick = `comm122ss` (no occupancy on contracted \(a\)). AMC: `learn/amc_tts/comm_tts/output/comm122tts_unred.tex`. \(Z=\mathrm{Eq1}-\mathrm{Eq2}\). Unreduced scalars: `MakeReduced` then the same kernel (no `comm122ss` fallthrough).

m-gold. All contracted \(m_a\); leftover-2b CG project vs tilde `GetTBME_J`:

\[
Z(m)_{ijkl}=\sum_{a,m_a}\bigl(
X_{ia}Y_{ajkl}+X_{ja}Y_{iakl}-Y_{ijal}X_{ak}-Y_{ijka}X_{al}
-(X\leftrightarrow Y)
\bigr),
\qquad
\tilde Z^{J}_{ijkl}
=\frac{1}{2J+1}
\sum_{m}
\mathrm{CG}(j_i m_i j_j m_j|J M)\,
\mathrm{CG}(j_k m_k j_l m_l|J M)\,
Z(m).
\]

AMC leftover 2b for this topology is \(\hat J_0^{-1}\) (not \(\hat J_0^{-2}\) like 222). Convert only \(\hat\lambda^{-1}\to\hat\lambda^{-2}\). Keep \(\hat J_2\). Keep AMC 6j and phases (no \(\lambda\) in the phase). \(\delta_{J_0 J_1}\). NAS \(1/\sqrt{2}\) for \(i=j\) or \(k=l\) on `AddToTBME`.

Eq1 T1: \((-1)^{J_0+j_i+j_j}\hat J_0^{-1}\hat J_2\hat\lambda^{-2}\)
6j \(\{J_0 J_2\lambda; j_a j_i j_j\}\)
`OneBody(i,a)`, `GetTBME_J(J2,J0,a,j,k,l)`.

Eq1 T2: \((-1)^{j_i}(-1)^{J_2+j_a}\hat J_0^{-1}\hat J_2\hat\lambda^{-2}\)
6j \(\{J_0 J_2\lambda; j_a j_j j_i\}\)
`OneBody(j,a)`, `GetTBME_J(J2,J0,i,a,k,l)`.

Eq1 T3: \(-(-1)^{j_l}(-1)^{J_2+j_a}\hat J_0^{-1}\hat J_2\hat\lambda^{-2}\)
6j \(\{J_2\lambda J_0; j_k j_l j_a\}\)
`GetTBME_J(J0,J2,i,j,a,l)`, `OneBody(a,k)`.

Eq1 T4: \(-(-1)^{J_0+j_k+j_l}\hat J_0^{-1}\hat J_2\hat\lambda^{-2}\)
6j \(\{J_2\lambda J_0; j_l j_k j_a\}\)
`GetTBME_J(J0,J2,i,j,k,a)`, `OneBody(a,l)`.

Eq2 = \(X\leftrightarrow Y\). Native antiherm 2b on \(X\), herm \(Y\) \(\to\) herm \(Z_2\).

### `comm231tts` (locked \(\lambda=0,1,2\); leftover 1b; 2b\(\times\)3b)

Wick = `comm231ss`. AMC: `learn/amc_tts/comm_tts/output/comm231tts_unred.tex` (raw terms, **not** the grouped `threebody_commutator` tex). Leftover = **111**. Occupancy \(n_a n_b\bar n_c\bar n_d\). Restore AMC’s omitted \(\tfrac14\). Keep AMC 6j and phases (no \(\lambda\) to drop). Unreduced \(\lambda=0\): `comm231ss`.

m-gold. All contracted \(m_a\ldots m_d\); leftover \(m\)-average:

\[
Z(m)_{ij}
=\frac14\sum_{abcd,m}
n_a n_b\bar n_c\bar n_d
\bigl[
(X_{abcd}Y_{cdiabj}-Y_{abicdj}X_{cdab})
-
(Y_{abcd}X_{cdiabj}-X_{abicdj}Y_{cdab})
\bigr],
\qquad
Z_{ij}=\frac{1}{2j_i+1}\sum_{m_i=m_j}Z(m).
\]

Two couplings (different 6j **and** different \(J\) superscripts on the particle-order swap). AMC \(\hat\lambda^{-1}\); code \(\hat\lambda^{-2}\); keep \(\hat\jmath_i^{-2}\) and \(\hat\jmath_0\hat\jmath_1\).

**A.** \(Y_{cdiabj}^{J_1 j_0,J_0 j_1}\): \(j_0\) couples \((J_1,i)\), \(j_1\) couples \((J_0,j)\).

\[
(-1)^{j_i+j_1+J_0}
\begin{Bmatrix}J_1&\lambda&J_0\\ j_1&j_i&j_0\end{Bmatrix}
\bigl(\tilde X^{J_0 J_1}_{abcd}Y^{J_1 j_0,J_0 j_1}_{cdiabj}
-\tilde Y^{J_0 J_1}_{abcd}X^{J_1 j_0,J_0 j_1}_{cdiabj}\bigr)
\]

```
GetTBME_J(J0, J1, a, b, c, d)                 // X_abcd / Y_abcd
Y3.GetME_pn(J1, j0, J0, j1, c, d, i, a, b, j) // Y_cdiabj
X3.GetME_pn(J1, j0, J0, j1, c, d, i, a, b, j) // X_cdiabj
```

**B.** \(Y_{abicdj}^{J_0 j_0,J_1 j_1}\): \(j_0\) couples \((J_0,i)\), \(j_1\) couples \((J_1,j)\).

\[
(-1)^{j_i+j_1+J_1}
\begin{Bmatrix}\lambda&J_0&J_1\\ j_i&j_1&j_0\end{Bmatrix}
\bigl(X^{J_0 j_0,J_1 j_1}_{abicdj}\tilde Y^{J_1 J_0}_{cdab}
-Y^{J_0 j_0,J_1 j_1}_{abicdj}\tilde X^{J_1 J_0}_{cdab}\bigr)
\]

```
GetTBME_J(J1, J0, c, d, a, b)                 // X_cdab / Y_cdab
Y3.GetME_pn(J0, j0, J1, j1, a, b, i, c, d, j) // Y_abicdj
X3.GetME_pn(J0, j0, J1, j1, a, b, i, c, d, j) // X_abicdj
```

Stretched leftover \(m\) is not a scalar for \(\lambda>0\), \(j>1/2\). Compute every leftover \((i,j)\) independently (like `comm231ss`); do not copy \(Z_{ji}\) from \(Z_{ij}\). For \(n_i\neq n_j\) the m-average is not \(i\leftrightarrow j\) symmetric.

### `comm132tts` (locked \(\lambda=0,1,2\); 1b\(\times\)3b \(\to\) leftover 2b)

Wick = `comm132ss`. AMC: `learn/amc_tts/comm_tts/output/comm132tts_unred.tex`. Occupancy \(n_a-n_b\). \(Z=\mathrm{Eq1}-\mathrm{Eq2}\). Unreduced \(\lambda=0\): `comm132ss`.

3b\(\times\)1b in the pair+third basis:

\[
\langle ij\,J_0\,b|\,O\,|kl\,J_0\,a\rangle\,X_{ab}
=\langle ab\,J_{ab}\,c|\,O\,|de\,J_{de}\,f\rangle\,\langle f|G|c\rangle
\]

with \((ab,c)=(ij,b)\), \((de,f)=(kl,a)\), \(G=X\). Code: `GetME_pn(J0, j0, J0, j1, i, j, b, k, l, a)` and `OneBody(a,b)`. \(j_0=(J_0,j_b)\), \(j_1=(J_0,j_a)\).

m-gold. All contracted \(m_a,m_b\); leftover-2b CG project vs tilde `GetTBME_J`:

\[
Z(m)_{ijkl}=\sum_{ab,m}(n_a-n_b)\bigl(X_{ab}Y_{ijbkla}-Y_{ab}X_{ijbkla}\bigr),
\qquad
\tilde Z^{J}_{ijkl}
=\frac{1}{2J+1}
\sum_{m}
\mathrm{CG}(j_i m_i j_j m_j|J M)\,
\mathrm{CG}(j_k m_k j_l m_l|J M)\,
Z(m).
\]

AMC leftover 2b is \(\hat J_0^{-2}\) (same as `comm132ss` / 222). Convert only \(\hat\lambda^{-1}\to\hat\lambda^{-2}\). Keep \(\hat J_0^{-2}\) and \(\hat\jmath_0\hat\jmath_1\). Phase \((-1)^{J_0+j_b+j_0}\). 6j \(\{j_a j_b\lambda; j_0 j_1 J_0\}\).

He4 \(e_{\max}=1\): use odd parity so \(n_a-n_b\) hp 1b (\(s\leftrightarrow p\)) fires. Native herm 3b (native antiherm disagrees with `GetMscheme` on leftover 2b).

### `comm232tts` (locked \(\lambda=0,1,2\); leftover 2b; 2b\(\times\)3b)

Wick = `comm232ss`. AMC: `learn/amc_tts/comm_tts/output/comm232tts_unred.tex`. Occupancy \(\mathcal N_{abc}=n_a n_b\bar n_c+\bar n_a\bar n_b n_c\). Restore AMC’s omitted \(\tfrac12\) and \((1-P_{ij}^{J_0})\) on terms 1&3, \((1-P_{kl}^{J_1})\) on terms 2&4. \(\delta_{J_0 J_1}\). Unreduced \(\lambda=0\): `comm232ss`.

m-gold. All contracted \(m_a,m_b,m_c\); leftover-2b CG project vs tilde `GetTBME_J`:

\[
Z(m)_{ijkl}=-\tfrac12\sum_{abc,m}\mathcal N_{abc}\bigl(
X_{icab}Y_{abjklc}-X_{jcab}Y_{abiklc}-Y_{ijcabl}X_{abkc}+Y_{ijcabk}X_{ablc}
-Y_{icab}X_{abjklc}+Y_{jcab}X_{abiklc}+X_{ijcabl}Y_{abkc}-X_{ijcabk}Y_{ablc}
\bigr),
\qquad
\tilde Z^{J}_{ijkl}
=\frac{1}{2J+1}
\sum_{m}
\mathrm{CG}(j_i m_i j_j m_j|J M)\,
\mathrm{CG}(j_k m_k j_l m_l|J M)\,
Z(m).
\]

AMC leftover 2b for this topology is \(\hat J_0^{-1}\) (same as `comm232ss`), not \(\hat J_0^{-2}\) like 222. Convert only \(\hat\lambda^{-1}\to\hat\lambda^{-2}\). Keep \(\hat J_0^{-1}\), \(\hat J_2\) (or \(\hat J_3\) on \(kl\) terms), \(\hat\jmath_0\hat\jmath_1\). Keep AMC 6j and phases.

Eq1: \((-1)^{j_i}\hat J_0^{-1}(-1)^{J_2+j_1}\hat J_2\hat\jmath_0\hat\jmath_1\hat\lambda^{-2}\)
6j \(\{j_c j_j J_2; j_i j_1 J_0\}\{J_3\lambda J_2; j_1 j_i j_0\}\)
`GetTBME_J(J2,J3,c,j,a,b)`, `GetME_pn(J3,j0,J0,j1,a,b,i,k,l,c)`.

Eq2: \((-1)^{j_k}\hat J_0^{-1}(-1)^{J_2+j_1}\hat J_3\hat\jmath_0\hat\jmath_1\hat\lambda^{-2}\)
6j \(\{j_c j_l J_3; j_k j_0 J_0\}\{J_3 J_2\lambda; j_1 j_0 j_k\}\)
`GetME_pn(J0,j0,J2,j1,i,j,c,a,b,k)`, `GetTBME_J(J2,J3,a,b,c,l)`.

Eq3/4 = \(-\)(\(X\leftrightarrow Y\)) of 1/2. \((1-P_{ij}^{J_0})\) uses \((-1)^{j_i+j_j-J_0}\), not a bare minus.

Native herm 3b (same leftover-2b rule as 132). Native antiherm 3b disagrees with `GetMscheme` on leftover 2b.

### `comm223tts` (locked \(\lambda=0,1,2\); leftover 3b; 2b\(\times\)2b)

Wick = `comm223ss`. AMC: `learn/amc_tts/comm_tts/output/comm223tts_unred.tex` (primitive product only). **Do not use the ss contraction** \(X_{ijna}Y_{kalm}\). AMC is \(X_{ijla}Y_{akmn}\) with a 6j×9j. No occupancy on contracted \(a\). \(\delta_{j_0 j_1}\). Unreduced scalars: `MakeReduced` then the same kernel (no `comm223ss` fallthrough).

Wick \((1-P_{ik}-P_{jk})(1-P_{lm}-P_{ln})\): bra ABC/CBA/ACB (\(1,P_{ik},P_{jk}\)); ket ABC/CBA/BAC (\(1,P_{ln},P_{lm}\)). ss uses ACB on the ket because its primitive is \(X_{ijna}\), not \(X_{ijla}\).

m-gold. All contracted \(m_a\); leftover-3b CG project vs stored `GetME_pn_ch`:

\[
Z(m)_{ijklmn}=(1-P_{ik}-P_{jk})(1-P_{lm}-P_{ln})\sum_{a,m_a}(X_{ijla}Y_{akmn}-Y_{ijla}X_{akmn}),
\qquad
Z^{J_{ab} j,\,J_{de} j}
=\frac{1}{2j+1}
\sum_{m}
\mathrm{CG}(j_i m_i j_j m_j|J_{ab} M_{ij})\,
\mathrm{CG}(J_{ab} M_{ij}\, j_k m_k|j M)\,
\mathrm{CG}(j_l m_l j_m m_m|J_{de} M_{lm})\,
\mathrm{CG}(J_{de} M_{lm}\, j_n m_n|j M)\,
Z(m).
\]

AMC leftover 3b has \(\hat J_1\) (ket pair), no \(\hat\jmath_0^{-2}\). Convert only \(\hat\lambda^{-1}\to\hat\lambda^{-2}\). Keep \(\hat J_1\hat J_2\hat J_3\hat J_4\). Keep AMC 6j, 9j, and phases. Overall minus vs leftover (AMC Eq1 is \(+XY\)).

Eq1: \((-1)^{J_0+j_l+j_m+j_n}(-1)^{J_4+j_a+\lambda}\hat J_1\hat J_2\hat J_3\hat J_4\hat\lambda^{-2}\)
6j \(\{j_n j_m J_4; j_l j_0 J_1\}\)
9j \(\{J_0 j_k j_0; \lambda J_3 J_4; J_2 j_a j_l\}\)
`GetTBME_J(J0,J2,i,j,l,a)`, `GetTBME_J(J3,J4,a,k,m,n)`.

Eq2 = \(X\leftrightarrow Y\). Code: \(-\)(Eq1 − Eq2).

After perm, \(i,j,k\to I_1 I_2 I_3\), \(l,m,n\to I_4 I_5 I_6\), \(J_0\to J_{1p}\), \(J_1\to J_{2p}\):
`GetTBME_J(J1p, J2, I1, I2, I4, a)`, `GetTBME_J(J3, J4, a, I3, I5, I6)`.

Native antiherm 2b on \(X\), herm \(Y\) → herm \(Z_3\).

### `comm222_pp_hhtts` (locked \(\lambda=0,1,2\); pp/hh ladder, tts not ttt)

Wick = `comm222_pp_hhss`. AMC: `learn/amc_tts/comm_tts/output/comm222_pp_hhtts_unred.tex`. Outer leftover = **111** (\(j_i\to J_0\)). TBME packaging = **220**.

m-gold. Wick (all \(m_a,m_b\); restore \(\tfrac12\)):

\[
Z(m)_{ijkl}
=\frac12\sum_{ab,m_a m_b}
(\bar n_a\bar n_b-n_a n_b)
\bigl(X_{ijab}Y_{abkl}-Y_{ijab}X_{abkl}\bigr).
\]

Project leftover onto integer \(J_0\) vs **tilde** `GetTBME_J`:

\[
\tilde Z^{J_0}_{ijkl}
=\frac{1}{2J_0+1}
\sum_{m}
\mathrm{CG}(j_i m_i j_j m_j|J_0 M)\,
\mathrm{CG}(j_k m_k j_l m_l|J_0 M)\,
Z(m).
\]

AMC as printed: \(\delta_{J_0 J_1}(-1)^{J_0+J_2+\lambda}\hat J_0^{-2}\hat\lambda^{-1}\) (no \(\tfrac12\)). Code \(\hat\lambda^{-2}\), restore \(\tfrac12\), **keep** \(\lambda\) in the phase:

\[
\tilde Z^{J_0}
=\frac12\sum_{ab J_2}
(\bar n_a\bar n_b-n_a n_b)
\frac{(-1)^{J_0+J_2+\lambda}}{(2J_0+1)\,\hat\lambda^{2}}
\bigl(
\tilde X^{J_0 J_2}_{ijab}\tilde Y^{J_2 J_0}_{abkl}
-(X\leftrightarrow Y)
\bigr),
\quad\triangle(J_0,J_2,\lambda).
\]

```
GetTBME_J(J0, J2, i, j, a, b)   // X_ijab
GetTBME_J(J2, J0, a, b, k, l)   // Y_abkl
```

Then NAS \(1/\sqrt{2}\) for \(i=j\) or \(k=l\) only when `AddToTBME`. Unreduced \(\lambda=0\): `comm222_pp_hhss`.

Dgemm (`Commutator::comm222_pp_hhtts`): same as `ConstructScalarMpp_Mhh`, intermediate channel \(J_2\) with \(\triangle(J_0,J_2,\lambda)\), pref as above (no extra \(\tfrac12\); NAS ket loop already has it).

### `comm222_phtts` (locked \(\lambda=0,1,2\); ph ladder, tts not ttt)

Scalar analog is `comm222_phss` (not `hpss`). Wick string = phss. Pipeline = χ^η Path B **plus** \((1-P_{ij})\). χ^η is not AS.

AMC: `learn/amc_tts/comm_tts/output/comm222_phtts_via_pandya.tex`. Direct Wick in ordinary indices does not recouple (Pandya topology).

1. **Tensor Pandya** of both \(X,Y\) (AMC `barO = -O`, scheme \(((1,-4),(3,-2))\); same 9j as `DoTensorPandya` / χ^η `amc_bar_omega`):
\[
\bar O^{J_0 J_1\lambda}_{ijkl}
=
-(-1)^{J_0+j_i+j_k+\lambda}\hat J_0\hat J_1
\sum_{J_2 J_3}(-1)^{J_2}\hat J_2\hat J_3
\begin{Bmatrix}\lambda&J_0&J_1\\ J_3&j_l&j_k\\ J_2&j_i&j_j\end{Bmatrix}
O^{J_2 J_3\lambda}_{ijkl}.
\]
CC store \(\langle il|ba\rangle=\bar O_{iabl}\).

2. **Couple in CC** (leftover-2b tts; AMC \(\hat\lambda^{-1}\), code \(\hat\lambda^{-2}\); keep \(\lambda\) in the phase):
\[
\bar Z^{J_0}
=
(-1)^{J_0}\hat J_0^{-2}
\sum_{ab J_2}(n_a-n_b)(-1)^{J_2+\lambda}\hat\lambda^{-2}
\bigl(
\bar X^{J_0 J_2}_{iabl}\bar Y^{J_2 J_0}_{bjka}
-(X\leftrightarrow Y)
\bigr).
\]
DGEMM: \(\bar Z[J_0]\mathrel{+}= w_L w_R\,(\bar X[J_0,J_2]\,R_Y-\bar Y[J_0,J_2]\,R_X)\) with occupancy on the contracted \(|ba\rangle\) rows.

3. **Inverse scalar Pandya** \(Z=-\bar Z\) (like `scalar_pandya_inv`; χ^η used \(+\bar\chi\)) then **\((1-P_{ij})\)**:
\[
\tilde Z^{J}
=
-(1-(-1)^{j_i+j_j-J}P_{ij})
\sum_{J'}(2J'+1)
\begin{Bmatrix}j_i&j_j&J\\ j_k&j_l&J'\end{Bmatrix}
\bar Z^{J'}_{ilkj}.
\]

m-gold: phss Wick (all \(m_a,m_b\)), leftover-2b CG project vs tilde `GetTBME_J`. Unreduced \(\lambda=0\): `comm222_phss`.

---

## What this is not

- **EOM norm** `GetVSEOM_Overlap_single`: \(X_{ai}X_{ai}/(2\lambda+1)\), no \(XY-YX\). Not \([X,Y]_0\).
- **Neithan eq. 279**: reduced coupled product, one \(\hat\lambda\). Overlap \(=\mathrm{Neithan}(O{=}X)/\hat\lambda\).
- **Neithan App. A** 22→0: after \(\tfrac14\) on 220 gold, gold \(=\) **−Neithan**.

---

## AMC vs code (how to compare)

This convert is **global**. Diagrams only name leftover vs contracted legs, occupations, and whether AMC dropped \(\tfrac12\) or \(\tfrac14\).

| AMC `.tex` | code / m-gold |
|---|---|
| \(\hat\lambda^{-1}\) (coupled product) | \(\hat\lambda^{-2}\) |
| leftover \(\hat\jmath_i^{-2}\) / \(\hat J_0^{-2}\) | keep \(1/(2j_i+1)\) / \(1/(2J_0+1)\) |
| omitted \(\tfrac12\) or \(\tfrac14\) | restore to match the ss Wick |
| \((-1)^{j+j+\lambda}\) / \((-1)^{J_0+J_2+\lambda}\) | keep \(\lambda\) (111, 220, 222 pp/hh and ph) |
| 221 \((-1)^{J_0+J_1+\lambda}\) | drop \(\lambda\) (221 only) |
| `reduce=true` on rank-0 \(Z\) | still \(\hat\lambda^{-1}\) vs \(\hat\lambda^{-2}\); name both sides |

Do **not** retune AMC input to force \(\hat\lambda^{-2}\).

---

## Fast path (`Commutator::comm220tts`)

Same formula as the reference. Stored `MatEl` (`ch_bra\le ch_ket`), NAS hh×pp. Same channel: \(\mathrm{accu}(n X_{\mathrm{hp}}\bar n\circ Y_{\mathrm{ph}}^{\mathrm{T}})-(X\leftrightarrow Y)\). Off-diagonal: reverse ME is \(h(-1)^{J_1-J_0}\) times the stored element, not a transpose.

---

## Derivation sketch: why \(\hat K^{-1}\) vs \(\hat K^{-2}\)

Edmonds WE (this code; comment on `GetMschemeMatrixElement_1b` has the CG **order** wrong, the **code** is right):

\[
\langle p m_p|A^{(K)}_M|q m_q\rangle
=\frac{\mathrm{CG}(j_q m_q,KM|j_p m_p)}{\hat j_p}\langle p\|A\|q\rangle.
\]

**Coupled** \(E^{(0)}\): insert \(1/\hat K\) and \((-1)^{K-M}\). Result \(\hat K^{-1}\) times RMEs. That is AMC.

**Trace / leftover average**: omit those. Result \(\hat K^{-2}\). That is code. Leftover 1b: extra \(\hat\jmath_i^{-2}\). Leftover 2b: extra \(\hat J_0^{-2}\) after pair CGs. 0-body has no leftover \(\hat J_0^{-2}\) (the \(m\)-sum already traces \(M\)).

Neithan’s \(A_M=\hat K^{-1}\sum\langle p\|A\|q\rangle[a_p^\dagger\times\tilde a_q]_M\) is already Edmonds WE: do **not** put another \(\hat K^{-1}\) on `OneBody` / `GetMscheme`.

Do not put extra \((-1)^{K-M}\) on the `GetMscheme` loop: that phase is spent in going from CGs to the J-scheme hats.
