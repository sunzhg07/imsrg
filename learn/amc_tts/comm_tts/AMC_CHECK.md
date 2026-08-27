# AMC regeneration: comm231tts / comm132tts / comm232tts

**Date:** 2026-08-26  
**Inputs:** `input/*_unred.txt` (production packaging) and `input/*_red.txt` (reference)  
**Outputs:** `output/*.tex`  
**Regenerate:** `./regenerate.sh`

---

## Packaging rules (CRITICAL)

AMC default (`amc/ast.py`): `reduce = reduce or not scalar`

| Operator | AMC declare | IMSRG storage | Notes |
|---|---|---|---|
| X, Y (tensor 2b/3b/1b) | `scalar=false` | WE-reduced | Always reduced in AMC |
| Z (scalar 1b or 2b) | `scalar=true`, **`reduce=false`** | `IsReduced()=false` | **Must set explicitly** |
| Z reduced reference | `scalar=true`, `reduce=true` | \(Z_{\mathrm{red}} = \hat J\, Z_{\mathrm{unred}}\) | For packaging checks only |

**Verified on comm231tts:** unred has \(\hat\jmath_i^{-2}\); red has \(\hat\jmath_i^{-1}\) (same 6j/hats otherwise).

Old `amc/examples/sample_input/comm*.txt` files had **no `reduce=false` on Z**. Outputs happened to show \(\hat\jmath_i^{-2}\) for 231tts (unreduced-like), but this was not documented or enforced. **Always set `reduce=false` for C++ compare.**

Convention: `--wet-convention wigner` (matches IMSRG `GetMscheme` / UnitTest).

---

## comm231tts (tensor 2b × tensor 3b → scalar 1b)

### m-scheme input (Eq 1 + Eq 2 blocks)

```
# Eq1: X(2b) Y(3b)
Z_ij = sum_abcd((n_a*n_b*nbar_c*nbar_d)*(X_abcd*Y_cdiabj - Y_abicdj*X_cdab));
# Eq2: X(3b) Y(2b)
Z_ij = sum_abcd((n_a*n_b*nbar_c*nbar_d)*(X_abicdj*Y_cdab - Y_abcd*X_cdiabj));
```

### Check vs `learn/threebody_commutator equations.tex` (lines 82–118)

| Item | AMC unred | threebody tex | C++ note |
|---|---|---|---|
| \(\delta_{j_j,j_i}(-1)^{j_i}\hat\jmath_i^{-2}\) | ✓ | ✓ | `diag_phase`, `diag_factor` |
| Occ \(n_a n_b \bar n_c \bar n_d\) | ✓ | ✓ | |
| Eq1 phase / 6j | \((-1)^{J_0+j_1}\), 6j\((J_1,\lambda,J_0;j_1,j_i,j_0)\) | ✓ | Term1 X2×Y3 |
| Eq1 exchange | Term2 (AMC sub-term) | same 6j/phase as Term1 in tex | AMC splits `-Y_abicdj X_cdab` into sub-term with **reordered** 6j index labels — algebraically same sum |
| Eq2 X3×Y2 | \((-1)^{J_1+j_1}\), 6j\((\lambda,J_0,J_1;j_i,j_1,j_0)\) | ✓ | |
| Overall **1/4** | **not in AMC** | explicit in tex | Add in C++ / UnitTest m AS |
| **`reduce=false` on Z** | \(\hat\jmath_i^{-2}\) | unreduced | **Required** |

**Status: PASS** (unreduced AMC direct = reference tex up to overall 1/4).

---

## comm132tts (tensor 1b × tensor 3b → scalar 2b)

### m-scheme input

```
# Eq1 (+): X(1b) Y(3b)
Z_ijkl = sum_ab((n_a-n_b)*X_ab*Y_ijbkla);
# Eq2 (−): Y(1b) X(3b)  — subtract this equation from Eq1
Z_ijkl = sum_ab((n_a-n_b)*Y_ab*X_ijbkla);
```

### Check vs threebody tex (lines 351–373)

| Item | AMC unred | threebody tex |
|---|---|---|
| \(\delta_{J_1,J_0}(-1)^{J_0}\hat J_0^{-2}\) | ✓ | ✓ |
| Occ \((n_a-n_b)\) | ✓ | ✓ (old input **omitted** occ) |
| Phase \((-1)^{j_b+j_0}\) | ✓ | ✓ |
| 6j \((j_a,j_b,\lambda;j_0,j_1,J_0)\) | ✓ | ✓ |
| \(X_{ab}^{\lambda}\), \(Y_{ijbkla}^{J_0 j_0 J_0 j_1 \lambda}\) | ✓ | ✓ |
| Full formula | **\(Z = \mathrm{Eq1} - \mathrm{Eq2}\)** | Eq2 has leading minus |

**Old bug:** `comm132tts.txt` had no \((n_a-n_b)\) and spurious `sum_abc`.

**Status: PASS** (with Eq1 − Eq2 convention).

---

## comm232tts (tensor 2b × tensor 3b → scalar 2b)

### m-scheme input (4 equations, N_abc included)

```
N = (n_a*n_b*nbar_c + nbar_a*nbar_b*n_c)
# Block A: X(2b) Y(3b)
Z_ijkl = sum_abc(N*X_cjab*Y_abiklc);
Z_ijkl = sum_abc(N*Y_ijcabk*X_abcl);
# Block B: X(3b) Y(2b)
Z_ijkl = sum_abc(N*Y_cjab*X_abiklc);
Z_ijkl = sum_abc(N*X_ijcabk*Y_abcl);
```

### Check vs threebody tex (lines 246–292)

| Item | AMC unred | threebody tex | Post-AMC |
|---|---|---|---|
| \(\delta_{J_1,J_0}\) | ✓ | ✓ | |
| \((-1)^{j_i}\) / \((-1)^{j_k}\) on [1,3] / [2,4] | ✓ | ✓ | |
| \(\hat J_0^{-1}\) (unred 2b scalar) | ✓ | \(\hat J_0^{-1}\) in tts | red version → \(\hat J_0^{0}\) |
| \(N_{abc}\) in sum | ✓ | ✓ | old input **omitted** |
| 6j content (2×6j terms) | ✓ | ✓ (index map in tex) | |
| Overall **1/2** | **missing** | explicit | multiply in code |
| **(1−P_ij^{J_0})** on [1,3] | **missing** | explicit | apply in code |
| **(1−P_kl^{J_1})** on [2,4] | **missing** | explicit | apply in code |

Index map (AMC ↔ code): \(J_0\leftrightarrow J_1^{\mathrm{code}}\), \(J_1\leftrightarrow J_2^{\mathrm{code}}\), \(J_2/J_3\leftrightarrow J_3/J_4\), \(j_0\leftrightarrow j_1\), \(j_1\leftrightarrow j_2\).

**Status: PASS** for angular / ME content; **1/2 and Pauli factors still manual** (same as comm232st workflow in threebody tex).

---

## What was wrong with old sample inputs

| File | Issues |
|---|---|
| `comm231_tts.txt` | No `reduce=false`; no doc; otherwise m-string OK |
| `comm132tts.txt` | No occ \((n_a-n_b)\); no `reduce=false` |
| `comm232tts.txt` | No \(N_{abc}\); no `reduce=false` |

None of the three had been m-locked against `GetMschemeMatrixElement_*`.

---

## Production files to use in code port

| Commutator | AMC input | AMC output |
|---|---|---|
| comm231tts | `input/comm231tts_unred.txt` | `output/comm231tts_unred.tex` |
| comm132tts | `input/comm132tts_unred.txt` | `output/comm132tts_unred.tex` |
| comm232tts | `input/comm232tts_unred.txt` | `output/comm232tts_unred.tex` |

Canonical copies synced to `amc/examples/sample_{input,output}/`.

---

## Next: m-scheme lock

1. Build m-gold loops (132/231: include occ; 231 λ≠0: \([\Omega\times\Omega]^{(0)}\) if needed for nested).
2. Compare `GetMschemeMatrixElement_*` after `comm*tts` to m-gold.
3. If ratios \(\sim \hat J^{\pm1}\) → packaging; if not → index windows / ME order.
