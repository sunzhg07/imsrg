# Status: \(\Gamma^{\mathrm{III}_c}\) / \(\chi^\theta\) (2026-07-29)

**Detailed notes:** [`factored_GIIIc/NOTES.md`](factored_GIIIc/NOTES.md) · [`GIIIc_详细说明.md`](GIIIc_详细说明.md)

## Locked chain

| Step | Status |
|---|---|
| m-scheme χ^θ (T×T→S) | **PASS** `test_chi_theta_mscheme.py` |
| AMC χ^θ / fold (regenerated) | `factored_GIIIc/input/` |
| m ≡ DIRECT ≡ Path B | **PASS** `test_tts_GIIIc_mscheme.py` |
| DIRECT ≡ Path A ≡ Path B ≡ AMC 6j/9j | **PASS** `test_tts_GIIIc_new_direct.py` |

## Packaging (critical)

Bare \(S=\sum(-1)^{J_0+J_2+\lambda}\hat\lambda^{-1}w\,\Omega\Omega\).
\(\chi_{\mathrm{red}}=S/\hat J\), \(\chi_{\mathrm{unred}}=S/\hat J^2\).
**Do not store bare \(S\) as reduced** (fixed 2026-07-29).

**AMC was correct** — the FAIL was packaging (code labeled bare \(S\) as reduced). Follow AMC `reduce` as printed; do not retune AMC. Details: [REDUCED_UNREDUCED.md](REDUCED_UNREDUCED.md), [factored_GIIIc/NOTES.md](factored_GIIIc/NOTES.md).

---

## 1. Physics (ground truth)

**Unfactorized m-scheme** (`diag1_compact_new copy_to reform.tex` ~350–352):

\[
\Gamma^{\mathrm{III}_c}_{ijkl}
= -\tfrac12 \sum_{abcd}
\bigl(\bar n_a\bar n_b n_c + n_a n_b\bar n_c\bigr)
(1-P_{ij})(1-P_{kl})
\Bigl(
\Omega_{abcl}\Omega_{idab}\Gamma_{cjkd}
+ \Omega_{icab}\Omega_{abdl}\Gamma_{djkc}
\Bigr).
\]

**Folded form** (same file ~442, 500–501):

\[
\Gamma^{\mathrm{III}_c}_{ijkl}
= -\tfrac12(1-P_{ij})(1-P_{kl})
\sum_{ab}\chi^\theta_{iabl}\,\Gamma_{bjka},
\]

\[
\chi^\theta_{ijkl}
= \sum_{ab}
\bigl[
f(a,b,k)+f(a,b,j)
\bigr]
\Omega_{ijab}\Omega_{abkl},
\quad
f(a,b,x)=n_a n_b\bar n_x + \bar n_a\bar n_b n_x.
\]

Split strips (same slots):

- \(\chi_k\): occupation on ket-first \(k\)
- \(\chi_j\): occupation on bra-second \(j\)
- \(\chi^\theta=\chi_k+\chi_j\) (Hermitian in ordinary \(J\)-scheme; **not** \(\chi_k+\chi_k^T\) via matrix transpose of one strip alone in the RME table — pack each strip, unreduce, then **add**)

**Storage:** \(\Omega\) (tensor) is **reduced**; scalar \(\Gamma\) and final \(Z=\Gamma^{\mathrm{III}_c}\) are **not reduced**. Intermediate \(\chi^\theta\) used in the fold must be **unreduced** (`MakeNotReduced` after packing RME \(\Omega\otimes\Omega\) products).

---

## 2. J-scheme fold (AMC)

Input: `learn/amc_tts/input/dc_G3c_chi_theta.txt`  
Output: `G3c_chi_theta.tex` (3×6j) / `G3c_chi_theta_ninej.tex` (9j).

Seed Term1 (plus \((1-P)(1-P)\) → four exchange terms):

\[
Z^{J_0}_{ijkl}
= -\tfrac12\,\delta_{J_0 J_1}
\sum_{ab J_2 J_3}
\hat J_2^2\,\hat J_3^2
\begin{Bmatrix}
j_i & j_a & J_2 \\
j_j & J_3 & j_b \\
J_0 & j_k & j_l
\end{Bmatrix}
\chi^{J_2}_{iabl}\,\Gamma^{J_3}_{bjka}
\]

(with \(\chi,\Gamma\) in **non-reduced** / AMC `reduce=false` convention). Equivalent 6j form has an extra mid-\(J_4\) and \(\hat J_4^2\).

---

## 3. Implementations (current code)

| Name | Entry | What it does | Status |
|---|---|---|---|
| **DIRECT (new)** | `ReferenceImplementations::comm223_232_tts_GIIIc` | Build `MakeNotReduced(χ_k)` + `MakeNotReduced(χ_j)` (same slots, **no** `M+=M.t()`), then AMC 9j + `(1-P)(1-P)` | **Reference** |
| **Path A** | `comm223_232_GIIIc` with `use_TypeGIIIc_factorized=false` (default) or `*_slow=true` | Delegates to DIRECT | ≡ DIRECT |
| **Path B** | `comm223_232_GIIIc` with `use_TypeGIIIc_factorized=true` | Same χ^θ packaging → Pandya × \(\bar\Gamma\) (IIe-style) → inv Pandya | ≡ DIRECT |
| **AMC 6j/9j debug** | `DebugDirectChiThetaNoPandya` | Same χ^θ Op × Γ via 6j and 9j | ≡ DIRECT; 6j ≡ 9j |
| **Old dual-strip** | resurrected only in `DebugOldDirectReducedGamma` | RME `ChiTab` × Γ with Term1/Term2 sixjs from unfactorized AMC `G3c.tex` | **Wrong** (see §5) |

Flags (`FactorizedDoubleCommutator_eths`):

- `SetUse_TypeGIIIc_factorized(bool)` — Path B vs Path A/DIRECT  
- `SetUse_TypeGIIIc_which_term(0|1|2)` — χ_k+χ_j / χ_k / χ_j  
- `SetUse_TypeGIIIc_slow` — force DIRECT  
- `use_TypeGIIIc_single_chi` — unused (old diagnostic)

χ build helpers: `FillChiThetaG3c` / `FillChiThetaG3c_DGEMM`, `ChiThetaToScalarOperator` (two-Op add after `MakeNotReduced`).

---

## 4. Benchmarks (agreeing family)

`run/test_tts_GIIIc_new_direct.py`, `DebugDirectChiThetaNoPandya`:

| Comparison | λ=0 | λ=2 |
|---|---|---|
| DIRECT ↔ Path A | PASS (0) | PASS (0) |
| DIRECT ↔ Path B | PASS (~1e-14) | PASS (~1e-14) |
| DIRECT ↔ AMC-9j | PASS (0 / ~1e-14) | PASS |
| AMC-6j ↔ AMC-9j | PASS (~1e-14) | PASS |
| χ packaging: two-Op add ↔ pack(χ_k+χ_j) | PASS (0) | PASS |
| χ DGEMM ↔ AMC orbit loop | PASS (0) | PASS |

Representative norms (one random seed): DIRECT/Path B \(\|Z\|_2 \approx 181\) (λ=0), \(\approx 223\) (λ=2).

---

## 5. Failed / obsolete paths

### 5.1 Old DIRECT / old Path A dual-strip

Unfactorized AMC `G3c.tex` Term1+Term2 sixjs on RME `ChiTab` × **unreduced** \(\Gamma\), writing unreduced \(Z\).

| vs Path B / new DIRECT | λ=0 ratio \(\|Z_{\mathrm{old}}\|/\|Z_{\mathrm{B}}\|\) |
|---|---|
| As originally coded | **0.549 FAIL** |
| `MakeReduced(Γ)` → old AMC → `MakeNotReduced(Z)` | **0.561 FAIL** |
| `MakeReduced(Γ)` → old AMC, leave \(Z\) as-written | **0.938 FAIL** |

**Conclusion:** Not fixed by a simple reduce/unreduce sandwich on \(\Gamma\) and \(Z\). Mixing reduced \(\Omega\) products with unreduced \(\Gamma\) in that sixj string is inconsistent, **and** that dual-strip kernel is not the same J-coupling as the folded χ^θ×Γ (9j / Pandya) path.

### 5.2 Wrong χ merges (historical)

| Merge | Result |
|---|---|
| Path A Term1 sixj only with \(x=\chi_k+\chi_j\) (no Term2 ladder) | FAIL (~0.4–1.2×) |
| Factorized outer-leg `CHI_IV` = \(\eta\eta_c+(\eta\eta_d)^T\) (different occ) | ≠ AMC χ^θ |
| `M+=M.t()` on χ_k alone as substitute for χ_k+χ_j | Wrong packaging; current code uses **two Ops + add** |

### 5.3 Pandya “wrong” era

Earlier Path B failed vs **old** DIRECT (~1.8×) because old DIRECT was the buggy reference. With new DIRECT, Path B **PASS**.

---

## 6. Practical recipe (production)

1. Build \(\chi_k,\chi_j\) from reduced \(\Omega\otimes\Omega\) (DGEMM or loop), with \(f(a,b,k)\) / \(f(a,b,j)\).  
2. Pack each as reduced scalar channel MEs → `MakeNotReduced` → **add** → \(\chi^\theta\).  
3. Contract with **unreduced** \(\Gamma\):
   - **DIRECT / Path A:** AMC 9j fold, or  
   - **Path B:** Pandya × \(\bar\Gamma\) → inv Pandya.  

Both give the same unreduced \(\Gamma^{\mathrm{III}_c}\).

---

## 7. Tests / AMC artifacts

| File | Role |
|---|---|
| `run/test_tts_GIIIc_new_direct.py` | DIRECT ↔ Path A ↔ Path B |
| `run/test_tts_GIIIc_chi_theta_amc.py` | AMC 6j/9j vs DIRECT |
| `run/test_tts_GIIIc_old_direct_reduced.py` | Old dual-strip ± reduce sandwich |
| `learn/amc_tts/input/dc_G3c_chi_theta.txt` | Folded m-scheme AMC input |
| `learn/amc_tts/output/G3c_chi_theta*.tex` | 6j / 9j J-scheme |
| `learn/amc_tts/output/G3c.tex` | Unfactorized Term1/Term2 (old DIRECT source) |

---

## 8. Open / not done

- No full m-scheme brute-force numerical cross-check of new DIRECT vs reform.tex loops (J-scheme internal consistency is solid).  
- Old dual-strip not deleted from history; only available as debug for the reduce/unreduce experiment.  
- Factorized scalar `CHI_IV` (outer-leg occ) remains a different object from AMC \(\chi^\theta\).
