# Factorized TTS (ethS) — implemented J-scheme equations

**Curated gold pack (equations + benches):** [tensor_pro_final/README.md](tensor_pro_final/README.md)
— \(f^{\mathrm{I,II,III_a,III_b}}\), \(\Gamma^{\mathrm{I,II,III_a,III_c,IV_a,IV_c}}\).

**Living document.** Update this file whenever a new diagram is validated or an
optimization changes the readable ↔ implemented mapping.

**Before implementing a new diagram:** read [LESSONS.md](LESSONS.md) §Playbook
(product rank → Path A oracle → λ=0 Factorized copy → dual-oracle flags → λ≠0).

**One-body finals (\(f^{\mathrm{I\ldots III}}\)) procedure** (m → AMC direct → Path B → DGEMM):
**[F_1B_PROCEDURE.md](F_1B_PROCEDURE.md)**.

| Status | Diagram | Intermediate | ethS entry | Benchmark |
|---|---|---|---|---|
| **done** | \(f^{\mathrm{I}}\) | \(\chi^\alpha\) (scalar 1b) | `comm223_231_chi1b_tensor` / `use_TypeI_1b` | `run/test_tts_fI.py` (DIRECT↔PathB) |
| **done** | \(f^{\mathrm{II}}\) | \(\chi^\beta\) (tensor 1b) | same / `use_TypeII_1b` | `run/test_tts_fII.py` (DIRECT↔PathB) |
| **done** | \(f^{\mathrm{III}_b}\) | \(\chi^\delta\) (scalar 2b, non-Hermitian) | `comm223_231_chi2b_tensor` / `use_TypeIII_1b` | `run/test_tts_fIIIb.py` (DIRECT↔PathB) |
| **done** | \(f^{\mathrm{III}_a}\) | \(\chi^\gamma\) **scalar**; \(W_1/W_2\) AMC (λ≠0) / Pandya+DGEMM (λ=0) | ethS / `use_TypeIIIa_1b` | `run/test_tts_fIIIa.py` |
| **done** | \(\Gamma^{\mathrm{I}}\) | \(\chi^\varepsilon\) (scalar 1b; Ω×Ω) **ordinary-channel DGEMM** | `comm223_232_chi1b_tensor` / `use_TypeGI_2b` | `run/test_tts_GI.py` (+ `test_tts_GI_dgemm.py`) |
| **done** | \(\Gamma^{\mathrm{II}}\) | \(\chi^\zeta\) (tensor 1b; Γ×Ω) — m ≡ DIRECT ≡ Path B ≡ ethS any λ | same / `use_TypeGII_2b` | `test_chi_zeta_mscheme.py`, `test_tts_GII_{pathB,direct}_mscheme.py`, `test_tts_GII_eths_pathB.py` |
| **done** | \(\Gamma^{\mathrm{III}_a}\) | χ^η: m ≡ DIRECT ≡ Path B; ladder m ≡ from_chi ≡ DGEMM | `comm223_232_GIIIa` | `test_chi_eta_mscheme.py`, `test_GIIIa_ladder_mscheme.py` |
| **done**\* | \(\Gamma^{\mathrm{III}_b}\) | Path B χ̄→RC(pack)→Γ̄·RC→Inv; m ≡ fold ≡ pack (λ=0…4); **code RC ≠ AMC CC print** | `comm223_232_GIIIb` / `use_TypeGIIIb_2b` | `test_G3b_pathB_{fold,pack}_mscheme.py`, `test_G3b_term{1,2}_pathA.py` |
| **done** | \(\Gamma^{\mathrm{III}_c}\) | χ^θ **T×T→S**; DIRECT Case-2 9j; Path B Pandya DGEMM; **m ≡ AMC** | `comm223_232_GIIIc` | `run/test_tts_GIIIc_mscheme.py`, `test_chi_theta_mscheme.py`, `test_tts_GIIIc_new_direct.py` |
| **done** | \(\Gamma^{\mathrm{IV}_a}\) | χ^κ: m ≡ AMC ≡ Path B ≡ ethS | `comm223_232_GIVa` | `test_chi_kappa_*`, `test_G4a_pathB_mscheme.py`, `test_tts_GIVa_eths_vs_pathB.py` |
| **done**\* | \(\Gamma^{\mathrm{IV}_b}\) | Path B χ→fold ≡ m; Fac RC λ=0; rectangular Fac RC λ≠0 open (fold gold) | `comm223_232_GIVb` | `test_G4b_pathB_fold_mscheme.py`, `test_G4b_factorized_fullZ.py` |
| **done**† | \(\Gamma^{\mathrm{IV}_c}\) | χ^λ: ethS Path B Pandya→DGEMM (≡ m ≡ AMC direct) | `comm223_232_GIVc` | `test_chi_lambda_mscheme.py`, `test_tts_GIVc_mscheme.py`, `test_tts_GIVc_pathB.py`, `test_tts_GIVc.py` |
| pending | — | — | — | — |

**Master status table (m ≡ AMC-direct ≡ Path B + notes):** [DIAGRAM_LOCK_STATUS.md](DIAGRAM_LOCK_STATUS.md).

**Remaining research (does not block gold):** AMC Pandya→CC print ≡ code RC (III_b / IV_b); rectangular Factorized RC at λ≠0 for IV_b.

**Code:** `src/FactorizedDoubleCommutator_eths.cc`  
**Unfactored reference:** `ReferenceImplementations::comm223_231_tts_*` / `comm223_232_tts_*`  
**Γ^II@λ=0 reference:** scalar `FactorizedDoubleCommutator` (TTS GII angular formula is for reduced tensor Ω)  
**Trusted M-scheme:** `learn/factorized_code_analyze.tex` §unfact (always unreduced).  
**M vs J probe:** `run/test_tts_f_mscheme.py` — compare literal m to
`GetMschemeMatrixElement_1b(Z_J)` (handles reduce/degeneracy only).
For Ω tensor → scalar \(f\), keep only \([\Omega^\lambda\times\Omega^\lambda]^{(0)}\)
(\(=\) AMC \(\hat\lambda^{-1}\)); see §1b notes below.  
**Agent pitfalls:** [LESSONS.md](LESSONS.md) · trust/AMC caveats: [NOTES.md](NOTES.md)  
**Γ^{III_a} experience:** [factored_GIIIa/NOTES.md](factored_GIIIa/NOTES.md) · [LESSONS.md](LESSONS.md) §\(\Gamma^{\mathrm{III}_a}\)  
**Γ^{IV} experience:** [factored_GIV/NOTES.md](factored_GIV/NOTES.md) · LESSONS §\(\Omega\times\Gamma\to\) tensor χ

---

## 1. Setup and conventions

### Operators

| Symbol | Code | Rank | Scalar? | Typical \(h\) |
|---|---|---|---|---|
| \(\Omega\) | `Eta` | \(\lambda =\) `GetJRank()` | no | \(h_\Omega=-1\) (anti-Hermitian) |
| \(\Gamma\) | `Gamma` | 0 | yes | \(h_\Gamma=+1\) |
| \(Z\) / \(f\) | `Z` | 0 | yes | \(h_Z=h_\Gamma\) |
| \(\chi\) from \(\Omega\times\Omega\) | — | **0** | **yes** | \([\Omega\times\Omega]^0\) only |

**Rule:** tensor×tensor → scalar intermediate. No \([\Omega\times\Omega]^{\lambda\neq 0}\). After χ exists, use the scalar Factorized ladder / Pandya / RC path.

Driver: `comm223_231_st` — tensors stay reduced; scalar `Gamma`/`Z` may be
`MakeNotReduced` for consistency with the scalar factorized path.

### Hats and phases

\[
\hat x = \sqrt{2x+1},\qquad
\hat\jmath_i^{-2} = \frac{1}{2j_i+1} = \frac{1}{j_{i,2}+1}.
\]

`modelspace->phase(n)` \(= (-1)^n\). Half-integer \(j\) enter as
`phase((j2_a + j2_b)/2 + …)`.

### Two-body access

- Tensor: `GetTBME_J(J_bra, J_ket, a,b,c,d)` (bra/ket \(J\) may differ).
- Scalar: `GetTBME_J(J, J, …)` only.
- **Normalized** storage: `GetMatrix(ch_bra,ch_ket)(ibra,iket)`.  
  Physical ME \(= \sqrt{(1+\delta_{ab})(1+\delta_{cd})}\times\) (phase)\(\times\) matrix entry.

### Hermitian storage (why transpose is subtle)

| Type | `MatEl` keys | Other triangle |
|---|---|---|
| Hermitian / anti-Hermitian | `ch_bra ≤ ch_ket` | conjugate + \(h\) (+ tensor \((-1)^{J_b-J_k}\)) |
| **Non-Hermitian** | do not rely on conjugate | must store the ME you need |

Never raw-`.t()` a **tensor** \(\Omega\) block (channel pair flips).  
For **scalar** products: `.t()` only on **square** same-channel blocks.

---

## 2. \(f^{\mathrm{I}}\) — \(\chi^\alpha\) (scalar) \(\times\) \(\Gamma\) (scalar)

**Flags:** `use_TypeI_1b`  
**AMC / comments:** unfactored `output/f1.tex`; factorized as below in code.

### Intermediate \(\chi^\alpha\) (scalar one-body)

Require \(j_d = j_e\) (not full `OneBodyChannel`).

\[
\boxed{
\begin{aligned}
\chi^\alpha_{de}
&=
\frac{1}{\hat\jmath_d^{2}}
\sum_{abc\,J_0 J_1}
w_\alpha\,
(-1)^{J_0+J_1+\lambda}\,\hat\lambda^{-1}\,
\Omega^{J_0 J_1\lambda}_{cdab}\,
\Omega^{J_1 J_0\lambda}_{abce},
\\
w_\alpha
&=
\bar n_a\bar n_b n_c n_d
-
n_a n_b\bar n_c\bar n_d.
\end{aligned}
}
\]

Overall \(\tfrac12\) is **not** in \(\chi^\alpha\); it sits in the \(f^{\mathrm{I}}\) write.

### Final \(f^{\mathrm{I}}\) (scalar one-body)

Require \(j_p = j_q\).

\[
\boxed{
\begin{aligned}
Z_{pq}
&\mathrel{+}=
\frac{1}{2\,\hat\jmath_p^{2}}
\sum_{de\,J}
\hat J^{2}\,
\chi^\alpha_{de}
\bigl(
\Gamma^{J}_{epdq}
+
\Gamma^{J}_{depq}
\bigr),
\\
Z_{qp}
&\mathrel{+}=
h_Z\,Z_{pq}
\quad(p\neq q).
\end{aligned}
}
\]

### Optimizations / readability traps (\(f^{\mathrm{I}}\))

| In code | Why | Pitfall |
|---|---|---|
| Loop \(J=0\ldots J_{\max}\) for both \(\Gamma\) strings | Second string \(\Gamma_{depq}\) has a different \(J\) window | Shared `Jmin`/`Jmax` from one string **drops** valid \(J\) (~3× too small) |
| \(j_d=j_e\) only, all orbits | Scalar intermediate, not OB channel | Restricting to `(l,j,tz)` drops same-\(j\) different-\(l\) |
| \(\tfrac12\) only on write | Match unfactored bookkeeping | Putting \(\tfrac12\) also in \(\chi\) double-counts |
| OpenMP over outer \(d\) / \(p\) | Speed | — |

---

## 3. \(f^{\mathrm{II}}\) — \(\chi^\beta\) (tensor) \(\times\) \(\Omega\) (tensor) \(\to\) scalar

**Flags:** `use_TypeII_1b`  
**AMC:** `factored_fII/output/chi_beta.tex`, `f2a_from_chi.tex`, `f2b_from_chi.tex`

### Intermediate \(\chi^\beta\) (tensor one-body, rank \(\lambda\))

Require \(\lvert j_d-j_e\rvert \le \lambda \le j_d+j_e\) (**not** \(j_d=j_e\)).  
Store as **reduced** 1b ME (no \(\hat\jmath_d^{-2}\)).

\[
\boxed{
\begin{aligned}
\chi^\beta_{de}{}^{\lambda}
&=
\frac12\,
(-1)^{j_e+\lambda}
\sum_{abc\,J_0 J_1}
w_\beta\,
(-1)^{J_0+j_c}\,
\hat J_0\hat J_1
\begin{Bmatrix}
\lambda & J_1 & J_0 \\
j_c & j_d & j_e
\end{Bmatrix}
\Gamma^{J_0}_{cdab}\,
\Omega^{J_0 J_1\lambda}_{abce},
\\
w_\beta
&=
n_a n_b\bar n_c\bar n_e
-
\bar n_a\bar n_b n_c n_e
\quad\text{(depends on free }e\text{, not }d\text{)}.
\end{aligned}
}
\]

Code folds \((-1)^{j_e+\lambda+J_0+j_c}\) into one `phase(...)`.

### Final \(f^{\mathrm{II}}\) (scalar)

Require \(j_p=j_q\). Assemble **\(f_{2a}+h_\Gamma f_{2b}\)** because AMC cannot orient \(\Omega_{diej}\); use \(\Omega_{ejdi}\) instead (see NOTES hermiticity note).

\[
\boxed{
\begin{aligned}
Z_{pq}
&\mathrel{+}=
\frac{1}{\hat\jmath_p^{2}}
\sum_{de\,J_3 J_4}
(-1)^{j_p+J_4+j_d}\,
\hat J_3\hat J_4\,\hat\lambda^{-1}
\begin{Bmatrix}
J_4 & \lambda & J_3 \\
j_e & j_p & j_d
\end{Bmatrix}
\chi^\beta_{de}{}^{\lambda}
\bigl(
\Omega^{J_3 J_4\lambda}_{eidj}
+
h_\Gamma\,
\Omega^{J_3 J_4\lambda}_{ejdi}
\bigr),
\\
Z_{qp}
&\mathrel{+}=
h_Z\,Z_{pq}
\quad(p\neq q).
\end{aligned}
}
\]

(\(\tfrac12\) already inside \(\chi^\beta\).)

### Optimizations / readability traps (\(f^{\mathrm{II}}\))

| In code | Why | Pitfall |
|---|---|---|
| Independent loops \(J_3,J_4=0\ldots J_{\max}\) | \(\Omega_{eidj}\) vs \(\Omega_{ejdi}\) windows differ | Shared \(J\) bound drops terms |
| \(\chi\) as `arma::mat`, not `Operator` | Temporary reduced tensor 1b | Easy to wrongly force \(j_d=j_e\) |
| Partner \(\Omega_{ejdi}\) + \(h_\Gamma\) | Yutsis / hermiticity | Coding \(\Omega_{diej}\) with wrong phase |
| \(\lambda\) fixed to `Eta.GetJRank()` | AMC may write a dummy sum over \(\lambda_0\) | Do not sum over \(\lambda\) |

---

## 4. \(f^{\mathrm{III}_b}\) — \(\chi^\delta\) (scalar 2b) \(\times\) \(\Gamma\) (scalar)

**Flags:** `use_TypeIII_1b`  
**AMC:** `factored_fIIIb/output/chi_delta.tex`  
**Pattern:** scalar `Chi_222_b` in `FactorizedDoubleCommutator.cc`

### Intermediate \(\chi^\delta\) / \(P\) (scalar two-body, **non-Hermitian**)

Only **same-channel** squares \(\chi(ch,ch)\). Fill **full** `(ibra,iket)` — occupations on the second \(\Omega\) break hermiticity (`SetNonHermitian()`).

Physical ladder (what the sum computes before normalization):

\[
\boxed{
\begin{aligned}
P^{J}_{ij,kl}
&=
\sum_{ab\,J_2}
w_\delta\,
(-1)^{J+J_2+\lambda}\,\hat\lambda^{-1}\,
\Omega^{J\,J_2\lambda}_{ijab}\,
\Omega^{J_2\,J\lambda}_{abkl},
\\
w_\delta
&=
n_a n_b\bar n_k\bar n_l
-
\bar n_a\bar n_b n_k n_l.
\end{aligned}
}
\]

**Storage:** write **normalized** matrix entries
\(P_{\mathrm{norm}} = P_{\mathrm{phys}}/\bigl(\sqrt{1+\delta_{ij}}\sqrt{1+\delta_{kl}}\bigr)\)
into `GetMatrix(ch,ch)`.

No \(\tfrac14\) inside \(P\) (applied in the contraction as below).  
No \(\hat J^{-2}\) in \(P\) for the TTS path we ship (net \(\hat J^{2}\) cancels vs unfactored).

### Contraction (matrix product → one-body)

For each scalar MatEl key \((ch_b,ch_k)\) with \(J_b=J_k\):

\[
\boxed{
\begin{aligned}
M
&=
\chi(ch_b)\,\Gamma(ch_b,ch_k)
\quad\text{(normalized matrices)},
\\
M
&\leftarrow
\begin{cases}
M+M^{T} & ch_b=ch_k \quad (A+A^{T},\ h_\Omega=-1),\\
M-\Gamma\,\chi(ch_k) & ch_b\neq ch_k \quad \text{(explicit partner; no }M^{T}\text{ into this block)}.
\end{cases}
\\
M
&\leftarrow
2\,M,
\\
Z_{pq}
&\mathrel{+}=
\frac{1}{4\,\hat\jmath_p^{2}}
\sum_{c\,J}
M^{J}_{cp,\,cq},
\qquad
Z_{qp}\mathrel{+}=h_Z Z_{pq}\ (p\neq q).
\end{aligned}
}
\]

Effective one-body content matches analyze
\(\sum(\chi\Gamma+\Gamma\chi)\) with channel-safe transpose, and unfactored
`comm223_231_tts_fIIIb` (validated).

### Optimizations / readability traps (\(f^{\mathrm{III}_b}\))

| In code | Why | Pitfall |
|---|---|---|
| `SetNonHermitian` + full `(ibra,iket)` | Occ breaks hermiticity | Upper-triangle + conjugate **wrong** |
| Only \(\chi(ch,ch)\) | pp ladder like scalar | No need for \(\chi(ch_b,ch_k)\) off-diagonal |
| Physical → normalized `/SQRT2` | Match `GetMatrix` / `GetTBME` | Writing physical MEs into the matrix double-counts \(\sqrt2\) |
| `M*=2` then `0.25` on 1b | \(A+A^{T}\) already doubles topologies; TTS has no extra \((2J+1)\) | Blind copy of scalar `*=4*(2J+1)` → exact **2×** error |
| Unequal-channel `M -= Γ χ_ket` | Rectangular block | `.t()` would be wrong shape / wrong `MatEl` key |
| OpenMP over channels / \(p\) | Speed | Race-free: each `(ch,ch)` / `(p,q)` owned by one thread |

---

## 5. Code map (where to look)

```
comm223_231_st
 ├─ use_1b_intermediates → comm223_231_chi1b_tensor
 │     ├─ use_TypeI_1b   → χ^α, f^I
 │     └─ use_TypeII_1b  → χ^β, f^II
 └─ use_2b_intermediates → comm223_231_chi2b_tensor
       ├─ use_TypeIII_1b   → χ^δ, f^III_b
       └─ use_TypeIIIa_1b  → χ^γ, f^III_a
```

Python flags: `Commutator.FactorizedDoubleCommutator_eths.SetUse_*`.

Unfactored extracts for diffs:  
`comm223_231_tts_fI`, `_fII`, `_fIIIa`, `_fIIIb` in `ReferenceImplementations.cc`.

---

## 6. Validation checklist (when adding a diagram)

1. Unfactored AMC from analyze §unfact → reference function.  
2. Factorized AMC with correct `scalar=` on intermediates / final.  
3. Paste unfactored loops into ethS → must match reference (wiring).  
4. Enable factorization; never share one \(J\) window for swapped legs.  
5. Check hermiticity class of every intermediate (esp. occ-weighted \(\chi\)).  
6. Prefactors: trust analyze rationals; take 6j/hat/phase from AMC.  
7. Add `run/test_tts_*.py` and a row in the status table above.

---

## 7. Changelog

| Date | Change |
|---|---|
| 2026-07-23 | Initial doc: \(f^{\mathrm{I}}\), \(f^{\mathrm{II}}\), \(f^{\mathrm{III}_b}\) as implemented in ethS |
| 2026-07-23 | \(f^{\mathrm{III}_a}\) / \(\chi^\gamma\): AMC closed form + ladder; dense AMC-index table (not `TwoBodyME` GetTBME) |
| 2026-07-23 | \(\Gamma^{\mathrm{III}_a}\): Path-A Term strips + ladder PASS; TTS ref exchange-phase fixes; lessons in `LESSONS.md` + `factored_GIIIa/NOTES.md` |
| 2026-07-24 | \(\Gamma^{\mathrm{IV}_{a,b,c}}\) λ≠0 Path A χ strips PASS vs TTS (GIVb ~60×); lessons in `factored_GIV/NOTES.md` + playbook §6 |
| 2026-07-29 | \(\chi^\eta\): m ≡ AMC direct ≡ Pandya Path B; \(\Gamma^{\mathrm{III}_a}\) ladder: m ≡ `G3a_from_chi` ≡ χ×Γ DGEMM (~1e-14) |

---

## Related notes (detail / AMC dumps)

- [NOTES.md](NOTES.md) — trust boundary, unfactored TTS, name map  
- [LESSONS.md](LESSONS.md) — debugging checklist (incl. \(\Gamma^{\mathrm{III}_a}\) experience)  
- [factored_fII/NOTES.md](factored_fII/NOTES.md) — \(\chi^\beta\) AMC  
- [factored_fIIIb/NOTES.md](factored_fIIIb/NOTES.md) — \(\chi^\delta\) storage / transpose  
- [factored_fIIIa/NOTES.md](factored_fIIIa/NOTES.md) — \(\chi^\gamma\) scalar; \(W_1/W_2\) + Pandya/DGEMM  
- [factored_GIIIa/NOTES.md](factored_GIIIa/NOTES.md) — \(\chi^\eta\) lock, ladder gold (m/from_chi/DGEMM), TTS Path A  

---

## Appendix: \(f^{\mathrm{III}_a}\) / \(\chi^\gamma\) (Pandya topology)

\(\chi^\gamma\) is always **scalar** (rank \(0\)). \(\lambda\) on \(\Omega\) may be nonzero for TTS.

### Slow (AMC) — matches `tts_fIIIa`

AMC closed form in normal coupling (`factored_fIIIa/output/chi_gamma.tex`) + ladder
(`fIIIa_from_chi.tex`). Dense table \(\chi[i,j,k,l,J]\) (AMC index order); do **not**
use `TwoBodyME::GetTBME` for arbitrary \((i,j,k,l)\).

**Production speed for \(\lambda\neq 0\):** precompute \(W_1,W_2\) (factored \(J_2J_3\) /
\(J_4J_5\) recoupling) then form \(\chi\) as an \(W_1 W_2\) product — ~10⁴× vs
unfactored at emax=2. See [factored_fIIIa/NOTES.md](factored_fIIIa/NOTES.md).

\[
\begin{aligned}
\chi^{\gamma\,J}_{ijkl}
&=
-(-1)^{j_j+j_l}
\sum_{ab J_2\ldots J_6 j_0}
(\mathrm{occ})\,(-1)^{J_3+J_5+\lambda}\,
\hat J_2\hat J_3\hat J_4\hat J_5\,\hat J_6^{2}\,\hat j_0^{2}\,\hat\lambda^{-1}
\\
&\quad\times(5~\mathrm{sixj})\,
\Omega^{J_2 J_3\lambda}_{ajkb}\,\Omega^{J_4 J_5\lambda}_{ibal},
\\
f_{pq}
&=
\hat\jmath_p^{-2}
\sum_{abc J}
\hat J^{2}
\bigl(
\Gamma^{J}_{cpab}\,\chi^{\gamma\,J}_{abcq}
-\chi^{\gamma\,J}_{pcab}\,\Gamma^{J}_{abqc}
\bigr).
\end{aligned}
\]

Occ: \(n_a\bar n_b n_j\bar n_k-\bar n_a n_b\bar n_j n_k\).

### Fast (Pandya + DGEMM) — Ω \(\lambda=0\) production

Paper B4c/B5c style (same structure as scalar Factorized IIa/IIc, implemented in
**ethS only**): Pandya(\(\Omega\)) → \(\bar\chi^\gamma=\bar\Omega\,(\mathrm{occ})\,\bar\Omega\) (DGEMM)
→ stay-in-Pandya \(\bar\chi\bar\Gamma\) → IIa/IIc 1b. No inverse Pandya.

Auto-fallback to AMC (\(W_1/W_2\)) when \(\Omega\) has \(\lambda\neq 0\). Mid-J NineJ
Pandya for Factorized convention is sketched but not AMC-calibrated.
Details: [factored_fIIIa/NOTES.md](factored_fIIIa/NOTES.md).

### Optimizations / readability traps (\(f^{\mathrm{III}_a}\))

| In code | Why | Pitfall |
|---|---|---|
| \(\chi^\gamma\) scalar | Always rank \(0\) | `lambda` in code is \(\Omega\)’s rank |
| Dense \(\chi[i,j,k,l,J]\) + \(W_1/W_2\) | AMC index order; kill \(J_2\ldots J_5\) nest | `TwoBodyME` GetTBME exchange phases ≠ AMC |
| Ladder \(\Gamma\chi-\chi\Gamma\) with \(\hat J^{2}\) | AMC `fIIIa_from_chi` | Blind copy of \(f^{\mathrm{III}_b}\) (`M*=2`, `0.25`, \(A+A^{T}\)) |
| Term2 uses \(\chi_{pcab}\Gamma_{abqc}\) | Not \((\chi\Gamma)_{cp,cq}\) | Index order on outer legs |
| Fast stay-in-Pandya + DGEMM (λ=0) | Paper B4c; large \(N\) | Do not edit scalar `FactorizedDoubleCommutator.cc` for TTS |

---

## 6. \(\Gamma^{\mathrm{I}}\) — \(\chi^\varepsilon\) (scalar) \(\times\) \(\Gamma\) (scalar)

**Flags:** `use_TypeGI_2b` (default **ordinary-channel DGEMM**); `use_TypeGI_slow` for orbit/`J` loops  
**Driver:** `comm223_232` (tensor \(\Omega\) stays reduced)  
**Ref:** `comm223_232_tts_GI`  
**AMC / m-scheme:** [factored_GI/NOTES.md](factored_GI/NOTES.md)  
**Benches:** `test_chi_epsilon_mscheme.py` (m≡AMC χ); `test_tts_GI_mscheme.py` (AMC PathB≡ethS); `test_tts_GI.py` (ethS≡tts)

### Intermediate \(\chi^\varepsilon\) — DGEMM (not Pandya)

\[
W_c = n_c\,(\bar n_a\bar n_b)+\bar n_c\,(n_a n_b)
\quad\Rightarrow\quad
T_{\mathrm{pp}}=\Omega\,W_{\mathrm{pp}}\,\Omega_{\mathrm{R}},\quad
T_{\mathrm{hh}}=\Omega\,W_{\mathrm{hh}}\,\Omega_{\mathrm{R}}.
\]

Per channel pair \((J_0,J_1)\): **normalized** `GetMatrix` / MatEl (same as scalar χ^α) →

\[
T = 2\cdot\mathrm{ang}\cdot\mathrm{Left}\,W\,\mathrm{Right},
\quad
\chi^\varepsilon_{pq}
\mathrel{+}=
\tfrac12\,\frac{1}{\hat\jmath^{2}}
\sum_c N_{cp}N_{cq}\bigl(n_c T_{\mathrm{pp}}+\bar n_c T_{\mathrm{hh}}\bigr)_{(cp),(cq)}.
\]

Leading **2** replaces the orbit-loop \(a\leftrightarrow b\) double count (channel stores \(a\le b\) once).  
Outer \(N=\sqrt{1+\delta}\) restores unnormalized products (`GetTBME` vs `GetTBME_norm`).  
Skip MatEl-missing \((J,\pi,T_z)\) pairs (do not call `GetMatrix.at`).  
\(\mathrm{ang}=(2J+1)\) (λ=0) or \((-1)^{J_0+J_1+\lambda}\hat\lambda^{-1}\) (tensor).  
Require \(j_p=j_q\) only. Insertions: `GetTBME_J` + √2 (unchanged).

---

## 7. \(\Gamma^{\mathrm{II}}\) — \(\chi^\zeta\) \(\times\) \(\Omega\) (any λ)

**Flags:** `use_TypeGII_2b`; same `use_TypeGI_slow`  
**Path:** χ^ζ (reduced tensor 1b) → \(W-V\) with \((1-P)\) restore (Path B / ethS)  
**AMC / m-scheme:** [factored_GII/NOTES.md](factored_GII/NOTES.md) — **PASS** m ≡ DIRECT ≡ Path B ≡ ethS any λ  
**Bench:** `test_chi_zeta_mscheme.py`; `test_tts_GII_pathB_mscheme.py`; `test_tts_GII_direct_mscheme.py`; `test_tts_GII_eths_pathB.py`
