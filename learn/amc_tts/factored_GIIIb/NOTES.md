# Factorized \(\Gamma^{\mathrm{III}_b}\) (IIb+IId)

## WRAP-UP LOCKED (2026-08) — full \(Z\) any λ (0…4)

**Gold chain:**
\[
m_{\text{4-index}}
\;\equiv\;
Z\bigl[\text{fold}(\chi^\eta_m)\bigr]
\;\equiv\;
Z\bigl[\text{Path\,B: }\bar\chi\to\mathrm{RC}(\mathrm{pack})\to\bar\Gamma\cdot\mathrm{RC}\to\mathrm{Inv}\bigr]
\]
with
\[
\chi^\eta_{m,\mathrm{red}}
\;\equiv\;
\chi^\eta_{\mathrm{AMC}}
\;\equiv\;
\chi^\eta_{\mathrm{Path\,B}}.
\]

| Check | Bench | Status |
|---|---|---|
| \(\chi^\eta\): \(m_{\mathrm{red}}\equiv\mathrm{AMC}\equiv\mathrm{Path\,B}\) | `run/test_chi_eta_mscheme.py` (GIIIa) | **PASS** all λ |
| Analyze fold ≡ 4-index \(m\) | `run/test_G3b_pathB_fold_mscheme.py` | **PASS** λ=0…4 |
| Path B pack→Inv→\((1-P)^2\) ≡ fold \(Z\) | `run/test_G3b_pathB_pack_mscheme.py` | **PASS** λ=0…4 |
| Factorized Fac Pandya χ̄ (λ=0 code path) ≡ \(m\) | `run/test_G3b_factorized_mscheme.py` | **PASS** λ=0 |
| Term1 Path B: RC₁→Γ̄→Inv ≡ \(W_{1,m}\) | `run/test_G3b_term1_pathA.py` | **PASS** λ=0…4 |
| Term2 Path B: RC₂; bare \(Z_2(ijkl)\equiv W_2(jilk)\); \((1-P)^2\) locks | `run/test_G3b_term2_pathA.py` | **PASS** λ=0…4 |
| Direct IIb+IId \(m\) ≡ BruteForce JT | [AMC_VS_BRUTEFORCE.md](AMC_VS_BRUTEFORCE.md) | **PASS** scalar |
| ethS Path B ≡ Python Path B (JT) | `run/test_tts_GIIIb.py` | **PASS** λ=0,2 (~1e-13) |

**Analyze fold** (arxiv \(\Gamma^{\mathrm{III}_b}\); \(P\) restored only at the end):
\[
W_{ijkl}
=-\sum_{ab}\Bigl(
  \chi^\eta_{bkai}\,\Gamma_{jbla}
  +\chi^\eta_{lajb}\,\Gamma_{aibk}
\Bigr),
\qquad
Z=(1-P_{ij})(1-P_{kl})W.
\]

**Why not literal IV_b `Z[χ_PathB_WE]≡Z[χ_m]`?**
\(\chi^\eta\) is always scalar (\([\Omega\times\Omega]^0\)) but **not antisymmetrized**.
JT-reduced \(\chi\) loses m-content the fold still uses; WE unpack ≠ \(\chi_m\), and
\(Z[W[\chi_{\mathrm{AS}}]]\neq Z[W[\chi_m]]\). IV_b’s \(\chi^\iota\) is a proper tensor (WE OK).
Any-λ gold here is therefore **fold(\(\chi_m\)) ≡ 4-index**, with Path B locked on \(\chi\) at JT.

**Code RC = JT primitive** (λ=0 Factorized dual only): do **not** block gold on AMC printing RC.

**Code RC** (`FactorizedDoubleCommutator.cc` ~L1828):
\[
\overline{\overline\chi}{}^{J}_{ab,cd}
=-\sum_{J'}(2J'+1)(-1)^{j_b+j_c+J'}
\begin{Bmatrix}j_a&j_b&J\\ j_c&j_d&J'\end{Bmatrix}
\bigl(\bar\chi^{J'}_{bc,ad}+\bar\chi^{J'}_{ad,bc}\bigr).
\]

Conventions: strip \(P\) for AMC; restore by hand; pack ≠ \((1-P)\);
\((1-P)^2\) only on final \(W\to Z\); no TTS as gold.

**Optional later (do not block):** AMC-legal composition ≡ code RC; IIb↔IId as pack split;
rectangular Factorized RC at λ≠0.

---

### AMC discipline (important)

**AMC only does angular-momentum coupling.** Index permutations
(\(P_{ij}\), \(ad\leftrightarrow bc\), reading \(W_{jlik}\) vs \(W_{ijkl}\), …) are
**not** AMC’s job:

1. **Strip** the permutation from the m-scheme equation.
2. Feed AMC the same-label scheme change / product.
3. **Add the permutation back by hand** on the printed JT equation.

Do **not** send permuted RHS indices into AMC (that causes Yutsis crashes and is
the wrong tool). Same rule as \(P\)-stripped direct m-scheme elsewhere.

Related: [AMC_VS_BRUTEFORCE.md](AMC_VS_BRUTEFORCE.md) (direct IIb/IId ≡ BruteForce JT).

---

## Why Factorized does Pandya \(\Gamma\) × RC \(\chi^\eta\)

| Object | Hermitian / AS? | Allowed move |
|---|---|---|
| \(\chi^\eta=\bar\Omega\,(\mathrm{occ}\,\bar\Omega)\) | no | **RC only** (AMC: Pandya→CC) |
| \(\Gamma\) | yes | bra↔ket, \(P\) → \(\pm1\), then Pandya |

IIb+IId m-scheme (P stripped), arxiv / BruteForce map:
\[
W_{ijkl}
=-\sum_{ab}\Bigl(
  \chi^\eta_{bkai}\,\Gamma_{jbla}
  +\chi^\eta_{lajb}\,\Gamma_{aibk}
\Bigr).
\]
(Equivalent 4-index form uses distinct IIb/IId occupations on \(\Omega\Omega\Gamma\).)

**Index prep for Factorized (λ=0 dual):**

1. \(\Gamma_{jbla}=-\Gamma_{bjla}\) (antisym bra).
2. Pandya \(((1,-4),(3,-2))\): \(\Gamma_{bjla}\leftrightarrow\bar\Gamma_{b\bar a\,l\bar j}\).
3. Plain Pandya of \(\chi_{bkai}\) is \(\bar\chi_{b\bar i\,a\bar k}\) — **does not** share the
   \((b,\bar a)\) channel with \(\bar\Gamma\). Need **RC** so \(\chi\) sits as
   \(\overline{\overline\chi}_{b\bar a\,i\bar k}\) (or \(a\bar b,\,i\bar k\) in oriented storage).
4. DGEMM over shared \(a\bar b\):
   \[
   W_{j\bar l\,i\bar k}
   =\sum_{ab}\bar\Gamma_{j\bar l\,a\bar b}\,
    \overline{\overline\chi}{}^\eta_{a\bar b\,i\bar k}.
   \]
5. \(W\) is **not** the standard Pandya of \(Z_{ijkl}\) (that would be \(i\bar l\,k\bar j\)).
   InvPandya must use the **modified** layout \(j\bar l\,i\bar k\) plus extra phases /
   \((1-P_{ij})(1-P_{kl})\).

IId is the second arxiv term (\(\chi_{lajb}\Gamma_{aibk}\)); Factorized packs both inside RC
(\(\bar\chi_{bc,ad}+\bar\chi_{ad,bc}\)).

---

## AMC step chain (scalar)

Inputs: `input/steps/` · Outputs: `output/steps/`.

| Step | M-scheme to AMC (**permute stripped**) | Add back by hand | Status |
|---|---|---|---|
| **1** Pandya \(\Omega,\Gamma\) | `barO_ijkl=-Omega_ijkl` | — (same labels) | **PASS** |
| **2** \(\bar\chi^\eta\) product | `sum occ·barO·barO` | — | **PASS** |
| **3** RC (code loop) | JT primitive L1828 (not AMC Pandya→CC) | pack \(\bar\chi_{bcad}+\bar\chi_{adbc}\) | **PASS** as code primitive |
| **4** GEMM | `W_jlik=Σ barG_jlab dbarChi_abik` | index wiring is the physics (kept) | **PASS** |
| **5** InvPandya | `Z_ijkl=-W_ijkl` (`05_inv_noperm.txt`) | read \(W\) at \(jlik\) (+ \(ik\) orientation, \((1-P)(1-P)\)) | **PASS** + restore |
| **6** Full IIb+IId \(Z\) | Factorized ≡ \(m\) | — | **PASS** `test_G3b_factorized_mscheme.py` |

### Step 1 — Pandya (AMC)

\[
\bar\Omega^{J}_{ijkl}
=-\delta_{JJ'}\sum_{J_2}\hat J_2^2
\begin{Bmatrix}j_j&j_k&J\\ j_l&j_i&J_2\end{Bmatrix}
\Omega^{J_2}_{ijkl}
\]
(same for \(\bar\Gamma\)). Code writes \(O_{adcb}\) explicitly; AMC keeps labels and puts the map in the 6j.

### Step 2 — \(\bar\chi^\eta\) product (AMC)

\[
\bar\chi^{\eta J}_{ijkl}
=\sum_{ab}
(\bar n_a n_b\bar n_k+n_a\bar n_b n_k)\,
\bar\Omega^{J}_{iabl}\,\bar\Omega^{J}_{bjka}.
\]

### Step 3 — RC (**code JT primitive**; AMC dual optional)

**Production definition** = Factorized loop (~L1828), overall \(-\), pack
\(\bar\chi_{bc,ad}+\bar\chi_{ad,bc}\):
\[
\overline{\overline\chi}{}^{J}_{ab,cd}
=-\sum_{J'}(2J'+1)(-1)^{j_b+j_c+J'}
\begin{Bmatrix}j_a&j_b&J\\ j_c&j_d&J'\end{Bmatrix}
\bigl(\bar\chi^{J'}_{bc,ad}+\bar\chi^{J'}_{ad,bc}\bigr).
\]

**AMC dual (documentation only)** — noperm (`03_RC_noperm.txt`) then restore swap by hand:
\[
\overline{\overline\chi}{}^{\eta J}_{ijkl}
=-\,(-1)^{J+j_k+j_l}
\sum_{J'}(-1)^{J'}\hat J'^2
\begin{Bmatrix}j_l&j_j&J\\ j_k&j_i&J'\end{Bmatrix}
\bar\chi^{\eta J'}_{ijkl},
\]
then \(\bar\chi_{ijkl}\to(\bar\chi_{bcad}+\bar\chi_{adbc})\). This print is **not** numerically
identical to the code loop; do not treat AMC RC as gold.

### Step 4 — GEMM (AMC)

M-scheme (wiring is physics, not a “permute to strip”):
\[
W_{jlik}=\sum_{ab}\bar\Gamma_{jlab}\,
\overline{\overline\chi}{}^\eta_{abik}.
\]
AMC → channel product. Code: `CHI_III_final = bar_Gamma * barCHI_III_RC`.

### Step 5 — InvPandya (AMC noperm + restore layout)

**To AMC** (`05_inv_noperm.txt`):
\[
Z_{ijkl}=-W_{ijkl}
\quad\text{(standard InvPandya).}
\]
**AMC print:**
\[
Z^{J}_{ijkl}
=-\sum_{J'}\hat J'^2
\begin{Bmatrix}j_l&j_k&J\\ j_j&j_i&J'\end{Bmatrix}
W^{J'}_{ijkl}.
\]
**Add back by hand** (code ~L2724–2849):

- evaluate \(W\) at \(j\bar l,\,i\bar k\) (and the three exchange slots), not \(i\bar l,\,k\bar j\);
- phases \((-1)^{J+j_i+j_j}\), \((-1)^{J'+j_i+j_k}\) (oriented \(ik\) vs \(ki\));
- full \((1-P_{ij})(1-P_{kl})\).

---

## Locked vs open

**Locked (WRAP-UP any λ = 0…4)**

- \(\chi^\eta\) \(m_{\mathrm{red}}\equiv\mathrm{AMC}\equiv\mathrm{Path\,B}\) (GIIIa)
- Analyze fold ≡ 4-index \(m\) (`test_G3b_pathB_fold_mscheme.py`)
- Path B pack→RC→Γ̄→Inv→\((1-P)^2\) ≡ fold (`test_G3b_pathB_pack_mscheme.py`)
- Term1 / term2 Path B duals (`test_G3b_term{1,2}_pathA.py`)
- Factorized Fac χ̄ full \(Z\) ≡ \(m\) at λ=0 (`test_G3b_factorized_mscheme.py`)
- Direct IIb+IId ≡ BruteForce JT (scalar)
- Steps 1–2, 4, 5 (noperm) via AMC; permute/layout restored **by hand**
- **Step 3 = code RC as JT primitive**

## Critical finding: code RC ≠ AMC Pandya→CC

After matching loops, comments, arxiv (185–188), and AMC:

| Transform | Formula |
|---|---|
| **AMC / diag1 / arxiv Pandya→CC** | \(\overline{\overline\chi}_{ijkl}=-\,(-1)^{J+j_k+j_l}\sum(-1)^{J'}\hat J'^2\{j_l j_j J;\,j_k j_i J'\}\,\bar\chi_{ijkl}\) |
| **Code RC** (L1752 comment + L1828 loop) | \(\overline{\overline\chi}^J_{ab,cd}=-\sum(2J'+1)(-1)^{j_b+j_c+J'}\{j_a j_b J;\,j_c j_d J'\}\,(\bar\chi^{J'}_{bc,ad}+\bar\chi^{J'}_{ad,bc})\) |

These are **different** (6j, phase, and code’s swap). Numeric coeff comparison does not collapse to a single phase. AMC schemes `((1,-2),(3,-4))` / `((1,-4),(2,-3))` **Yutsis-crash**.

**Decision (mirror IV_b):** treat code RC as the production JT definition. AMC noperm+hand-restore remains documentation dual, not the gold path.

### Self-consistent production chain

| Step | Source of JT equation | AMC? |
|---|---|---|
| 1 Pandya \(\Omega,\Gamma\) | m-scheme `bar=-O` | **yes** |
| 2 \(\bar\chi^\eta\) product | m-scheme Pandya product | **yes** |
| 3 RC \(\bar\chi\to\) RC | **code loop** (JT primitive) | **no** (≠ AMC CC) |
| 4 GEMM \(\bar\Gamma\times\) RC | m-scheme `W_jlik=Σ…` | **yes** |
| 5 InvPandya from \(jl,ki\) | m-scheme `Z=-W` + restore layout | **yes** + hand |

**Optional (do not block gold)**

1. AMC-legal composition ≡ code RC (still open).
2. IIb↔IId as pack split (like IV_b term1/term2).
3. Rectangular Factorized RC at λ≠0 (fold is gold).

---

## File map

```
input/steps/01_pandya_fwd.txt
input/steps/02_chi_eta_pandya_product.txt
input/steps/03_RC_noperm.txt              # strip swap; restore by hand
input/steps/04_gemm_IIb_wiring.txt
input/steps/05_inv_noperm.txt             # strip jlik; restore by hand
input/steps/06_IIb_factorized_chain.txt
output/steps/*.tex
```
