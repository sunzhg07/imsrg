# Final check: diagram Path B vs m-scheme / AMC-direct

**Date:** 2026-08-27 (G^{III_c} / G^{IV_c} extracts locked: m ≡ AMC-direct ≡ Path B)  
**Scope:** Tensor factorized double commutator in `FactorizedDoubleCommutator_eths`, locked against m-scheme and AMC-direct J-scheme.

**C++ snapshot:** \(\Gamma^{\mathrm{III}_a}/\Gamma^{\mathrm{III}_b}\) Path B is restored (`BuildChiEtaPathB` / `InvChiEtaRed`, 2n×2n scalar \(\chi^\eta\)). \(\Gamma^{\mathrm{IV}_a}\) Path B is any λ (\(W\) by DGEMM). \(\Gamma^{\mathrm{IV}_b}\) C++ is still Factorized λ=0 only.

**Three layers**

| Layer | Meaning | Typical code |
|---|---|---|
| **A** m-scheme | Physical ⟨m…⟩ equation (golden) | Python in `run/test_*mscheme*.py` via `GetMschemeMatrixElement_*` |
| **B** AMC-direct | Unfactorized J-scheme from AMC | `ReferenceImplementations::comm223_*_tts_*` |
| **C** Path B (ethS) | Factorized Pandya / RC / DGEMM / inverse | `FactorizedDoubleCommutator_eths` |

**Status vocabulary (this file)**

| Status | Meaning |
|---|---|
| **LOCKED** | A ≡ B ≡ C confirmed by a live bench (four \(f\) diagrams, \(\Gamma^{\mathrm{I,II}}\), \(\Gamma^{\mathrm{III}_{a,b,c}}\), \(\Gamma^{\mathrm{IV}_{a,c}}\) as of 2026-08-27). Layer B is the AMC/Path-B formula that matches m — **not** automatically `tts_*`. \(\Gamma^{\mathrm{III}_c}\) is the exception: `tts_GIIIc` DIRECT **is** gold |
| **PARTIAL** | Some links verified; missing a direct A↔B, A↔C, or B↔C, or production not wired to the clean extract |
| **OPEN / FAIL** | Known disagreement or no live bench |
| **DEAD BENCH** | Test exists but dies on a missing pybind toggle (`SetUse_Type*_slow` etc.) |

Line numbers refer to the tree as of this write. Re-grep if the file moves.

**Implementing a new diagram:** see [Skills: locking a Path B diagram](#skills-locking-a-path-b-diagram) (procedure, ratio diagnostics, G^{III_a} pitfalls, checklist).

---

## Summary

| Diagram | χ | ethS Path B (C) | AMC-direct (B) | Primary benches | Status |
|---|---|---|---|---|---|
| \(f^{\mathrm{I}}\) | \(\chi^\alpha\) | `…_eths.cc` **254–427** | `ReferenceImplementations.cc` **5309** | `test_tts_f_mscheme.py`, `test_tts_fI.py` | **LOCKED** |
| \(f^{\mathrm{II}}\) | \(\chi^\beta\) | **428–655** | **5413** | `test_tts_f_mscheme.py`, `test_tts_fII.py` | **LOCKED** |
| \(f^{\mathrm{III}_b}\) | \(\chi^\delta\) | **656–873** | **5803** | `test_tts_f_mscheme.py`, `test_tts_fIIIb.py` | **LOCKED** |
| \(f^{\mathrm{III}_a}\) | \(\chi^\gamma\) | **882–1157** (slow oracle **874–881**) | **5793** | `test_tts_fIIIa_mscheme.py`, `test_tts_fIIIa_pathB_cc.py`, `test_tts_fIIIa.py` | **LOCKED** (see note) |
| \(\Gamma^{\mathrm{I}}\) | \(\chi^\varepsilon\) | **1242–1485** | **6012** | `test_tts_GI.py`, `test_tts_GI_mscheme.py`, `test_chi_epsilon_mscheme.py` | **LOCKED** |
| \(\Gamma^{\mathrm{II}}\) | \(\chi^\zeta\) | **1498–1741** | **6195** | `test_tts_GII_pathB_mscheme.py`, `test_tts_GII_eths_pathB.py`, `test_chi_zeta_mscheme.py` | **LOCKED** |
| \(\Gamma^{\mathrm{III}_a}\) | \(\chi^\eta\) | **4338–4440** `comm223_232_GIIIa`; `chi2b` calls it when `use_TypeGIIIa_2b` | AMC `G3a_from_chi` / Chi_AS (not a `tts_*` twin) | `test_chi_eta_mscheme.py`, `test_GIIIa_ladder_mscheme.py`, `test_tts_GIIIa.py` | **LOCKED** |
| \(\Gamma^{\mathrm{III}_b}\) | \(\chi^\eta\)→RC | **4441–4790** `comm223_232_GIIIb` | fold / Path B pack (not `tts_GIIIb`) | `test_tts_GIIIb.py`, `test_G3b_normal_to_RC.py`, `test_G3b_pathB_*` | **LOCKED** (extract; production flag unused) |
| \(\Gamma^{\mathrm{III}_c}\) | \(\chi^\theta\) | **4791–5119** `comm223_232_GIIIc` | **6804 / 6828** `tts_GIIIc` DIRECT (gold) | `test_chi_theta_mscheme.py`, `test_tts_GIIIc_mscheme.py` | **LOCKED** (extract; production flag unused) |
| \(\Gamma^{\mathrm{IV}_a}\) | \(\chi^\kappa\) | Path B `comm223_232_GIVa` (any λ) | AMC analyze / Wbra (not `tts_GIVa`) | `test_chi_kappa_*`, `test_G4a_pathB_mscheme.py`, `test_tts_GIVa_eths_vs_pathB.py` | **LOCKED** (extract) |
| \(\Gamma^{\mathrm{IV}_b}\) | \(\chi^\iota\) | extract **6322–6670** (λ≠0 no-op) | **7412** | `test_tts_GIVb.py`, `test_chi_iota_*`, `test_G4b_*` | **PARTIAL** / **DEAD** benches |
| \(\Gamma^{\mathrm{IV}_c}\) | \(\chi^\lambda\) | Path B **6671–6974** / wrapper **6979** (λ=0 no-op) | ring fold / AMC `G4c` (not `tts_GIVc`) | `test_chi_lambda_mscheme.py`, `test_tts_GIVc_mscheme.py`, `test_tts_GIVc_pathB.py` | **LOCKED** (extract; λ≠0; production flag unused) |

Smoke pack: `learn/amc_tts/tensor_pro_final/run_all_gold.sh` (18 benches).

---

## Shared entry points

| Role | File | Lines |
|---|---|---|
| One-body driver | `src/FactorizedDoubleCommutator_eths.cc` | `comm223_231_st` **191–…** → chi1b **227** + chi2b **630** |
| Two-body driver | same | `comm223_232` **1159–1211** → chi1b **1206** + chi2b **1750** |
| Universal Pandya kernel | same | `TensorPandya` **133–…** (GIIIb extract; GIIIa uses AMC same-label Pandya) |
| Python bindings | `src/pyIMSRG.cc` | ethS submodule **1055–1138** |
| AMC-direct declarations | `src/ReferenceImplementations.hh` | **106–131** |
| AMC-direct implementations | `src/ReferenceImplementations.cc` | see per-diagram |

**Production wiring caveat:** `comm223_232` calls `comm223_232_chi1b_tensor` and `comm223_232_chi2b`. \(\Gamma^{\mathrm{III}_a}\) Path B is wired: `chi2b` calls `comm223_232_GIIIa` when `use_TypeGIIIa_2b`. \(\Gamma^{\mathrm{III}_b}\), \(\Gamma^{\mathrm{III}_c}\), and \(\Gamma^{\mathrm{IV}_c}\) extracts are locked but **not** redirected inside the monolith (`use_TypeGIIIb/c_2b` / `use_TypeGIVc_2b` do not call those extracts; IIb+IId still builds Factorized χ̄ from Ω). `GIVa/b` extracts are also not redirected.

---

## One-body diagrams (confirmed LOCKED 2026-08-27)

### \(f^{\mathrm{I}}\) / \(\chi^\alpha\)

| Piece | Location |
|---|---|
| ethS Path B | `FactorizedDoubleCommutator_eths.cc` **254–427** (`if (use_TypeI_1b)`) |
| Driver | `comm223_231_st` → `comm223_231_chi1b_tensor` |
| AMC-direct | `ReferenceImplementations.cc` **5309** `comm223_231_tts_fI` |
| Bindings | `pyIMSRG.cc` `SetUse_TypeI_1b` **1074**, `comm223_231_st` **1059** |

**Benches**

| Test | What it checks | Result (emax=1 λ=2 / emax=2 for TTS) |
|---|---|---|
| `run/test_tts_f_mscheme.py` | A vs B vs C (`run_pathB` → ethS) | PASS (max\|m−J\| ~ 2e-13) |
| `run/test_tts_fI.py` | B vs C norms | PASS (‖Δ‖ ~ 1e-12 @ emax=2) |

**Status: LOCKED** — m ≡ AMC-direct ≡ ethS Path B.

---

### \(f^{\mathrm{II}}\) / \(\chi^\beta\)

| Piece | Location |
|---|---|
| ethS Path B | `FactorizedDoubleCommutator_eths.cc` **428–655** (`if (use_TypeII_1b)`) |
| AMC-direct | `ReferenceImplementations.cc` **5413** `comm223_231_tts_fII` |

**Benches**

| Test | What it checks | Result |
|---|---|---|
| `run/test_tts_f_mscheme.py` | A vs B vs C | PASS (~1e-15) |
| `run/test_tts_fII.py` | B vs C | PASS (~1e-13 @ emax=2) |
| `run/test_chi_beta_mscheme.py` | χ^β m ≡ J | PASS |

**Status: LOCKED.**

---

### \(f^{\mathrm{III}_b}\) / \(\chi^\delta\)

| Piece | Location |
|---|---|
| ethS Path B | `FactorizedDoubleCommutator_eths.cc` **656–873** (`if (use_TypeIII_1b)`) |
| AMC-direct | `ReferenceImplementations.cc` **5803** `comm223_231_tts_fIIIb` |

**Benches**

| Test | What it checks | Result |
|---|---|---|
| `run/test_tts_f_mscheme.py` | A vs B vs C | PASS (~4e-14) |
| `run/test_tts_fIIIb.py` | B vs C | PASS (~5e-13 @ emax=2) |

**Status: LOCKED.**

---

### \(f^{\mathrm{III}_a}\) / \(\chi^\gamma\)

| Piece | Location |
|---|---|
| ethS slow oracle | `FactorizedDoubleCommutator_eths.cc` **874–881** (`use_TypeIIIa_slow`: W1/W2 + ladder) |
| ethS Path B (production) | **882–1157** (CC Pandya → DGEMM χ̄×Γ̄ → one-body CC trace; **no inverse Pandya**) |
| AMC-direct | `ReferenceImplementations.cc` **5793** `comm223_231_tts_fIIIa` |
| Flag | `SetUse_TypeIIIa_slow` bound at `pyIMSRG.cc` **1086** |

**Benches**

| Test | What it checks | Result |
|---|---|---|
| `run/test_tts_fIIIa_mscheme.py` | **A vs B** (m vs `tts_fIIIa`) | PASS (~2.5e-14) |
| `run/test_tts_fIIIa_pathB_cc.py` | **C vs gold** (CC Path B vs W1/W2 / AMC-direct ladder) λ=0,1,2 | PASS (~1e-14) |
| `run/test_tts_fIIIa.py` | B vs C | PASS (~7e-14 @ emax=2) |
| `run/test_chi_gamma_pathB_amc.py` | χ^γ: m ≡ AMC-direct ≡ AMC Path-B equations | PASS (~3e-15) |

**Status: LOCKED.**

**Note on “transitive”:** unlike \(f^{\mathrm{I/II/III_b}}\), there is no single script that does `GetMschemeMatrixElement_1b(ethS Path B)` against the literal m-scheme sum. The lock is:

1. m ≡ `tts_fIIIa` (`test_tts_fIIIa_mscheme.py`)
2. ethS CC Path B ≡ W1/W2 gold ≡ same AMC-direct ladder family (`test_tts_fIIIa_pathB_cc.py`, `test_tts_fIIIa.py`)

So m ≡ ethS follows from (1)+(2) with the gold path identified with AMC-direct, not from one three-way unpack in a single test. Numerically the chain is closed at machine precision.

**Do not use** `run/test_chi_gamma_mscheme.py` for status — it pins a superseded AMC print and is expected to FAIL.

---

## Two-body diagrams

### \(\Gamma^{\mathrm{I}}\) / \(\chi^\varepsilon\)  — LOCKED

Driver: `comm223_232_chi1b_tensor` **1206–1744**, flag `use_TypeGI_2b` (`do_GI`).

| Piece | Location in `FactorizedDoubleCommutator_eths.cc` |
|---|---|
| Block | **1242–1485** `if (do_GI)` |
| χ^ε DGEMM | **1242–1365** — `T = Ω W Ω` per \((J_0,J_1)\) channel pair, then trace spectator |
| Fold χ×Γ → \(Z\) | **1367–1484** — loop two-body channels, contract 1b χ with 2b Γ (same as scalar `CHI_I × Gamma`) |
| AMC-direct | `ReferenceImplementations.cc` **6012** `comm223_232_tts_GI` |

**Benches (2026-08-27)**

| Test | Check | Result |
|---|---|---|
| `run/test_chi_epsilon_mscheme.py` | χ: m ≡ AMC-direct | PASS (~4e-15) |
| `run/test_tts_GI_mscheme.py` | Path-B χ×Γ fold ≡ ethS | PASS (~1e-14) |
| `run/test_tts_GI.py` | ethS ≡ `tts_GI` | PASS (~1e-11 @ emax=2) |

**Status: LOCKED** — χ by DGEMM, Γ^I by channel 1b×2b contraction.

### \(\Gamma^{\mathrm{II}}\) / \(\chi^\zeta\)  — LOCKED

| Piece | Location in `FactorizedDoubleCommutator_eths.cc` |
|---|---|
| Block | **1498–1741** `if (do_GII)` |
| χ^ζ DGEMM | **1507–1618** — `T = Γ^{J_0} W_{bc} Ω^{J_0 J_1}` then 6j + trace spectator \(a\) |
| Fold χ×Ω → \(Z\) | **1625–1736** — channel loop, \(W-V\) with \((1-P)\); 6j is tensor 1b×2b recoupling (scalar analogue: `CHI_II × Eta`) |
| AMC-direct | `ReferenceImplementations.cc` **6195** `comm223_232_tts_GII` |

**Benches (2026-08-27, after χ^ζ DGEMM)**

| Test | Check | Result |
|---|---|---|
| `run/test_chi_zeta_mscheme.py` | χ: m ≡ AMC-direct | PASS (~2e-15) |
| `run/test_tts_GII_pathB_mscheme.py` | m ≡ Path B \(W-V\) | PASS (~9e-16) |
| `run/test_tts_GII_direct_mscheme.py` | AMC-direct ≡ Path B | PASS |
| `run/test_tts_GII_eths_pathB.py` | ethS ≡ Path B | PASS (~9e-16 @ λ=2; ~7e-15 @ λ=0) |

**Status: LOCKED** — χ by DGEMM, Γ^II by channel 1b×2b contraction (same philosophy as scalar).

### \(\Gamma^{\mathrm{III}_a}\) / \(\chi^\eta\)  — LOCKED

χ^η is scalar (\(T\times T\to S\)), **not** AS. Gold chain:

**m-scheme (unfactorized) ≡ AMC-direct (`G3a_from_chi` / χ Path A) ≡ Path B** (Pandya→inv χ^η, then Chi_AS×Γ DGEMM).

There is no separate `tts_GIIIa` twin and no Factorized IIa n×n path.

**Path B** in `comm223_232_GIIIa` (**4338–4440**); production `chi2b` calls it when `use_TypeGIIIa_2b` (**2304–2306**):

1. **Same-label fwd Pandya** of Ω into 2n CC (`amc_bar_omega`; AMC `chi_eta_via_pandya` eq1).
2. **Mid DGEMM** \(\bar\chi^{J}=\hat J^{-1}(-1)^{J}\sum_{J'}(-1)^{J'+\lambda}\hat\lambda^{-1}\,w\,\bar\Omega^{JJ'}(\mathrm{occ}\odot\bar\Omega^{J'J})\) with \(w=\bar n_a n_b\bar n_k+n_a\bar n_b n_k\).
3. **Inv Pandya** to ordinary \(\chi_{\mathrm{red}}\) (AMC: \(\chi=\hat J_0\sum\hat J'\,6j\,\bar\chi\); **no** \((1-P)\) on χ).
4. **Ladder** \(Z_{\mathrm{red}}=-\mathrm{Chi}_{AS}\Gamma_{\mathrm{unred}}-\Gamma_{\mathrm{unred}}\mathrm{Chi}_{AS}^{T}\), store \(Z_{\mathrm{unred}}=Z_{\mathrm{red}}/\hat J\) (√2 for \(i=j\)/\(k=l\)).  
   \(\mathrm{Chi}_{AS}[ij,ab]=\chi_{ijab}-(-1)^{J+j_i+j_j}\chi_{jiab}\).

| Piece | Location |
|---|---|
| Path B | `FactorizedDoubleCommutator_eths.cc` **4338–4440** `comm223_232_GIIIa` |
| Production | `comm223_232_chi2b` **2304–2306** `if (use_TypeGIIIa_2b) comm223_232_GIIIa(…)` |
| Bindings | `pyIMSRG.cc` `comm223_232_GIIIa` **1113**, `SetUse_TypeGIIIa_2b` **1095** |
| AMC χ / ladder | `learn/amc_tts/factored_GIIIa/input/chi_eta_{direct,via_pandya}.txt`, `G3a_from_chi.txt` |

**Packaging:** Ω reduced tensor; χ and AMC `G3a_from_chi` are reduced; Γ and stored Z are unreduced. \(Z_{\mathrm{red}}=\chi_{\mathrm{red}}\Gamma_{\mathrm{unred}}\) (no extra \(\hat J\) between factors). Bench compares `GetTBME*\hat J` to \(Z_{\mathrm{red}}\). See `REDUCED_UNREDUCED.md`.

**Benches (2026-08-27, emax=1 seed=11)**

| Test | Check | Result |
|---|---|---|
| `run/test_chi_eta_mscheme.py` | χ: m ≡ AMC direct ≡ Pandya Path B | PASS (~1e-15) |
| `run/test_GIIIa_ladder_mscheme.py` | m ≡ AMC from_chi ≡ Chi_AS DGEMM | PASS (~1e-14) |
| `run/test_tts_GIIIa.py` | ethS Path B ≡ ladder DGEMM | PASS (~2e-14 @ λ=2; ~3e-14 @ λ=0) |

**Wall (extract only, OMP 4):** emax=1 λ=2 ~0.4 ms; emax=2 λ=2 ~7 ms.

**Status: LOCKED** — m ≡ AMC-direct ≡ Path B. Production uses the same extract.

### \(\Gamma^{\mathrm{III}_b}\)

| Piece | Location |
|---|---|
| Clean extract | **4441–4790** `comm223_232_GIIIb` |
| Shared χ^η Path B | `BuildChiEtaPathB` / `InvChiEtaRed` (same as GIIIa) |
| AMC-direct | **6569** `tts_GIIIb` — **not gold** |

**Path B (locked):** normal χ^η (not AS / not Hermitian) → Fac Pandya of \(\chi_{\mathrm{unnorm}}=\hat J\chi_{\mathrm{red}}\) (all 2n, no hermiticity) → code RC → \(\bar\Gamma\) DGEMM → Inv + \((1-P_{ij})(1-P_{kl})\).

**Benches:** `test_tts_GIIIb.py` — ethS ≡ Python Path B (JT) **PASS** (~1e-13, λ=0 and λ=2); `test_G3b_normal_to_RC.py` — Python Path B ≡ m; `test_G3b_pathB_fold_mscheme.py` / `test_G3b_pathB_pack_mscheme.py` — fold ≡ m. `test_GIIIb_mscheme.py` still compares m vs `tts_GIIIb` (FAIL; ignore).

**Status: LOCKED** (extract) — m ≡ fold ≡ Python Path B ≡ ethS. Production `chi2b` is not yet redirected to this extract.

### \(\Gamma^{\mathrm{III}_c}\) / \(\chi^\theta\)  — LOCKED (extract)

χ^θ is scalar (\(T\times T\to S\)). Gold chain:

**m-scheme ≡ AMC fold / `tts_GIIIc` DIRECT (Case-2) ≡ Path B** (ordinary-channel χ DGEMM → Fac Pandya × Γ̄ → inv).

Unlike \(\Gamma^{\mathrm{III}_{a,b}}\) and \(\Gamma^{\mathrm{IV}_c}\), **`tts_GIIIc` DIRECT is gold** here. Ignore `run/test_tts_GIIIc.py` (dead `SetUse_TypeGIIIc_slow`).

**Path B** in `comm223_232_GIIIc` (**4791–5119**):

1. **χ DGEMM** `FillChiThetaG3c_DGEMM` (**3883**): pair mats \(T_{\mathrm{pp,hh}}+=(-1)^{J+J'+\lambda}\hat\lambda^{-1}\,\Omega^{JJ'}(\mathrm{occ}\odot\Omega^{J'J})\). Two occ-sided tables: \(\chi_k\) weights ket \(k\), \(\chi_j\) weights bra \(j\). ChiTab stores bare \(S\).
2. **Pack** `ChiThetaToScalarOperator`: \(\chi_{\mathrm{red}}=S/\hat J\) then `MakeNotReduced` → \(\chi_{\mathrm{unred}}=S/\hat J^2\); add the two Ops (same slots, no transpose).
3. **Fac Pandya** of \(\chi^\theta\) and scalar Γ, **DGEMM** \(\bar\chi^\theta\cdot\bar\Gamma\), inv Pandya (IIe-style).

Fold: \(Z=-\tfrac12(1-P_{ij})(1-P_{kl})\sum\chi_{iabl}\Gamma_{bjka}\).

| Piece | Location |
|---|---|
| χ DGEMM | `FactorizedDoubleCommutator_eths.cc` **3883** `FillChiThetaG3c_DGEMM` |
| Path B | **4791–5119** `comm223_232_GIIIc` |
| Bindings | `pyIMSRG.cc` `comm223_232_GIIIc` **1119**, `SetUse_TypeGIIIc_2b` **1101** |
| AMC-direct | `ReferenceImplementations.cc` **6804** `tts_GIIIc`, **6828** `tts_GIIIc_tensor_red` |
| AMC χ / fold | `learn/amc_tts/factored_GIIIc/input/chi_theta{,_reduced}.txt`, `G3c_from_chi*` |

**Packaging:** Ω reduced tensor; Γ and stored Z unreduced. Bare \(S=\sum(-1)^{J_0+J_2+\lambda}\hat\lambda^{-1}w\,\Omega\Omega\). \(\chi_{\mathrm{red}}=S/\hat J\), \(\chi_{\mathrm{unred}}=S/\hat J^2\). Do **not** compare bare \(S\) to AMC unreduced (fake \(\hat J^{\pm2}\)). See `REDUCED_UNREDUCED.md`.

**Benches (2026-08-27, emax=1 λ=2 seed=11)**

| Test | Check | Result |
|---|---|---|
| `run/test_chi_theta_mscheme.py` | χ: m ≡ AMC red & unred | PASS (~5e-15, 120/120) |
| `run/test_tts_GIIIc_mscheme.py` | m ≡ `tts_GIIIc` DIRECT ≡ ethS Path B | PASS (~1e-14; ‖D−B‖~1e-13) |

**Status: LOCKED** (extract) — m ≡ AMC-direct ≡ Path B. Production `chi2b` is not redirected to this extract.

### \(\Gamma^{\mathrm{IV}_a}\) / \(\chi^\kappa\)  — LOCKED (extract)

χ^κ is tensor (\(T\times S\to T\)). Gold chain (any λ, including 0):

**m-scheme ≡ AMC analyze ≡ Path B** (Pandya → VI_II DGEMM → invPlus → \(W=-\chi\Omega\) pair-channel DGEMM → \((1-P)\) on W).

Ω is **WE-reduced at every λ** (λ=0 is the equal-\(J\) limit, not a Factorized CHI_VI fork). **Do not use `tts_GIVa` as layer B.**

**Path B** in `comm223_232_GIVa`:

1. Scalar Pandya Γ (6j) + tensor Pandya Ω (9j, IMSRG `adcb`, **no extra `scale`**).
2. **DGEMM** \(\bar\chi^{J_0 J_1}=h_\Omega(-1)^{J_0+J_1}(\mathrm{occ}_{ABbarD}\odot\bar\Omega^{J_1 J_0})^{T}\bar\Gamma^{J_1}\).
3. InvPlus (AMC Eq4 **without** printed leading minus).
4. \(W_{\mathrm{red}}=-\,(-1)^{J_0}\hat J_0^{-1}\hat\lambda^{-1}\sum_{J_2}\chi^{J_0 J_2}\Omega^{J_2 J_0}(db;kl)\) as pair-channel DGEMM; \((1-P)\) **on W only**; Hermitian \(W+W_{klij}\) (no extra \(h_\Omega\)). Store \(Z_{\mathrm{unred}}=Z_{\mathrm{red}}/\hat J\).

**Benches (2026-08-27, emax=1 seed=11)**

| Test | Check | Result |
|---|---|---|
| `run/test_chi_kappa_m_vs_amc.py` | χ: m ≡ AMC | PASS |
| `run/test_chi_kappa_pathB_vs_direct.py` | χ: Path B ≡ AMC | PASS all λ |
| `run/test_G4a_pathB_mscheme.py` | m ≡ Path B full Z | PASS |
| `run/test_tts_GIVa_eths_vs_pathB.py` | ethS ≡ Path B | PASS (~7e-15 @ λ=2; ~1e-14 @ λ=0) |

**Status: LOCKED** (extract) — m ≡ AMC-direct ≡ Path B, any λ. Production `chi2b` is not redirected to this extract.

### \(\Gamma^{\mathrm{IV}_b}\) / \(\chi^\iota\)

| Piece | Location |
|---|---|
| Extract | **6322–6670** — **returns immediately if λ≠0** |
| AMC-direct | **7412** |

**Benches:** `test_tts_GIVb.py` **DEAD** (`SetUse_TypeGIVb_slow`); `test_chi_iota_m_vs_amc.py` (χ PASS).  
**Status: PARTIAL.**

### \(\Gamma^{\mathrm{IV}_c}\) / \(\chi^\lambda\)  — LOCKED (extract)

χ^λ is tensor (\(T\times S+S\times T\to T\)). Gold chain (λ≠0):

**m-scheme ≡ AMC χ / ring fold ≡ Path B** (χ DGEMM → tensor Pandya → mid-J DGEMM → inv → fermionic AS).

There is no `tts_GIVc` twin for this fold. **Do not use `tts_GIVc` as layer B** (TTS DIRECT ≠ χ-fold). λ=0 Path B is a no-op (`comm223_232_GIVc_pathB` returns; scalar CHI_VII is a separate fork).

**Path B** in `comm223_232_GIVc_pathB` (**6671–6974**); wrapper `comm223_232_GIVc` (**6979**):

1. **χ DGEMM** `FillChiLambdaG4c_DGEMM` (**3993**): pair mats
   \(T_{\mathrm{pp,hh}}=\Gamma^{J_0}\mathrm{diag}(w^{J_0})\Omega^{J_0 J_1}\),
   \(U_{\mathrm{pp,hh}}=\Omega^{J_0 J_1}\mathrm{diag}(w^{J_1})\Gamma^{J_1}\),
   \(\chi=n_l T_{\mathrm{pp}}+\bar n_l T_{\mathrm{hh}}+n_j U_{\mathrm{pp}}+\bar n_j U_{\mathrm{hh}}\).
   AMC `chi_lambda.tex` (unreduced Γ, no 6j).
2. **IMSRG tensor Pandya** (`adcb`): AMC \(\bar\chi(p,b,a,r)=\) IMSRG\((p,r,a,b)\); \(\bar\Omega(a,q,s,b)=\) IMSRG\((a,b,s,q)\). Hats × NineJ (no extra `scale` — that is the locked convention that matches m).
3. **Mid-J DGEMM** \(\bar X^{J}\mathrel{+}=\hat\lambda^{-1}(-1)^{J'+\lambda}\bar\chi^{JJ'}\bar\Omega^{J'J}\).
4. **Inv Pandya** (drop AMC-sample minus) then \(Z=\tfrac12(1-P_{ij})(1-P_{kl})X\). Store \(Z_{\mathrm{unred}}=Z_{\mathrm{red}}/\hat J\).

| Piece | Location |
|---|---|
| χ DGEMM | `FactorizedDoubleCommutator_eths.cc` **3993** `FillChiLambdaG4c_DGEMM` |
| Path B | **6671–6974** `comm223_232_GIVc_pathB` |
| Wrapper | **6979** `comm223_232_GIVc` |
| Bindings | `pyIMSRG.cc` `comm223_232_GIVc` **1137**, `SetUse_TypeGIVc_2b` **1128** |
| AMC χ / fold | `learn/amc_tts/factored_GIV/input/chi_lambda.txt`, `G4c_from_chi{,_reduced}.txt` |
| AMC-direct TTS | `ReferenceImplementations.cc` **7643** `tts_GIVc` — **not gold** |

**Packaging:** χ^λ is always WE-reduced (tensor). AMC χ with unreduced Γ has no extra 6j. Fold: \(S=\sum\mathrm{CG}\,\mathrm{CG}\,Z(m)\), \(Z_{\mathrm{red}}=S/\hat J\), \(Z_{\mathrm{unred}}=S/\hat J^2\). Do **not** compare bare \(S\) to AMC unreduced.

**Benches (2026-08-27, emax=1 λ=2 seed=11)**

| Test | Check | Result |
|---|---|---|
| `run/test_chi_lambda_mscheme.py` | χ: m ≡ AMC WE | PASS (~4e-15, 4000/4000) |
| `run/test_tts_GIVc_mscheme.py` | m ≡ ring fold ≡ ethS Path B | PASS (~4e-15) |
| `run/test_tts_GIVc_pathB.py` | ethS ≡ ring JT | PASS (~8e-15) |

**Status: LOCKED** (extract, λ≠0) — m ≡ AMC-direct (ring / AMC χ) ≡ Path B. Production `chi2b` is not redirected to this extract.

---

## Skills: locking a Path B diagram

Lessons from locking \(\Gamma^{\mathrm{III}_a}\) (and the same class of mistakes on earlier \(f\) / \(\Gamma\) pieces). Longer agent playbook: `learn/amc_tts/LESSONS.md` and `REDUCED_UNREDUCED.md`.

### Procedure (do not skip)

1. **Classify the product before writing loops.** \(T\times T\to S\) → scalar χ (Ω’s λ only enters Pandya of Ω). \(T\times S\to T\) → tensor χ of rank λ. Never invent \(T\times T\to T\).
2. **Name packaging on every tensor** before any numeric compare. AMC `reduce=true` / `IsReduced()=true` = reduced. `GetMschemeMatrixElement_*` = always physical (unreduced). \(X_{\mathrm{red}}=\hat J\,X_{\mathrm{unred}}\). If ratios cluster near \(\hat J^{\pm1}\) or \(1/\sqrt{2J+1}\), stop — packaging, not 6j.
3. **Name the gold.** Layer A is m-scheme. Layer B is the AMC formula that matches m (often `from_chi` / Path B print), **not** automatically `tts_*`. Layer C is the ethS extract. A dead bench (`SetUse_Type*_slow` unbound) is not a lock.
4. **Lock χ alone, then the fold, then C++.** Isolate χ vs m and AMC (Path A direct and Path B Pandya→mid→inv). Isolate the ladder/fold from that locked χ (Python DGEMM is enough). Only then implement C as those two steps with CC storage + DGEMM. Debugging χ and the fold at once wastes days.
5. **Implement C as the locked equations, not as a Factorized twin.** Scalar Factorized IIa/IIc is a *storage* trick for AS operators. Copying it for a non-AS χ (or mixing Factorized χ̄ with AMC inv / AMC gold) is the usual FAIL. After χ is scalar, reuse the *algebra* of the scalar ladder, not n×n `GetMatrix` packing.
6. **Bench C against the same gold as steps 3–4.** For \(\Gamma^{\mathrm{III}_a}\) that is Chi_AS×Γ (`test_tts_GIIIa.py`), not a strip `tts_*`.
7. **Wire production to the extract.** \(\Gamma^{\mathrm{III}_a}\): `chi2b` calls `comm223_232_GIIIa` when `use_TypeGIIIa_2b`. \(\Gamma^{\mathrm{III}_{b,c}}\) and \(\Gamma^{\mathrm{IV}_c}\) flags still do not redirect the monolith.

### Diagnostic: what the ratios mean

| Pattern | Likely cause | What to do |
|---|---|---|
| All ratios \(\approx\hat J^{\pm1}\) or \(1/\sqrt{2J+1}\) | Reduced vs unreduced mixed | Name both sides; fix store / compare (see `REDUCED_UNREDUCED.md`) |
| Scattered ratios (here ~3–15, not one \(\hat J\)) | Wrong χ̄ / occ / inv / legs — **not** packaging | Stop retuning hats. Diff Pandya formula, occupation placement, and inv 6j against the locked χ bench |
| Ratios \(\approx-1\) | Overall sign / \((1-P)\) phase | Check AMC overall minus and \(P_{ij}\) phase \((-1)^{J+j_i+j_j}\) |
| PASS on χ, FAIL on Z | Fold / √2 / GetTBME vs GetMatrix | Keep χ; fix Chi_AS×Γ and identical-orbit √2 |

### \(\Gamma^{\mathrm{III}_a}\) — what actually locked (2026-08-27)

| Piece | Correct | Wrong (do not retry) |
|---|---|---|
| χ^η | AMC same-label Pandya → occ DGEMM → inv; \(w=\bar n_a n_b\bar n_k+n_a\bar n_b n_k\); inv \(\chi=+\hat J_0\sum\hat J'\,6j\,\bar\chi\) | Factorized IIa χ̄: adcb Pandya, `occ_AbarBC` (n on a contracted CC index), inv \(-\sum(2J'+1)6j\) |
| Layout | 2n ordinary χ (not AS). Apply \((1-P)\) **after** building χ (`Chi_AS`, or DGEMM then \(1-P\)) | Pack \((1-P)\) into n×n `GetMatrix` like scalar IIa — drops \(i>j\) / \(a>b\) pieces of χ |
| Fold | \(Z_{\mathrm{red}}=-\mathrm{Chi}_{AS}\Gamma_{\mathrm{unred}}-\Gamma\mathrm{Chi}_{AS}^{T}\); store \(Z_{\mathrm{unred}}=Z_{\mathrm{red}}/\hat J\) | χ_unred×Γ with extra hats; missing overall minus |
| Gold | m ≡ AMC `G3a_from_chi` ≡ Chi_AS DGEMM ≡ ethS Path B | Factorized IIa n×n `GetMatrix`; strip `tts_*` twin |

Tensor Ω only affects step 1 (Pandya of Ω). After χ is scalar, the ladder is ordinary-channel DGEMM — same philosophy as \(\Gamma^{\mathrm{I,II}}\).

**CC product occupation:** AMC \(w\) depends on contracted \((a,b)\) **and** the outer ket index \(k\). That is not `occ_AbarBC` on the left CC pair. Split \(n_k\) / \(\bar n_k\) onto columns of the right \(\bar\Omega\) factor.

**NineJ cache:** fill \(\bar\Omega\) blocks serially (`GetNineJ` must not populate its cache under OMP); DGEMM and inv may be parallel.

### \(\Gamma^{\mathrm{III}_b}\) — what actually locked (2026-08-27)

| Piece | Correct | Wrong (do not retry) |
|---|---|---|
| χ^η | Same Path B as GIIIa (not AS). \(\chi_{\mathrm{unnorm}}=\hat J\chi_{\mathrm{red}}\) | Factorized Ω̄@Ω̄ χ̄ (`occ_AbarBC`); packing χ into n×n GetMatrix |
| Fac Pandya of χ | Invert CC \(\bar\chi^\eta\) **on the fly** (`InvChiEtaRed`) for every CC-oriented \((a,d,c,b)\) | Look up χ from ordinary-channel Chi2 — CC legs are not always ordinary kets; \(\lvert\bar\chi\rvert\) comes out short; ratios scatter |
| Γ̄ | Scalar Fac Pandya (hermiticity OK) | Tensor Pandya of Γ (Γ is scalar) |
| Fold | code RC \(\bar\chi_{bc,ad}+\bar\chi_{ad,bc}\) → \(\bar\Gamma\cdot\mathrm{RC}\) → Inv + \((1-P)^2\) | GIIIa Chi_AS×Γ ladder; `tts_GIIIb` |
| Gold | m ≡ fold ≡ Python Path B ≡ ethS (`test_tts_GIIIb.py`) | `tts_GIIIb`; Factorized IIb χ̄-from-Ω |

### \(\Gamma^{\mathrm{III}_c}\) — what actually locked (2026-08-27)

| Piece | Correct | Wrong (do not retry) |
|---|---|---|
| χ^θ | Ordinary-channel DGEMM `FillChiThetaG3c_DGEMM`; \(\chi=\chi_k+\chi_j\) (add two unreduced Ops) | Store bare \(S\) as reduced; \(\chi_k+\chi_k^T\) via one-strip transpose |
| Packaging | \(\chi_{\mathrm{red}}=S/\hat J\), \(\chi_{\mathrm{unred}}=S/\hat J^2\); pack then `MakeNotReduced` | Compare \(\sum\mathrm{CG}\,\mathrm{CG}\,\chi(m)\) to AMC unreduced |
| Fold | Fac Pandya of χ and Γ → \(\bar\chi\cdot\bar\Gamma\) DGEMM → inv | Orbit-nested χΩΩΓ; dead `*_slow` Path A |
| Gold | m ≡ `tts_GIIIc` DIRECT ≡ ethS (`test_tts_GIIIc_mscheme.py`) | `test_tts_GIIIc.py` (unbound slow toggle) |

### \(\Gamma^{\mathrm{IV}_c}\) — what actually locked (2026-08-27)

| Piece | Correct | Wrong (do not retry) |
|---|---|---|
| χ^λ | Pair-channel DGEMM `FillChiLambdaG4c_DGEMM` (ΓΩ + ΩΓ, occ on outer \(j,l\)) | 8-deep orbit nest; `tts_GIVc` χ |
| Pandya | IMSRG tensor `adcb` hats × NineJ (**no** extra `scale`) | Universal-kernel `scale` (that is a different convention; this diagram matches m without it) |
| Fold | Mid-J \(\bar\chi\cdot\bar\Omega\) DGEMM → inv (no AMC-sample minus) → \(\tfrac12(1-P)^2\) | `tts_GIVc` DIRECT; retune ring to Path A strips |
| Gold | m ≡ AMC χ / ring fold ≡ ethS (`test_tts_GIVc_mscheme.py`, `test_tts_GIVc_pathB.py`) | `tts_GIVc`; λ=0 Path B (no-op — CHI_VII is separate) |

### Checklist for the next diagram

- [ ] Product rank named (\(T\times T\to S\) vs \(T\times S\to T\))
- [ ] Reduced vs unreduced named on Ω, χ, Γ, Z
- [ ] χ locked: m ≡ AMC direct ≡ Path B Pandya (own bench)
- [ ] Fold locked from that χ (own bench) — gold written on the test, not “tts by default”
- [ ] C++ Path B = those two steps (CC + DGEMM), not a silent Factorized copy
- [ ] Extract bench vs that gold; production wire called out if still on the monolith
- [ ] No dead `*_slow` toggle claimed as a lock

---

## Dead benches (AttributeError on missing bindings)

These call toggles that are **not** in `pyIMSRG.cc` / ethS:

`SetUse_TypeGI_slow`, `SetUse_TypeGIIIb_slow`, `SetUse_TypeGIIIc_slow` / `_factorized` / `_single_chi` / `_which_term`, `SetUse_TypeGIVa_slow`, `SetUse_TypeGIVb_slow`, `SetUse_TypeGIVc_slow`

Affected examples: `run/test_tts_GIIIc.py` (+ several GIIIc variants), `test_tts_GIVb.py`, `test_tts_GI_dgemm.py`, `verify_fact_vs_tts_all.py`.

(`test_tts_GIIIa.py` / `test_tts_GIIIb.py` / `test_tts_GIIIc_mscheme.py` / `test_tts_GIVc_pathB.py` — ethS Path B extracts vs locked χ-fold gold; no slow toggle.)

**Not gold:** `tts_GIVc` (DIRECT ≠ ring fold). Use `test_tts_GIVc_mscheme.py` / `test_tts_GIVc_pathB.py`.

Only live `_slow` toggle for ethS today: **`SetUse_TypeIIIa_slow`**.

---

## How to re-confirm the four \(f\) locks and \(\Gamma^{\mathrm{III}_{a,c}}\) / \(\Gamma^{\mathrm{IV}_c}\) Path B

```bash
cd /path/to/imsrg
export PYTHONPATH=build
OMP_NUM_THREADS=4 python3 -B run/test_tts_f_mscheme.py 1 2          # fI, fII, fIIIb: A≡B≡C
OMP_NUM_THREADS=4 python3 -B run/test_tts_fIIIa_mscheme.py 1 2       # fIIIa: A≡B
OMP_NUM_THREADS=4 python3 -B run/test_tts_fIIIa_pathB_cc.py 1 0,1,2  # fIIIa: C≡gold
OMP_NUM_THREADS=4 python3 -B run/test_tts_fI.py                      # fI: B≡C @ emax=2
OMP_NUM_THREADS=4 python3 -B run/test_tts_fII.py
OMP_NUM_THREADS=4 python3 -B run/test_tts_fIIIb.py
OMP_NUM_THREADS=4 python3 -B run/test_tts_fIIIa.py
# Γ^{III_a} Path B (m ≡ AMC χ ≡ ladder ≡ ethS extract)
OMP_NUM_THREADS=4 python3 -B run/test_chi_eta_mscheme.py 1 2
OMP_NUM_THREADS=4 python3 -B run/test_GIIIa_ladder_mscheme.py 1 2
OMP_NUM_THREADS=4 python3 -B run/test_tts_GIIIa.py 1 2
OMP_NUM_THREADS=4 python3 -B run/test_tts_GIIIa.py 1 0
# Γ^{III_c} Path B (m ≡ tts_GIIIc DIRECT ≡ ethS)
OMP_NUM_THREADS=4 python3 -B run/test_chi_theta_mscheme.py 1 2
OMP_NUM_THREADS=4 python3 -B run/test_tts_GIIIc_mscheme.py 1 2
# Γ^{IV_c} Path B (m ≡ ring fold ≡ ethS; not tts_GIVc)
OMP_NUM_THREADS=4 python3 -B run/test_chi_lambda_mscheme.py 1 2
OMP_NUM_THREADS=4 python3 -B run/test_tts_GIVc_mscheme.py 1 2
OMP_NUM_THREADS=4 python3 -B run/test_tts_GIVc_pathB.py 1 2
```

Optional pack: `./learn/amc_tts/tensor_pro_final/run_all_gold.sh 1 2`

---

## Related docs

- `learn/amc_tts/LESSONS.md` — longer agent playbook (G^{III_a} Path B subsection there is older than this file)
- `learn/amc_tts/DIAGRAM_LOCK_STATUS.md` — longer design notes (may be more optimistic than this file on Γ)
- `learn/amc_tts/factored_GIIIa/NOTES.md` — \(\Gamma^{\mathrm{III}_a}\) / \(\chi^\eta\) Path B
- `learn/amc_tts/GIIIc_STATUS.md` — \(\Gamma^{\mathrm{III}_c}\) / \(\chi^\theta\)
- `learn/amc_tts/factored_GIV/NOTES.md` — \(\Gamma^{\mathrm{IV}_c}\) / \(\chi^\lambda\)
- `learn/amc_tts/REDUCED_UNREDUCED.md` — packaging rules
- `learn/amc_tts/F_1B_PROCEDURE.md` — one-body procedure
- `run/test_pandya_lambda0_reduction.py` — proves which Pandya convention reduces to scalar at λ=0
- `run/capture_pathB_baseline.py` — bit-exact fingerprints for GIII/GIV extracts
