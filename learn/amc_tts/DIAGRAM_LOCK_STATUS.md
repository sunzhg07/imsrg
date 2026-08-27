# Diagram lock status: m-scheme ≡ AMC-direct ≡ Path B

**Updated:** 2026-08-14  
**Meaning of columns**

| Column | Meaning |
|---|---|
| **m** | Literal analyze / 4-index m-scheme (physical MEs) |
| **AMC-direct** | AMC on the normal / unfactored string (or BruteForce JT twin), same packaging |
| **Path B** | Factorized dual: Pandya / mid-J / RC / DGEMM / fold as documented per diagram |
| **≡** | Numeric lock on emax=1 benches (typically λ=0…2 or 0…4) |

**Caveats that apply globally**

1. Name **reduced vs unreduced** on both sides ([REDUCED_UNREDUCED.md](REDUCED_UNREDUCED.md)). Ratios \(\sim\hat J^{\pm1}\) → packaging, not 6j.
2. Strip \(P\) for AMC; restore \((1-P)\) by hand. Pack inside RC ≠ external \((1-P)\).
3. **Code RC** (III_b / IV_b Factorized) is a **JT primitive**. It is **not** required to equal an AMC Pandya→CC print of RC.

Cross-refs: [FACTORIZED_TTS_IMPLEMENTED.md](FACTORIZED_TTS_IMPLEMENTED.md) · [F_1B_PROCEDURE.md](F_1B_PROCEDURE.md) · [tensor_pro_final/README.md](tensor_pro_final/README.md)

---

## Summary table

| Diagram | χ | m ≡ AMC-direct ≡ Path B | λ covered | Notes |
|---|---|---|---|---|
| \(f^{\mathrm{I}}\) | \(\chi^\alpha\) (S, 1b) | **YES** | any | ethS TypeI; [F_1B_PROCEDURE.md](F_1B_PROCEDURE.md) |
| \(f^{\mathrm{II}}\) | \(\chi^\beta\) (T, 1b) | **YES** | any | ethS TypeII; [factored_fII/NOTES.md](factored_fII/NOTES.md) |
| \(f^{\mathrm{III}_a}\) | \(\chi^\gamma\) (S, 2b NH) | **YES** | any | Path B = CC Pandya+DGEMM, χ̄ folded with \(\bar\Gamma\) **in Pandya** (no inverse, no `MakeNotReduced`); ethS ≡ gold λ=0…4; [factored_fIIIa/NOTES.md](factored_fIIIa/NOTES.md) |
| \(f^{\mathrm{III}_b}\) | \(\chi^\delta\) (S, 2b NH) | **YES** | any | RME+DGEMM same-ch; non-Hermitian χ — no `.t()` fill; [factored_fIIIb/NOTES.md](factored_fIIIb/NOTES.md) |
| \(\Gamma^{\mathrm{I}}\) | \(\chi^\varepsilon\) (S, 1b) | **YES** | any | ordinary-channel DGEMM; [factored_GI/NOTES.md](factored_GI/NOTES.md) |
| \(\Gamma^{\mathrm{II}}\) | \(\chi^\zeta\) (T, 1b) | **YES** | any | m ≡ DIRECT ≡ Path B ≡ ethS; [factored_GII/NOTES.md](factored_GII/NOTES.md) |
| \(\Gamma^{\mathrm{III}_a}\) | \(\chi^\eta\) (S, 2b) | **YES** | any | χ locked; ladder m ≡ AMC from_chi ≡ DGEMM; [factored_GIIIa/NOTES.md](factored_GIIIa/NOTES.md) |
| \(\Gamma^{\mathrm{III}_b}\) | \(\chi^\eta\) → RC pack | **YES**\* | 0…4 | \*Full \(Z\): \(m_{\text{4-index}}\equiv Z_{\text{fold}}\equiv Z[\text{Path B pack}]\). χ: \(m\equiv\mathrm{AMC}\equiv\mathrm{Path\,B}\). **Code RC ≠ AMC Pandya→CC print.** Term2 bare \(Z_2(ijkl)\equiv W_2(jilk)\). Fac λ=0 χ̄ path separate. [factored_GIIIb/NOTES.md](factored_GIIIb/NOTES.md) |
| \(\Gamma^{\mathrm{III}_c}\) | \(\chi^\theta\) (S, 2b) | **YES** | any | Packaging: \(\chi_{\mathrm{red}}=S/\hat J\) (do not store bare \(S\)); [GIIIc_STATUS.md](GIIIc_STATUS.md) |
| \(\Gamma^{\mathrm{IV}_a}\) | \(\chi^\kappa\) (T) | **YES** | any | One Path B any λ; ethS λ≠0; λ=0 Factorized fork = packaging only; [factored_GIV/NOTES.md](factored_GIV/NOTES.md) |
| \(\Gamma^{\mathrm{IV}_b}\) | \(\chi^\iota\) (T) | **YES**\* | fold: 0…2† | \*Full \(Z\): \(m\equiv\mathrm{AMC\,direct}\equiv\mathrm{Path\,B\,fold}\). Fac Pandya→RC→×Ω̄ at **λ=0** locked; rectangular Fac RC at λ≠0 **open** (fold is gold, does not block). [factored_GIV/NOTES.md](factored_GIV/NOTES.md) |
| \(\Gamma^{\mathrm{IV}_c}\) | \(\chi^\lambda\) (T) | **YES** | any | Path B Pandya→DGEMM; +½ in CC path; [factored_GIV/NOTES.md](factored_GIV/NOTES.md) |

† IV_b fold bench documented through λ=0,1,2; χ alone all λ. Extend fold λ=3,4 if needed — not a known FAIL.

---

## Legend for “AMC-direct” per family

| Family | What “AMC-direct” means in the ≡ chain |
|---|---|
| One-body \(f\) | AMC / TTS on unfactored or from_χ string |
| \(\Gamma^{\mathrm{I,II,III_a,III_c,IV_a,IV_c}}\) | AMC direct / Case-2 TTS twin of the diagram (or χ+ladder) |
| \(\Gamma^{\mathrm{III}_b}\) | Direct = **4-index / BruteForce-mapped m** + analyze **fold**; AMC used for χ / mid-J bars. **Not** AMC print of RC |
| \(\Gamma^{\mathrm{IV}_b}\) | AMC bare χ + analyze fold (Path B); Factorized RC is λ=0 speed twin |

---

## Still open (do not block gold)

| Item | Diagram | Status |
|---|---|---|
| AMC Pandya→CC print ≡ code RC | III_b, IV_b | Open / documentation only |
| Rectangular Factorized RC×Ω̄×Inv at λ≠0 | IV_b | Open; fold is gold |
| Fac equal-J χ̄ path at λ≠0 | III_b | Open; Path B mid-J χ̄ is gold |
| C++ ethS production audit vs new III_b benches | III_b | Bench dual locked; wire ethS when ready |

---

## Primary benches (smoke)

| Diagram | Bench(es) |
|---|---|
| \(f^{\mathrm{I\ldots III}}\) | `run/test_tts_fI.py`, `test_tts_fII.py`, `test_tts_fIIIa.py`, `test_tts_fIIIa_pathB_cc.py`, `test_tts_fIIIb.py` |
| \(\Gamma^{\mathrm{I}}\) | `run/test_tts_GI.py`, `test_tts_GI_mscheme.py`, `test_tts_GI_dgemm.py` |
| \(\Gamma^{\mathrm{II}}\) | `run/test_tts_GII_{pathB,direct}_mscheme.py`, `test_tts_GII_eths_pathB.py` |
| \(\Gamma^{\mathrm{III}_a}\) | `run/test_chi_eta_mscheme.py`, `test_GIIIa_ladder_mscheme.py` |
| \(\Gamma^{\mathrm{III}_b}\) | `run/test_G3b_pathB_fold_mscheme.py`, `test_G3b_pathB_pack_mscheme.py`, `test_G3b_term{1,2}_pathA.py` |
| \(\Gamma^{\mathrm{III}_c}\) | `run/test_tts_GIIIc_mscheme.py`, `test_tts_GIIIc_new_direct.py`, `test_chi_theta_mscheme.py` |
| \(\Gamma^{\mathrm{IV}_a}\) | `run/test_G4a_pathB_mscheme.py`, `test_chi_kappa_*`, `test_tts_GIVa_eths_vs_pathB.py` |
| \(\Gamma^{\mathrm{IV}_b}\) | `run/test_G4b_pathB_fold_mscheme.py`, `test_chi_iota_*`, `test_G4b_factorized_fullZ.py` (λ=0) |
| \(\Gamma^{\mathrm{IV}_c}\) | `run/test_tts_GIVc_mscheme.py`, `test_tts_GIVc_pathB.py`, `test_chi_lambda_mscheme.py` |

Optional pack: `./learn/amc_tts/tensor_pro_final/run_all_gold.sh`
