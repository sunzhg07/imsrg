# One-body finals: implementation procedure

**Living document** for \(f^{\mathrm{I}}\), \(f^{\mathrm{II}}\), \(f^{\mathrm{III}_a}\), \(f^{\mathrm{III}_b}\).

This is the procedure that reflects the locked verification chain and the
optimization ladder. Per-diagram equations live in `factored_*/NOTES.md`;
flags/tests in [FACTORIZED_TTS_IMPLEMENTED.md](FACTORIZED_TTS_IMPLEMENTED.md);
agent pitfalls in [LESSONS.md](LESSONS.md).

**Code:** `src/FactorizedDoubleCommutator_eths.cc`  
**TTS refs:** `ReferenceImplementations::comm223_231_tts_*`

---

## Universal chain (do not skip steps)

```text
m-scheme gold  →  AMC DIRECT (Path A)  →  AMC Path B (factorized)
                                              ↓
                                    Pandya / RME / DGEMM / ladder
                                              ↓
                                    ethS production path
```

| Step | What | Pass criterion |
|---|---|---|
| **0. Packaging** | Name reduced vs unreduced on every object | [REDUCED_UNREDUCED.md](REDUCED_UNREDUCED.md) |
| **1. m-scheme** | Literal unfactored MEs + \([\Omega^\lambda\times\Omega^\lambda]^{(0)}\) when \(T\times T\to S\) | Physical `GetMschemeMatrixElement_*` |
| **2. AMC DIRECT** | AMC on the unfactored / normal-scheme string | `m ≡ DIRECT` (same packaging) |
| **3. AMC Path B** | Pandya → mid product → inv (AMC printed eqs) | `DIRECT ≡ Path B` |
| **4. Code TTS** | `comm223_231_tts_*` = AMC DIRECT literally | `TTS ≡ DIRECT` |
| **5. ethS slow** | Dense AMC / Path A table (`*_slow`) | `slow ≡ TTS` |
| **6. ethS fast** | Pandya + DGEMM / Factorized layout | `fast ≡ slow` (or dual-oracle flag) |

**Do not** jump to DGEMM before step 3 is locked.  
**Do not** use raw Neithan as gold for \(\chi^\gamma\) / \(f^{\mathrm{III}_a}\) — Neithan ≠ m/AMC (~λ̂); see [factored_fIIIa/OMEGA_TT_TO_SCALAR.md](factored_fIIIa/OMEGA_TT_TO_SCALAR.md).

---

## Status summary

| Diagram | χ | m ≡ DIRECT | DIRECT ≡ Path B | ethS fast | Notes / detail |
|---|---|---|---|---|---|
| \(f^{\mathrm{I}}\) | \(\chi^\alpha\) scalar 1b | **PASS** | **PASS** | DGEMM / intermediates | [FACTORIZED_TTS_IMPLEMENTED.md](FACTORIZED_TTS_IMPLEMENTED.md) §f^I |
| \(f^{\mathrm{II}}\) | \(\chi^\beta\) tensor 1b | **PASS** | **PASS** | ethS TypeII | [factored_fII/NOTES.md](factored_fII/NOTES.md) |
| \(f^{\mathrm{III}_b}\) | \(\chi^\delta\) scalar 2b NH | **PASS** | RME+DGEMM same-ch | ethS TypeIII | [factored_fIIIb/NOTES.md](factored_fIIIb/NOTES.md) |
| \(f^{\mathrm{III}_a}\) | \(\chi^\gamma\) scalar 2b NH | **PASS** (χ+ladder) | **PASS** (AMC Path B, drop sample inv minus) | **PASS** — CC Pandya+DGEMM, χ̄ folded with \(\bar\Gamma\), **no inverse Pandya** (λ=0…4) | [factored_fIIIa/NOTES.md](factored_fIIIa/NOTES.md), [OMEGA_TT_TO_SCALAR.md](factored_fIIIa/OMEGA_TT_TO_SCALAR.md) |

---

## Per-diagram procedure

### \(f^{\mathrm{I}}\) / \(\chi^\alpha\)

| Item | Detail |
|---|---|
| m | \(\chi^\alpha\sim\sum\Omega\Omega\) (scalar 1b); \(f\sim\sum\Gamma\chi\) |
| AMC | `learn/amc_tts/input/` / reduced partners as used in TTS |
| DIRECT | `comm223_231_tts_fI` |
| Path B / fast | ethS `use_TypeI_1b` — 1b intermediates / DGEMM where wired |
| Flags | `use_TypeI_1b` |
| Benches | `run/test_tts_fI.py`, `run/test_tts_f_mscheme.py` (shared) |

### \(f^{\mathrm{II}}\) / \(\chi^\beta\)

| Item | Detail |
|---|---|
| m | \(\chi^\beta\sim\sum\Gamma\Omega\) (tensor 1b rank λ); \(f\sim\sum\chi\Omega\) with \(T\times T\to S\) |
| AMC | `factored_fII/input/chi_beta.txt`, `f2a/b_from_chi.txt` |
| DIRECT | `comm223_231_tts_fII` |
| Path B / fast | ethS `use_TypeII_1b` |
| Flags | `use_TypeII_1b` |
| Benches | `run/test_tts_fII.py`, `run/test_chi_beta_mscheme.py` |
| Pitfall | Do not double-count \((-1)^{j_i}\) (outer phase once) |

### \(f^{\mathrm{III}_b}\) / \(\chi^\delta\)

| Item | Detail |
|---|---|
| m | \(\chi^\delta\sim\sum\Omega\Omega\) (scalar 2b, **non-Hermitian** occ) |
| AMC | `factored_fIIIb/input/` |
| DIRECT | `comm223_231_tts_fIIIb` |
| Path B / fast | Same-channel RME: \(\bar\chi\sim\hat\lambda^{-1}(-1)^{J+J'+\lambda}\Omega\,W\,\Omega\) → DGEMM → \(\Gamma\) fold like scalar `Chi_222_b` |
| Flags | `use_TypeIII_1b` |
| Benches | `run/test_tts_fIIIb.py` |
| Pitfall | Full NH fill; no hermiticity conjugate; no `.t()` on tensor Ω blocks |

### \(f^{\mathrm{III}_a}\) / \(\chi^\gamma\)

| Item | Detail |
|---|---|
| m | \(\chi^\gamma_{ijkl}=\sum_{ab}w\,\Omega_{ajkb}\Omega_{ibal}\) + \([\Omega\times\Omega]^{(0)}\); then ladder \(\hat\jmath^{-2}\sum\hat J^{2}(\Gamma\chi-\chi\Gamma)\) |
| AMC DIRECT | `factored_fIIIa/input/chi_gamma_direct.txt` (`reduce=true`) → 5×6j / W1·W2 |
| AMC Path B | `chi_gamma_via_pandya.txt`: Pandya → mid \(\bar\chi\) → inv; **`chi=barChi`** (drop AMC-sample minus) |
| DIRECT code | `comm223_231_tts_fIIIa` = AMC W1/W2 + ladder (full \(j_0\) range) |
| ethS slow | `use_TypeIIIa_slow` → same W1/W2 |
| ethS fast | CC Pandya(Ω_red) + mid-\(J\) DGEMM χ̄ + \((2J+1)\bar\chi\bar\Gamma\) + CC trace. **No inverse Pandya / no `MakeNotReduced`** — the scalar II_a/II_c structure |
| Flags | `use_TypeIIIa_1b`, `use_TypeIIIa_slow` |
| Benches | `run/test_chi_gamma_pathB_amc.py` (m≡DIRECT≡Path B χ), `run/test_tts_fIIIa_mscheme.py`, `run/test_tts_fIIIa.py`, `run/test_tts_fIIIa_pathB_cc.py` (CC≡gold, λ=0…4) |
| **Not gold** | `neithan.tex` TT→0 (≠ m/AMC by ~λ̂) |
| Ring cousin | [factored_fIIIa/tts_ring.md](factored_fIIIa/tts_ring.md) — Path A correct; sample Path B = −A |

---

## Optimization ladder (after equations match)

1. **Dense AMC / TTS** — correctness oracle (`*_slow` or `tts_*`).
2. **W1/W2 factorization** — pull angular sums that depend on one Ω only (fIIIa slow).
3. **Pandya + CC matrices** — move product to ph / CC layout.
4. **DGEMM mid product** — \(\bar\chi = \bar\Omega\,(\mathrm{occ})\,\bar\Omega\) or mid-\(J\) rectangular product.
5. **Inv Pandya + scalar ladder** — reuse Factorized / IMSRG inv where packaging matches.
6. **Do not** edit scalar `FactorizedDoubleCommutator.cc` for TTS experiments — keep ports in ethS.

---

## Quick test map

| Goal | Command |
|---|---|
| χ^γ m ≡ AMC direct ≡ Path B | `PYTHONPATH=build python3 run/test_chi_gamma_pathB_amc.py 2 1` |
| Cross no-occ m ≡ direct | `PYTHONPATH=build python3 run/test_omega_cross_mscheme.py 2 1` |
| Ring m ≡ Path A | `PYTHONPATH=build python3 run/test_z_ring_mscheme_sign.py 2 1` |
| f^I / f^II / f^IIIb ethS vs TTS | `run/test_tts_fI.py`, `test_tts_fII.py`, `test_tts_fIIIb.py` |
| f^IIIa ethS vs TTS | `run/test_tts_fIIIa.py` (set `*_slow` as needed) |
| f^IIIa m vs DIRECT | `run/test_tts_fIIIa_mscheme.py` |

---

## Related docs

| Doc | Role |
|---|---|
| [LESSONS.md](LESSONS.md) | Agent playbook + pitfalls |
| [REDUCED_UNREDUCED.md](REDUCED_UNREDUCED.md) | Packaging before any m↔J compare |
| [FACTORIZED_TTS_IMPLEMENTED.md](FACTORIZED_TTS_IMPLEMENTED.md) | Checklist of all ethS pieces (incl. Γ) |
| [factored_DGEMM_STATUS.md](factored_DGEMM_STATUS.md) | Γ DGEMM vs Path A status |
| `factored_fII/NOTES.md` etc. | Per-diagram equations |
