# Factorized TTS \(f^{\mathrm{III}_a}\) / \(\chi^\gamma\)

All code: **`FactorizedDoubleCommutator_eths.cc`** (`comm223_231_chi2b_tensor`).

## Ranks / reduce packaging (locked)

| Object | Rank | Reduce | Notes |
|---|---|---|---|
| \(\Omega\) / Eta | \(\lambda\) | **always reduced** | tensors cannot `MakeNotReduced` |
| \(\Gamma\) | \(0\) | **unreduced** in ethS Path B ladder | TTS oracle `MakeReduced(Γ)` instead |
| \(\chi^\gamma\) after inv | \(0\) | written as **reduced**, then **`MakeNotReduced`** | required before Γ fold |
| \(Z\) | \(0\) | unreduced | |

## Path B (production) — stay in Pandya, no inverse (locked 2026-08-26)

**\(\chi^\gamma\) never leaves the Pandya (cross-coupled) representation.** It is
built there, multiplied by \(\bar\Gamma\) there, and traced there. This is
literally the scalar Factorized structure — `FactorizedDoubleCommutator.cc`
`comm223_231_chi2b`, diagram II_a/II_c, Eq. B5c:

```
IntermediateTwobody[ch_cc] = (2 J_cc + 1) * Eta_bar * Eta_bar_nnnn * Gamma_bar
f_pq = ĵ_p^{-2} Σ_{e,J_cc} Intermediate[(pe),(qe)] − partner
```

1. **Forward Pandya(\(\Omega_{\mathrm{red}}\))** — IMSRG legs \((a,d,c,b)\),
   \(\bar\Omega^{J_bJ_k}=-\sum\hat J_1\hat J_2\hat J_b\hat J_k(-1)^{j_b+j_d+J_k+J_2}\,9j\,\Omega^{J_1J_2}\).
   \(\lambda=0\) is the equal-\(J\) limit \(-\sum\hat J_1\hat J\,6j\,\Omega^{J_1}_{\mathrm{red}}\).
2. **Mid-\(J\) DGEMM** — scalar RME, one square block per CC channel:
   \[
   \bar\chi^{\gamma\,J}=\frac{\hat\lambda^{-1}}{2J+1}\sum_{J'}(-1)^{J+J'+\lambda}\,
   \bar\Omega^{JJ'}\,\bigl(w_\gamma\odot\bar\Omega^{J'J}\bigr),
   \]
   \(w_\gamma(a,b;c,d)=n_a\bar n_b\bar n_c n_d-\bar n_a n_b n_c\bar n_d\)
   (= the scalar `occ_factor`, AMC ket legs \(c=k,\;d=j\)).
3. **\(\bar\Gamma\)** — scalar Pandya of **unreduced** \(\Gamma\), same channel.
4. **\((2J+1)\,\bar\chi\,\bar\Gamma\)** DGEMM, then the scalar CC trace.

**No** inverse Pandya, **no** `MakeNotReduced`, **no** dense `chi_tab`, **no**
ordinary-channel ladder. \(\Omega\) stays reduced for every \(\lambda\).

### Why \((2J+1)^{-1}\) in step 2

At \(\lambda=0\), \(\bar\Omega^{JJ}(\Omega_{\mathrm{red}})=\hat J\,\bar\Omega_{\mathrm{scalar}}(\Omega_{\mathrm{unred}})\).
The \(\hat J^{2}\) that produces cancels \((2J+1)^{-1}\), so step 2 collapses to the
scalar `bar_Eta * nnnbar_Eta` **bit-exact**. One code path, all \(\lambda\).
(Do **not** copy GIIIb's per-element `scale = sqrt((2J1+1)/(2Jbra+1))`; it does
not reduce to the scalar convention here.)

| Bench (emax=1,2,3,4) | \(\lambda\) | Result |
|---|---|---|
| CC Path B vs W1/W2 gold | 0…4 | **PASS** ~1e-16 rel |
| speed vs W1/W2 ladder | 2 | ~9× faster @ emax=4 |

Bench: `run/test_tts_fIIIa_pathB_cc.py`. Oracle flag: `SetUse_TypeIIIa_slow(True)`.

### Old (broken) route — do not revive

\(\bar\chi\to\) inverse Pandya \(\to\) `ForceScalarMakeNotReduced` \(\to\) dense
`chi_tab` \(\to\) \(\hat\jmath^{-2}\sum\hat J^{2}(\Gamma_{cpab}\chi_{abcq}-\chi_{pcab}\Gamma_{abqc})\).
The ladder itself is fine (it is the W1/W2 oracle), but routing \(\bar\chi\)
through the inverse introduced the unresolvable `chi = ±barChi` / \(\hat J\)
packaging question, and it has no DGEMM. Its concrete bugs were: no
`scale`/hat correction on the tensor Pandya, mid-\(J\) phase \((-1)^{J'+\lambda}\)
instead of \((-1)^{J+J'+\lambda}\), and a \(\lambda=0\) branch that applied the
unreduced scalar Pandya formula to a **reduced** \(\Omega\)
(`comm223_231_st` un-reduces only \(\Gamma\)).

## Unfactorized / oracle packaging (different route, same physics)

| Path | \(\Omega\) | \(\Gamma\) | \(\chi\) | Comment |
|---|---|---|---|---|
| **TTS** `comm223_231_tts_fIIIa` | reduced | **`MakeReduced`** | inside AMC reduced formula | Case-2 AMC `f3a.tex`; \(Z\) stored unreduced |
| **ethS slow** `use_TypeIIIa_slow` | reduced | unreduced (as passed) | AMC \(W_1W_2\) dense χ → ladder | debug; intended to match TTS |

Do **not** mix: reduced-Γ AMC formula with unreduced Γ (or the reverse) without converting χ.

## Debug flags

| Flag | Role |
|---|---|
| `use_TypeIIIa_1b` | Enable \(f^{\mathrm{III}_a}\) (default = CC Path B) |
| `use_TypeIIIa_slow` | Force AMC \(W_1/W_2\) dense χ + ladder oracle |

## Files

- ethS: `use_TypeIIIa_1b` branch in `comm223_231_chi2b_tensor`
- Scalar template: `FactorizedDoubleCommutator.cc` `comm223_231_chi2b` §II_a/II_c
- Ref: `comm223_231_tts_fIIIa` (= `build_chi_gamma_tab` + `fold_fIIIa_ladder`)
- Bench: `run/test_tts_fIIIa_pathB_cc.py`, `run/test_tts_fIIIa.py`
- AMC: `output/chi_gamma.tex`, `fIIIa_from_chi.tex`
