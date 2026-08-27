# Γ^{III_c} / χ^θ

See [../REDUCED_UNREDUCED.md](../REDUCED_UNREDUCED.md). χ^θ is **T×T→S** (tensor ladder → scalar 2b).

## Packaging (locked 2026-07-29)

Bare product (code ChiTab / AMC numerator):

\[
S^{J_0}_{ijkl}
=
\sum_{ab J_2}
(-1)^{J_0+J_2+\lambda}\,\hat\lambda^{-1}
w\,\Omega^{J_0 J_2\lambda}_{ijab}\Omega^{J_2 J_0\lambda}_{abkl}
\]

| Object | Formula | AMC flag |
|---|---|---|
| \(\chi_{\mathrm{red}}\) | \(S/\hat J_0\) | `reduce=true` |
| \(\chi_{\mathrm{unred}}\) | \(S/\hat J_0^2\) | default |
| m | \(\sum w\,[\Omega\times\Omega]^{(0)}\) CG | — |

**Do not** store bare \(S\) as reduced (old bug: fake \(\hat J\) mismatch vs m).

### AMC was correct

The m ↔ code FAIL was **not** a wrong AMC angular formula. AMC’s `reduce` flag is the rule:

- default print = unreduced (\(\hat J^{-2}\) on \(S\) for χ^θ);
- `reduce=true` = reduced (\(\hat J^{-1}\) on \(S\)).

Code stored bare \(S\) labeled reduced → ratios \(\sim\hat J^{\pm1}\). Fix: store \(\chi_{\mathrm{red}}=S/\hat J\) (DIRECT) / pack to \(\chi_{\mathrm{unred}}=S/\hat J^2\) (Path B). **No AMC retune, no extra hats.**

See [../REDUCED_UNREDUCED.md](../REDUCED_UNREDUCED.md) §“AMC is correct”.

Code:
- DIRECT: store \(\chi_{\mathrm{red}}=S/\hat J\), Case-2 9j, write \(Z_{\mathrm{unred}}=Z_{\mathrm{red}}/\hat J\)
- Path B: pack \(\chi_{\mathrm{red}}\) then `MakeNotReduced` → \(\chi_{\mathrm{unred}}\) for Pandya

## AMC inputs

| File | Role |
|---|---|
| `input/chi_theta.txt` | χ^θ unreduced |
| `input/chi_theta_reduced.txt` | χ^θ `reduce=true` |
| `input/G3c_direct.txt` | unfactored ΩΩΓ with `P(i/j)P(k/l)` |
| `input/G3c_from_chi.txt` | fold unreduced (hats \(\hat J_2^2\hat J_3^2\)) |
| `input/G3c_from_chi_reduced.txt` | fold Case-2 (hats \(\hat J_0\hat J_2\hat J_3\)) |

## Status

| Piece | Status |
|---|---|
| χ^θ m ≡ AMC | **PASS** `run/test_chi_theta_mscheme.py` |
| m ≡ DIRECT ≡ Path B | **PASS** `run/test_tts_GIIIc_mscheme.py` |
| DIRECT ≡ Path A ≡ Path B ≡ AMC 6j/9j | **PASS** `run/test_tts_GIIIc_new_direct.py` |
| ethS Path B DGEMM (χ) + Pandya | **on** (`SetUse_TypeGIIIc_factorized(True)`) |

## Re-run

```bash
amc -o output/chi_theta.tex input/chi_theta.txt
amc -o output/chi_theta_reduced.tex input/chi_theta_reduced.txt
amc --collect-ninejs -o output/G3c_from_chi_ninej.tex input/G3c_from_chi.txt
amc --collect-ninejs -o output/G3c_from_chi_reduced_ninej.tex input/G3c_from_chi_reduced.txt
amc -o output/G3c_direct.tex input/G3c_direct.txt
```
