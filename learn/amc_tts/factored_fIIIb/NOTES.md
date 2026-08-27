# Factorized TTS \(f^{\mathrm{III}_b}\) / \(\chi^\delta\)

## Storage (critical)

| Object | Type | Storage rule |
|---|---|---|
| \(\Omega\) | tensor, usually anti-Hermitian | `ch_bra ≤ ch_ket` + tensor phase |
| \(\chi^\delta\) | scalar, **non-Hermitian** (occ) | full \(\chi(ch,ch)\) square; **no** hermiticity fill |
| \(\Gamma\), \(Z\) | scalar Hermitian | usual `ch_bra ≤ ch_ket` |

Do **not** `.t()` tensor \(\Omega\) blocks. Do **not** infer \(\chi_{kl,ij}\) from \(\chi_{ij,kl}\).

## Build \(\chi^\delta\) (RME + DGEMM)

Same-channel only, full `ibra,iket` (non-Hermitian). Production path:

\[
T_{\mathrm{pp/hh}}
=
\sum_{J_2}
(-1)^{J+J_2+\lambda}\hat\lambda^{-1}
\,\Omega^{J J_2}\,W_{\mathrm{pp/hh}}\,\Omega^{J_2 J},
\qquad
\chi
=
\bar n_k\bar n_l\,T_{\mathrm{hh}}
-
n_k n_l\,T_{\mathrm{pp}}.
\]

Store **normalized** matrix entries (`GetMatrix` convention). No \(\tfrac14\) in \(\chi\) (applied as \(4\times 0.25\) with \(\Gamma\)).

**Case 2 storage:** \(\Omega\) tensor reduced; \(\Gamma\) unreduced in ethS fold \(\equiv G_{\mathrm{red}}/\hat J\) in TTS; \(Z\) unreduced (\(\div\hat\jmath^2\)).

## Contract (scalar `Chi_222_b` structure)

```
M = χ(ch_bra) * Γ(ch_bra, ch_ket)
if ch_bra == ch_ket:  M += M.t()          # A+Aᵀ, square OK
else:                 M -= Γ * χ(ch_ket)  # explicit partner; no .t() into block
M *= 2
Z_pq += 0.25 * sum_c M(c,p;c,q) / (2jp+1)
```

`*=2` (not scalar’s `*=4*(2J+1)`): \(A+A^{T}\) already brings the second topology; TTS unfactored has no net \(\hat J^{2}\). Overall \(\tfrac14\) via \(0.25\).
## Files

- ethS: `comm223_231_chi2b_tensor` (`use_TypeIII_1b`)
- Reference: `comm223_231_tts_fIIIb`
- Benchmark: `run/test_tts_fIIIb.py`
- AMC: `output/chi_delta.tex`, `fIIIb_from_chi.tex`
