# AMC sources for the factorized tensor double commutator

AMC inputs (`*/input/*.txt`) and their generated J-scheme equations
(`*/output/*.tex`) behind `learn/factorized_code_analyze.tex`, Sec. `jscheme`.
The doc writes out the equations the code implements; this directory is the
provenance for them. Code: `src/FactorizedDoubleCommutator_eths.cc`.

## Regenerate / check

```sh
AMC_PYTHON=/path/to/python learn/factorized_amc/regenerate.sh --check   # diff against stored outputs
AMC_PYTHON=/path/to/python learn/factorized_amc/regenerate.sh           # rewrite outputs
```

`--check` reproduces all 61 stored outputs byte-for-byte (AMC from
`/Users/wolf/Apps/norm`, checked 2026-11-11).

## Conventions (apply to every file)

- Always `--wet-convention wigner` (matches `UnitTest::GetMscheme*`).
- `<stem>_ninej.tex` = same input with `--collect-ninejs`; `leftover132/` and
  `nested_tts/` are always run with `--collect-ninejs`.
- AMC default is `reduce = reduce or not scalar`:

  | Object | AMC declaration | IMSRG storage |
  |---|---|---|
  | tensor Ω, W (1b/2b/3b) | `scalar=false` (always reduced) | WE-reduced RME |
  | scalar Γ | `scalar=true, reduce=false` | `IsReduced()=false` |
  | scalar output Z | `scalar=true, reduce=false` (**must be explicit**) | unreduced |
  | reduced reference (`*_reduced.txt`) | `scalar=true, reduce=true` | Z_red = Ĵ Z_unred |

  Unreduced one-body outputs carry ĵ⁻²; reduced ones ĵ⁻¹.
- AMC strips the permutation operators P and the overall ¼ / ½ Pauli factors;
  the code restores them by hand. `*_noperm.txt` inputs make this explicit.
- `--collect-ninejs` can drop a sign (−1)^{2j₀} when it builds a 9j from three
  6j with a half-integer dummy j₀. For `GIVb/G4b_leftover_pandya` only the
  plain (6j) output is stored; the 9j form used in the code is derived by hand
  in the doc (`eq:js-GIVb-bars`).

## Files per diagram

"direct" = AMC on the unfactored string; "Path B" = intermediate χ then fold
(the form the code implements). Equation labels refer to the doc.

| Dir | Inputs | Doc |
|---|---|---|
| `fI` | `f1` (direct) | `eq:js-chia`, `eq:js-fI` |
| `fII` | `f2a`, `f2b` (direct); `chi_beta`, `f2a_from_chi`, `f2b_from_chi` (Path B) | `eq:js-chib`, `eq:js-fII` |
| `fIIIa` | `f3a`, `chi_gamma_direct` (direct); `chi_gamma`, `chi_gamma_via_pandya`, `fIIIa_from_chi` (Path B) | `eq:js-fIIIa` |
| `fIIIb` | `f3b`, `chi_delta`, `fIIIb_from_chi` | `eq:js-chid`, `eq:js-fIIIb` |
| `GI` | `G1_direct`; `chi_epsilon`, `G1_from_chi` | `eq:js-chie` |
| `GII` | `G2_direct`; `chi_zeta`, `chi_omega_gamma`, `G2_Wbra_noperm`, `G2_Wbra_OG_noperm`, `G2_Wket_noperm` | `eq:js-GII` |
| `GIIIa` | `chi_eta_direct`, `chi_eta_via_pandya`, `G3a_from_chi` | `eq:js-chieta`, `eq:js-GIIIa` |
| `GIIIb` | `G3b_direct(_reduced)`; `IIb_leftover`, `IId_leftover` | `eq:js-GIIIb` (RC is a code primitive, no AMC print) |
| `GIIIc` | `G3c_direct(_reduced)`; `chi_theta(_reduced)`, `G3c_from_chi(_reduced)` | `eq:js-GIIIc` |
| `GIVa` | `chi_kappa_analyze`, `chi_kappa_via_pandya`, `G4a_Wbra_noperm`, `G4a_Wklij_noperm` | Sec. `js-GIVa` |
| `GIVb` | `chi_iota_analyze`, `G4b_W1_leftover`, `G4b_W2_leftover`, `G4b_leftover_pandya` | `eq:js-GIVb` … `eq:js-GIVb-inv` |
| `GIVc` | `chi_lambda`, `chi_lambda_T1`, `chi_lambda_T2`, `G4c_direct(_reduced)`, `G4c_chi_omega_pandya` | Sec. `js-GIVc` |
| `leftover132` | `comm223_132tts_{ladder,onebodyA,onebodyB,onebodyB2,cross}_seed` | `eq:js-132-lad`, Sec. `js-132` |
| `nested_tts` | `comm231tts_unred`, `comm132tts_unred`, `comm232tts_unred`, `comm223tts_unred` | unfactorized reference (`ReferenceImplementations::comm*tts`) |

`GII/output/G2_from_chi_analyze_alt.tex` (the four-term Γ^II fold quoted in the
doc) has no runnable input: it was produced by a separate analysis step, and
`G2_from_chi.txt` fails in AMC ("Inconsistent M projections"), so neither
that input nor `f2.txt` / `chi_eta_RC.txt` (same failure) is kept here.
