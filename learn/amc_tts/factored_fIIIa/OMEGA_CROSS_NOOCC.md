# Ω×Ω cross-coupled diagram — no occupation (simplified)

## AMC conventions (locked)

```text
declare Omega { mode=4, scalar=false }              # tensor → always reduced
declare chi   { mode=4, scalar=true, reduce=true }  # scalar reduced
# Path B: barChi also scalar=true, reduce=true
```

## m-scheme

\[
\chi_{ijkl} = \sum_{ab} \Omega_{ajkb}\,\Omega_{ibal}
\]

## AMC outputs

| Path | Input | Output |
|---|---|---|
| A direct | `input/omega_cross_noocc_direct.txt` | `output/omega_cross_noocc_direct_{plain,ninej}.tex` |
| B Pandya | `input/omega_cross_noocc_pandya.txt` | `output/omega_cross_noocc_pandya_{plain,ninej}.tex` |

With `reduce=true` on χ (vs previous unreduced):

| Piece | Unreduced χ | Reduced χ |
|---|---|---|
| Path A body | … | same + overall \(\hat J_0\) |
| \(\bar\chi\) RME | \(\hat J_0^{-2}\) | \(\hat J_0^{-1}\) |
| inv Pandya | \(\sum (2J_p+1)\,6j\,\bar\chi\) | \(\hat J_0\sum \hat J_p\,6j\,\bar\chi\) |

**Sign note:** `barO=-Ω` twice + AMC `χ=-barχ` ⇒ Path B = −Path A at m-scheme.
Numeric Path B **drops** the last minus.

## Benchmark

```bash
python3 run/test_omega_cross_noocc.py 0 1   # λ=0: A≡B PASS
python3 run/test_omega_cross_noocc.py 2 1   # λ≠0: AMC A≠B FAIL (see below)
python3 run/test_tt_to_scalar.py 2 1        # Neithan A≡B PASS
```

## λ≠0 (TT→scalar): use Neithan, not raw AMC

AMC printed Path A/B disagree for λ≠0. The equation that closes is **Neithan**:

\[
\chi_{\mathrm{Neithan}} = -\hat\lambda\,(-1)^{J_0+\lambda}\,\chi_{\mathrm{AMC\_B}}
\]

with Path A = expansion of Neithan Path B. Details: [OMEGA_TT_TO_SCALAR.md](OMEGA_TT_TO_SCALAR.md).
