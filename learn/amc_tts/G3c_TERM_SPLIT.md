# G3c / Γ^{III_c} — m-scheme split → AMC → two DIRECT terms

## M-scheme (before AMC)

\[
\Gamma^{\mathrm{III}_c}_{ijkl}
=-\tfrac12\sum_{abcd}f(a,b,c)\,(1-P_{ij})(1-P_{kl})
\bigl(T_1+T_2\bigr)
\]

\[
f=\bar n_a\bar n_b n_c+n_a n_b\bar n_c
\]

| Term | M-scheme seed | Fact strip |
|---|---|---|
| **Term1** | \(\Omega_{abcl}\,\Omega_{idab}\,\Gamma_{cjkd}\) | \(\chi_k\) |
| **Term2** | \(\Omega_{icab}\,\Omega_{abdl}\,\Gamma_{djkc}\) | \(\chi_j\) |

AMC input (split):

```none
# Term1 only
G3c_T1_ijkl = -1/2 * sum_abcd(
    (nbar_a*nbar_b*n_c + n_a*n_b*nbar_c)
    * P(i/j) * P(k/l)
    * Omega_abcl*Omega_idab*Gamma_cjkd);

# Term2 only
G3c_T2_ijkl = -1/2 * sum_abcd(
    (nbar_a*nbar_b*n_c + n_a*n_b*nbar_c)
    * P(i/j) * P(k/l)
    * Omega_icab*Omega_abdl*Gamma_djkc);
```

AMC J-scheme output: `learn/amc_tts/output/G3c.tex` (Term1 / Term2 blocks).

## Code

| API | Meaning |
|---|---|
| `comm223_232_tts_GIIIc(Eta,Gamma,Z, which_term=0)` | DIRECT; `0` both, `1` Term1, `2` Term2 |
| `comm223_232_tts_GIIIc_term1/2` | wrappers |
| `SetUse_TypeGIIIc_which_term(t)` | Path A fact strip select |
| `DebugCompareGIIIcTerms(Eta,Gamma)` | print DIRECT vs fact per term |

## Bench

```bash
PYTHONPATH=build python3 run/test_tts_GIIIc_terms.py
```
