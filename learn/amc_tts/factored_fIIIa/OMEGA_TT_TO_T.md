# Tensor × tensor → tensor χ (AMC)

## m-scheme

\[
\chi_{ijkl}=\sum_{ab}\Omega^{(1)}_{ajkb}\,\Omega^{(2)}_{ibal}
\qquad(\lambda_1,\lambda_2\to\lambda_\chi\neq 0)
\]

## AMC status

| Path | Input | Result |
|---|---|---|
| **Direct** | `input/tt_to_t_direct.txt` | **FAIL** — Yutsis “Bigger than square not implemented yet” |
| **Pandya** | `input/tt_to_t_pandya.txt` | **OK** — `output/tt_to_t_pandya_{plain,ninej}.tex` |

So for TT→T, AMC can only deliver the Pandya pipeline, not a single Ω×Ω direct formula.

## Pandya equations (`--collect-ninejs`)

1. **Forward Pandya** \(\bar\Omega^{(1)}\), \(\bar\Omega^{(2)}\): 1×9j each (same as tensor Pandya).
2. **Product (RME)** — mid-\(J\) + rank 6j:
\[
\bar\chi^{J_0 J_1\lambda}
=(-1)^{J_0+J_1}\,\hat\lambda
\sum_{ab\,J_2\,\lambda_1\lambda_2}
(-1)^{\lambda_1+\lambda_2}
\begin{Bmatrix}\lambda&\lambda_1&\lambda_2\\ J_2&J_0&J_1\end{Bmatrix}
\bar\Omega^{(1)\,J_2 J_1\lambda_1}_{ajkb}\,
\bar\Omega^{(2)\,J_0 J_2\lambda_2}_{ibal}
\]
(AMC prints \(\sj{\lambda}{\lambda_1}{\lambda_2}{J_2}{J_0}{J_1}\).)3. **Inverse Pandya** \(\chi=-\bar\chi\): 1×9j (drop overall minus in numerics if both forwards used `bar=-Ω`).

Plain (no `--collect-ninejs`): Pandya fwd/inv as 3×6j; product still 1×6j (rank).

```bash
amc --collect-ninejs -o output/tt_to_t_pandya_ninej.tex input/tt_to_t_pandya.txt
amc -o output/tt_to_t_pandya_plain.tex input/tt_to_t_pandya.txt
```
