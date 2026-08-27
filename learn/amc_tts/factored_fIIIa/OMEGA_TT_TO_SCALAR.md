# Tensor × tensor → scalar χ (λχ = 0) — AMC Path A / Path B

## Packaging (locked)

| Object | AMC declare |
|---|---|
| \(\Omega\) | `scalar=false` (always reduced) |
| \(\chi\), \(\bar\chi\) | `scalar=true, reduce=true` |

Compare to m via \(X^J_{\mathrm{red}}=\hat J^{-1}\sum_m\mathrm{CG}\,\mathrm{CG}\,X(m)\).
See [REDUCED_UNREDUCED.md](../REDUCED_UNREDUCED.md).

## m-scheme

\[
\chi_{ijkl}=\sum_{ab}\Omega_{ajkb}\,\Omega_{ibal}
\qquad(+\;\chi^\gamma\text{ occ weights when present})
\]

with \([\Omega^\lambda\times\Omega^\lambda]^{(0)}\) CG.

## AMC equations (re-run)

Inputs: `input/omega_cross_noocc_{direct,pandya}.txt`, `input/chi_gamma_{direct,via_pandya}.txt`.

| Path | Role |
|---|---|
| **A (direct)** | 5×6j / W1·W2 on \(\Omega\Omega\) |
| **B (Pandya)** | tensor Pandya → mid RME \(\bar\chi\) → inv Pandya |

### Path B correction (locked by m)

AMC’s **sample** inv uses `chi = - barChi`, which prints an overall minus on eq 3.
Numerically that makes

\[
\chi_{\mathrm{AMC\_B\,(sample)}} = -\,\chi_{\mathrm{m}} = -\,\chi_{\mathrm{AMC\,direct}}.
\]

**Correct Path B:** `chi_ijkl = barChi_ijkl` (drop the sample minus). Then AMC Path B ≡ m ≡ direct.

| Comparison (emax=1, λ=2, reduce=true) | Result |
|---|---|
| m vs AMC direct | **PASS** |
| m vs AMC Path B (sample `chi=-barChi`) | **FAIL** (exactly −1) |
| m vs AMC Path B (corrected `chi=barChi`) | **PASS** |

Same sign bug as the ring (`tts_ring.md`: Path B printed = −Path A).

### Corrected Path B (reduce=true)

1. \(\bar\Omega\) — AMC eq 1 (tensor Pandya, scheme `((1,-4),(3,-2))`)
2. \(\bar\chi^{J}=\displaystyle\frac{(-1)^{J}}{\hat J}\sum_{ab J_2}(-1)^{J_2+\lambda}\hat\lambda^{-1}\,w\,\bar\Omega^{J_2 J}_{ajkb}\,\bar\Omega^{J J_2}_{ibal}\)
3. \(\chi^{J_0}=\hat J_0\sum_{J}\hat J\begin{Bmatrix}j_l&j_k&J_0\\ j_j&j_i&J\end{Bmatrix}\bar\chi^{J}\)
   (**no** overall minus)

## Bench

```bash
PYTHONPATH=build python3 run/test_omega_cross_mscheme.py 2 1   # m ≡ direct
# Path B literal vs m: see run/test_chi_gamma_pathB_amc.py
```

## Takeaway for ethS

- Trust **AMC direct** (Path A) vs m with reduced packaging + full \(j_0\) range.
- Trust **AMC Path B** only after dropping the sample inv minus.
- Do **not** use the old Neithan \(\times(-\hat\lambda(-1)^{J+\lambda})\) story as gold — that was compensating for wrong packaging / wrong inv sign.
