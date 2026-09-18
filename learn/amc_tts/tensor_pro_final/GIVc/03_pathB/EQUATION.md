# $\Gamma^{\mathrm{IV}_c}$ / $\chi^\lambda$ — Path B

Path B: Pandya → DGEMM → inv (`G4c=barG4c`, drop sample minus).

Input: `G4c_chi_omega_pandya.txt`. Printed `G4c_from_chi*_ninej.tex` is **wrong** vs m — do not use as gold.

## Packaging (locked 2026-09-18, λ=3)

| Object | Flag | Meaning |
|---|---|---|
| \(\chi^\lambda\), \(\Omega\) | tensor | **WE-reduced** (ChiTab ≡ `GetTBME_J`; `Eta.IsReduced()`) |
| \(\bar\chi\), \(\bar\Omega\) | AMC Eq1/2 | same-label print; **IMSRG** `adcb(i,l,k,j)` ≡ AMC `bar(i,j,k,l)` |
| \(\bar G\), \(Z\) | `reduce=true` | \(Z_{\mathrm{red}}=S/\hat J\); ethS stores \(Z_{\mathrm{unred}}=Z_{\mathrm{red}}/\hat J\) |

### 2n×2n full storage — normalized vs unnormalized

Do **not** mix these in the Pandya pack:

| API / layout | Convention |
|---|---|
| `GetTBME_J` / ChiTab | **unnormalized** (\(\sqrt2\times\) stored if \(a=b\) or \(c=d\)) |
| `GetTBME_norm` / `GetMatrix` | **normalized** (identical already \(\div\sqrt2\)) |
| Factorized `CHI_VII` / `DoTensorPandya` 2n | pack **unnormalized**; \((p,q)\) and \((q,p)\) separate if \(p\neq q\); \(p=q\) once |
| Occupancy-AS \(\chi^\lambda\) | same unnormalized 2n (no fermionic AS) |
| Write `AddToTBME` | fold is unnormalized → \(\div\sqrt2\) only on store |

Bench: `run/test_givc_pandya_2n_lambda3.py` — fwd Pandya map + 2n×2n CC pack **PASS**.

Do **not** use `PandyaChiIalb` / custom \((ia)(lb)\) helpers for the speed path — they disagree with AMC/adcb.
## Status

| Step | Status |
|---|---|
| 1. AMC Path B Eq1–4 rerun | done (`G4c_chi_omega_pandya*_ninej.tex`) |
| 2. Fwd Pandya χ/Ω → 2n×2n (adcb) | **PASS** |
| 3. Reduced/unreduced named | locked (table above) |
| 4. Mid-J ≡ CG bare | **PASS** — gold is **tts_ring** Path A \(X_{pqsr}=\sum\chi_{p\bar a rb}\Omega_{aqsb}\), not printed `G4c_from_chi` |
| 5. Path B adcb mid+inv ≡ Path A | **PASS** (Python); ethS `kind=1` uses mid \((-1)^{J_p}/\hat J_p\) + corrected inv |

Bench: `run/test_givc_midj_bare_vs_cg.py` (FoldAS=False, T1, λ=3).

