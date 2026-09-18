# $\Gamma^{\mathrm{IV}_c}$ / $\chi^\lambda$ — Path B

**Production = Path B only** (Pandya → DGEMM → inv). CG leftover removed.

**Read first if debugging again:** [LESSON.md](LESSON.md) — why the long debug,
wrong AMC print, adcb map, missing \((-1)^{J_p}/\hat J_p\).

Gold mid-J is **tts_ring** \(X_{pqsr}=\sum\chi_{p\bar a rb}\Omega_{aqsb}\),
**not** printed `G4c_from_chi*_ninej.tex`.

## Packaging (locked 2026-09-18, λ=3)

| Object | Flag | Meaning |
|---|---|---|
| \(\chi^\lambda\), \(\Omega\) | tensor | **WE-reduced** (ChiTab ≡ `GetTBME_J`; `Eta.IsReduced()`) |
| \(\bar\chi\), \(\bar\Omega\) | AMC Eq1/2 | **IMSRG** `adcb(i,l,k,j)` ≡ AMC `bar(i,j,k,l)` |
| \(\bar G\), \(Z\) | `reduce=true` | \(Z_{\mathrm{red}}=S/\hat J\); ethS stores \(Z_{\mathrm{unred}}=Z_{\mathrm{red}}/\hat J\) |

Mid factor: \(\bar G^{J_p}\mathrel{+}=(-1)^{J_p}/\hat J_p\cdot\hat\lambda^{-1}(-1)^{J_{ab}+\lambda}\,\bar\chi\,\bar\Omega\).
Inv: drop AMC-sample overall minus.

### 2n×2n — normalized vs unnormalized

| API / layout | Convention |
|---|---|
| `GetTBME_J` / ChiTab | **unnormalized** |
| `GetTBME_norm` / `GetMatrix` | **normalized** |
| Pandya 2n pack | **unnormalized**; ÷√2 only on `AddToTBME` |

## Status

| Step | Status |
|---|---|
| Fwd Pandya χ/Ω → 2n (adcb) | **PASS** |
| Mid-J tts_ring ≡ m/CG bare | **PASS** |
| Path B adcb mid+inv ≡ Path A | **PASS** |
| ethS production leftover | Path B only |

Benches: `run/test_givc_midj_bare_vs_cg.py`, `run/test_givc_pandya_2n_lambda3.py`.
