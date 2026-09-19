# $\Gamma^{\mathrm{IV}_c}$ / $\chi^\lambda$ — Path B

**Production = Path B only** (Pandya → DGEMM → inv). CG leftover removed.

**Read first if debugging again:**
- [LESSON.md](LESSON.md) — traps, adcb map, mid factor
- [PARITY_CHANNELS.md](PARITY_CHANNELS.md) — **odd-π** gold vs ring natural
  feed into the same opp-π CC slots (equation rethink, not mid scratch)

Gold mid-J is **tts_ring** \(X_{pqsr}=\sum\chi_{p\bar a rb}\Omega_{aqsb}\)
**only when χ is accidentally AS** (even π). Odd-π T1 needs gold Pandya
\(\chi_{ialb}\Omega_{bjak}\) into those same CC slots.

## Packaging (locked 2026-09-18, λ=3)

| Object | Flag | Meaning |
|---|---|---|
| \(\chi^\lambda\), \(\Omega\) | tensor | **WE-reduced** (ChiTab ≡ `GetTBME_J`; `Eta.IsReduced()`) |
| \(\bar\chi\), \(\bar\Omega\) (ring) | AMC Eq1/2 | **IMSRG** `adcb(i,l,k,j)` ≡ AMC `bar(i,j,k,l)` |
| \(\bar\chi\), \(\bar\Omega\) (gold) | Path B odd | \(\mathrm{Pandya}_{(ia)(lb)\to(il)(ab)}\), \(\mathrm{Pandya}_{(bj)(ak)\to(ab)(kj)}\) |
| \(\bar G\), \(Z\) | `reduce=true` | \(Z_{\mathrm{red}}=S/\hat J\); ethS stores \(Z_{\mathrm{unred}}=Z_{\mathrm{red}}/\hat J\) |

Mid factor (ring, locked even π): \(\bar G^{J_p}\mathrel{+}=(-1)^{J_p}/\hat J_p\cdot\hat\lambda^{-1}(-1)^{J_{ab}+\lambda}\,\bar\chi\,\bar\Omega\).
Inv: drop AMC-sample overall minus.

CC filter (both π): \((\pi_b+\pi_k)\bmod 2=\pi_\eta\) — same-π even, opp-π odd.

### 2n×2n — normalized vs unnormalized

| API / layout | Convention |
|---|---|
| `GetTBME_J` / ChiTab | **unnormalized** |
| `GetTBME_norm` / `GetMatrix` | **normalized** |
| Pandya 2n pack | **unnormalized**; ÷√2 only on `AddToTBME` |

## Status

| Step | Status |
|---|---|
| Fwd Pandya χ/Ω → 2n (adcb / ring) | **PASS** (even π) |
| Mid-J tts_ring ≡ m/CG bare | **PASS** (even π; odd T1 FAIL) |
| Path B adcb mid+inv ≡ Path A | **PASS** (even π) |
| Odd-π gold Pandya into opp-π CC | **PASS** — AMC scheme `((1,-3),(2,-4))` / `((3,-1),(4,-2))` |
| ethS production leftover | Path B gold fills (all π); ≡ ring on even |

Benches: `run/test_givc_midj_bare_vs_cg.py`, `run/test_givc_pandya_2n_lambda3.py`.
