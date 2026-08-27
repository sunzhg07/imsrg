# Factorized TTS \(\Gamma^{\mathrm{III}_a}\) / \(\chi^\eta\)

All extract code: **`FactorizedDoubleCommutator_eths.cc`**
(`comm223_232_GIIIa`). Living pitfalls: [../LESSONS.md](../LESSONS.md).

## χ^η isolation (locked 2026-07-29)

\[
\chi^{\eta}_{ijkl}
=\sum_{ab}
\bigl(\bar n_a n_b\bar n_k+n_a\bar n_b n_k\bigr)
\,\Omega_{iabl}\,\Omega_{bjka}
\qquad(T\times T\to S).
\]

| Link | Result |
|---|---|
| m ≡ AMC direct | **PASS** ~1e-15 |
| m ≡ Pandya→mid→inv | **PASS** ~1e-15 |
| Path A ≡ Path B | **PASS** |

AMC: `input/chi_eta_{direct,via_pandya}.txt` · Bench: `run/test_chi_eta_mscheme.py`.  
Inv: `chi=+barChi` (drop sample minus). Mid RME = DGEMM step.

**Verdict:** χ^η gold chain fully locked (m ≡ DIRECT ≡ Path B).

## Γ^{III_a} ladder from χ^η (locked 2026-07-29)

\[
\Gamma^{\mathrm{III}_a}_{ijkl}
=-\sum_{ab}\Bigl\{
  (1-\hat P_{ij})\,\chi^\eta_{ijab}\Gamma_{abkl}
 +(1-\hat P_{kl})\,\Gamma_{ijab}\chi^\eta_{klab}
\Bigr\}.
\]

| Link | Result |
|---|---|
| m ≡ AMC `G3a_from_chi` | **PASS** ~1e-14 |
| m ≡ χ^η×Γ channel DGEMM | **PASS** ~1e-14 |
| AMC ≡ DGEMM | **PASS** ~1e-15 |

AMC: `input/G3a_from_chi.txt` (reduce=true → \(Z_{\mathrm{red}}=-\chi_{\mathrm{red}}\Gamma_{\mathrm{unred}}\)).  
Bench: `run/test_GIIIa_ladder_mscheme.py`.  
DGEMM: \(Z=-\mathrm{Chi}_{AS}\Gamma-\Gamma(\mathrm{Chi}_{AS})^{T}\) with
\(\mathrm{Chi}_{AS}[ij,ab]=\chi_{ijab}-(-1)^{J+j_i+j_j}\chi_{jiab}\).

**Gold** = m / from_chi / χ×Γ — **ladder fully locked** (m ≡ AMC ≡ DGEMM).

## Naming

| This doc / AMC | Scalar Factorized label | Topology |
|---|---|---|
| \(\Gamma^{\mathrm{III}_a}\) | often IIa / IIc | Ladder \(\chi\Gamma\) / \(\Gamma\chi\) |
| \(\Gamma^{\mathrm{III}_b}\) | often IIb / IId | Pandya → RC → DGEMM → inv |

\(\chi^\eta\) is always **scalar** (\([\Omega\times\Omega]^0\)). \(\lambda\) is \(\Omega\)’s rank only.

## Path B extract (locked 2026-08-27)

`comm223_232_GIIIa`:

1. **χ̄^η** 2n CC: AMC same-label Pandya(Ω) → occ-weighted DGEMM (`w= n̄_a n_b n̄_k + n_a n̄_b n_k`; λ≠0 mid-J).
2. **Inv Pandya** to ordinary χ_red (AMC: \(\chi=\hat J_0\sum \hat J'\,6j\,\bar\chi\); no \((1-P)\) on χ).
3. Expand Γ to 2n (GetTBME).
4. \(Z_{\mathrm{red}}=-\mathrm{Chi}_{AS}\Gamma-\Gamma\mathrm{Chi}_{AS}^{T}\); store \(Z_{\mathrm{unred}}=Z_{\mathrm{red}}/\hat J\) (√2 for \(i=j\)/\(k=l\)).

Do **not** pack χ into n×n `GetMatrix` (scalar IIa storage trick; drops non-AS pieces).
Do **not** use Factorized IIa Pandya (adcb / `occ_AbarBC`) for this extract — that χ̄ ≠ locked AMC χ.

**Bench:** `run/test_tts_GIIIa.py` vs AMC ladder gold. ethS ≡ ladder **PASS** (~2e-14 @ emax=1 λ=2).

## Files

| Item | Path |
|---|---|
| ethS extract | `src/FactorizedDoubleCommutator_eths.cc` (`comm223_232_GIIIa`) |
| Unfactored gold | AMC `G3a_from_chi` / χ Path A (m-scheme benches) — no `tts_GIIIa` |
| χ^η lock | `run/test_chi_eta_mscheme.py` |
| Ladder lock | `run/test_GIIIa_ladder_mscheme.py` |
| ethS ≡ ladder | `run/test_tts_GIIIa.py` |
| Agent lessons | [../LESSONS.md](../LESSONS.md) §\(\Gamma^{\mathrm{III}_a}\) |

## Next

1. Shared χ^η helper with GIIIb (optional dedupe).
