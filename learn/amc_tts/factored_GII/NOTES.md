# Γ^II / χ^{ΩΓ} + χ^ζ — m-scheme → AMC → Path B

**Discipline:** strip \(P\) for AMC → restore by hand.

Unfactored bra is \(\Omega\times\Gamma\times\Omega\). Ket is \(\Gamma\times\Omega\times\Omega\).
\(\chi^\zeta=\Gamma\times\Omega\) on the **bra** matches only at \(\lambda=0\).

See also [../REDUCED_UNREDUCED.md](../REDUCED_UNREDUCED.md).

---

## 1. M-scheme (locked)

### χ^ζ (ket; \(T\times S\to T\))

\[
\chi^\zeta_{ij}
=
\frac12\sum_{abc}
w\,
\Gamma_{ciab}\,\Omega_{abcj}
=
\frac12\sum_{abc}
w\,
\Gamma_{aibc}\,\Omega_{bcaj},
\qquad
w=n_a n_b\bar n_c+\bar n_a\bar n_b n_c.
\]

### χ^{ΩΓ} (bra; \(T\times S\to T\))

\[
\chi^{\Omega\Gamma}_{ij}
=
\frac12\sum_{abc}
w\,
\Omega_{ciab}\,\Gamma_{abcj}.
\]

Do **not** relabel this to \(\Omega_{aibc}\Gamma_{bcaj}\) while keeping the same
\(w(a,b,c)\) — unlike χ^ζ, that rename is **not** equivalent.

### Γ^II

\[
\boxed{
\Gamma^{\mathrm{II}}_{ijkl}
=
-\,(1-\hat P_{ij})\,\chi^{\Omega\Gamma}_{ja}\,\Omega_{iakl}
\;-\;
(1-\hat P_{kl})\,\chi^\zeta_{ak}\,\Omega_{ijal}.
}
\]

Leftover CG on both folds when \(\lambda\neq 0\).

---

## 2. AMC Path B

### χ^{ΩΓ}

**Input** `input/chi_omega_gamma_analyze.txt`.

**Output** `output/chi_omega_gamma_analyze.tex`:

\[
\chi^{\Omega\Gamma\,\lambda}_{ij}
=
\tfrac12(-1)^{j_j+\lambda}
\sum_{abc J_0 J_1}
w\,(-1)^{J_0+j_c}\,\hat J_0\hat J_1
\begin{Bmatrix}J_0 & J_1 & \lambda \\ j_j & j_i & j_c\end{Bmatrix}
\Omega^{J_0 J_1\lambda}_{ciab}\,
\Gamma^{J_1}_{abcj}.
\]

χ^ζ analyze: `input/chi_zeta_analyze.txt` (Γ_ciab Ω_abcj). Older
`chi_zeta.txt` (Γ_aibc Ω_bcaj, same w labels) is **not** equivalent.

### Folds (same 6j as before)

Wbra with \(\chi^{\Omega\Gamma}_{ja}\Omega_{iakl}\): `input/G2_Wbra_OG_noperm.txt`
(identical 6j to old Wbra). Ket: `input/G2_Wket_noperm.txt`.

Assemble \(\Gamma^{\mathrm{II}}=-W-V\). Restore \((1-P)\) by hand.

---

## 3. ethS

```text
1. Build Chi_OG  (AMC χ^{ΩΓ}) and Chi_zeta
2. W ← (1-P_ij) χOG_ja Ω_iakl
3. V ← (1-P_kl) χζ_ak Ω_ijal
4. Z += −W − V
```

---

## 4. Status

| Piece | Status |
|---|---|
| χ^{ΩΓ} m ≡ AMC J | `run/test_tts_GII_pathB_mscheme.py` |
| χ^ζ m ≡ AMC J | `run/test_chi_zeta_mscheme.py` |
| m ≡ Path B ≡ `Mscheme_fact_GII` | `run/test_tts_GII_pathB_mscheme.py` |
| ethS ≡ Path B | `run/test_tts_GII_eths_pathB.py` |
| nested wick 232 ≡ Σ G* | `run/test_nested_tts_vs_mscheme.py` |

### Re-run AMC

```bash
cd learn/amc_tts/factored_GII
amc -o output/chi_omega_gamma_analyze.tex input/chi_omega_gamma_analyze.txt
amc -o output/G2_Wbra_OG_noperm.tex input/G2_Wbra_OG_noperm.txt
amc -o output/chi_zeta.tex input/chi_zeta.txt
amc -o output/G2_Wket_noperm.tex input/G2_Wket_noperm.txt
```
