# Γ^II / χ^ζ — m-scheme → AMC → Path B

**Discipline:** strip \(P\) for AMC → restore by hand. Do **not** feed
permuted indices into AMC. Code computes **one bra topology** and **one ket
topology** (each with \((1-P)\)), not four expanded AMC copies.

See also [../factored_GI/NOTES.md](../factored_GI/NOTES.md),
[../REDUCED_UNREDUCED.md](../REDUCED_UNREDUCED.md).

---

## 1. M-scheme (locked)

### χ^ζ (tensor 1b, rank \(\lambda_\Omega\))

\[
\chi^\zeta_{ij}
=
\frac12\sum_{abc}
w\,
\Gamma_{aibc}\,\Omega_{bcaj},
\qquad
w=n_a n_b\bar n_c+\bar n_a\bar n_b n_c.
\]

(Equivalent to analyze \(\Gamma_{ciab}\Omega_{abcj}\) by rename. χ is
**non-Hermitian**.)

### Γ^II — two topologies; sum is Hermitian

\[
\boxed{
\Gamma^{\mathrm{II}}_{ijkl}
=
\sum_a\Bigl[
\underbrace{(1-\hat P_{ij})\,\chi^\zeta_{aj}\,\Omega_{iakl}}_{\text{bra }W}
\;-\;
\underbrace{(1-\hat P_{kl})\,\chi^\zeta_{ak}\,\Omega_{ijal}}_{\text{ket }V}
\Bigr].
}
\]

| Piece | Hermitian alone? |
|---|---|
| Bra \(W=(1-P_{ij})\chi_{aj}\Omega_{iakl}\) | **no** |
| Ket \(V=(1-P_{kl})\chi_{ak}\Omega_{ijal}\) | **no** |
| \(W-V\) | **yes** |

Do **not** replace ket by \(h_Z W^{T}\). That is the `comm122ss` trick when the
1b factor is Hermitian; here χ^ζ is not, so \(W+hW^{T}\neq W-V\).

**AMC orientation note:** analyze writes \(\chi_{aj}\). AMC Yutsis needs
\(\chi_{ja}\Omega_{iakl}\) (swap χ legs). Same physics; angular print uses
that orientation.

---

## 2. AMC DIRECT — strip → expand χ → restore

Unfactored expand of χ into \(\Gamma\,\Omega\,\Omega\), same strip/restore as Path B.

| Topology | Input (P stripped) | Output |
|---|---|---|
| Bra \(W\) | `input/G2_Wdirect_noperm.txt` | `output/G2_Wdirect_noperm.tex` |
| Ket \(V\) | `input/G2_Vdirect_noperm.txt` | `output/G2_Vdirect_noperm.tex` |

Restore \((1-P_{ij})\) / \((1-P_{kl})\) by hand; assemble \(\Gamma^{\mathrm{II}}=W-V\).

AMC tensor overall signs are unreliable — numerical lock is
`run/test_tts_GII_direct_mscheme.py` (DIRECT ≡ Path B).

Optional P-kept expand: `input/G2_direct.txt`.

---

## 3. AMC Path B — strip → couple → restore

### Step A — χ^ζ alone

**Input** `input/chi_zeta.txt` (no \(P\)).

**Output** `output/chi_zeta.tex` — reduced tensor 1b RME.

### Step B — bra topology, \(P\) stripped

**Input** `input/G2_Wbra_noperm.txt`:

```text
W_ijkl = sum_a(chizeta_ja * Omega_iakl);
```

**Output** `output/G2_Wbra_noperm.tex` — **one** RME (unreduced \(Z\), reduced
χ/Ω, \(\lambda=\mathrm{rank}(\Omega)\)):

\[
W^{J}_{ijkl}
=
(-1)^{j_i}\hat J^{-1}
\sum_{a J_2}
(-1)^{J_2+j_a}\,\hat J_2\,\hat\lambda^{-1}
\begin{Bmatrix}J & J_2 & \lambda \\ j_a & j_j & j_i\end{Bmatrix}
\chi^\zeta_{ja}{}^\lambda\,
\Omega^{J_2 J\lambda}_{iakl}.
\]

**Restore by hand** \((1-\hat P_{ij})\) (fermionic bra exchange / `Ket::Phase`):
adds the exchange copy with \(\chi_{ia}\Omega_{jakl}\) and phase
\((-1)^{J+j_i}\) (same SixJ with \(i\leftrightarrow j\)).

### Step C — ket topology, \(P\) stripped

**Input** `input/G2_Wket_noperm.txt`:

```text
V_ijkl = sum_a(chizeta_ak * Omega_ijal);
```

**Output** `output/G2_Wket_noperm.tex` — **one** RME:

\[
V^{J}_{ijkl}
=
(-1)^{j_l}\hat J^{-1}
\sum_{a J_2}
(-1)^{J_2+j_a}\,\hat J_2\,\hat\lambda^{-1}
\begin{Bmatrix}J & J_2 & \lambda \\ j_a & j_k & j_l\end{Bmatrix}
\chi^\zeta_{ak}{}^\lambda\,
\Omega^{J J_2\lambda}_{ijal}.
\]

**Restore by hand** \((1-\hat P_{kl})\): exchange copy with \(\chi_{al}\Omega_{ijak}\)
and phase \((-1)^{J+j_l}\).

### Step D — assemble

\[
\Gamma^{\mathrm{II}\,J}=W_{\mathrm{restored}}-V_{\mathrm{restored}}.
\]

(\(\tfrac12\) lives inside χ^ζ, not here.)

### Optional: AMC with \(P\) expanded

`input/G2_from_chi_RME.txt` still has \(P(i/j)\), \(P(k/l)\) and prints four
terms — useful as a check that restore matches, **not** the implementation
target.

---

## 4. Path B / ethS code pattern

Mirror `comm122ss` **structure** (one insertion loop + channel factors), but
keep **both** topologies:

```text
1. Build Chi_zeta (AMC χ^ζ, any λ)     — reduced tensor 1b
2. W  ← (1-P_ij) χ_ja Ω_iakl           — flipphase / SQRT2 like GI bra
3. V  ← (1-P_kl) χ_ak Ω_ijal           — same on ket
4. Z  += W - V                           — Hermitian as a sum
5. AddToTBME (upper triangle OK)
```

λ=0 is the same formulas (\(\hat\lambda^{-1}=1\), equal-\(J\) SixJ); no separate
cheap scalar path.

---

## 5. Status

| Piece | Status |
|---|---|
| χ^ζ m ≡ AMC | **PASS** `run/test_chi_zeta_mscheme.py` (any λ) |
| m ≡ Path B (\(W-V\)) | **PASS** `run/test_tts_GII_pathB_mscheme.py` (any λ) |
| DIRECT ≡ Path B (\(W-V\)) | **PASS** `run/test_tts_GII_direct_mscheme.py` (any λ) |
| ethS ≡ Path B | **PASS** `run/test_tts_GII_eths_pathB.py` (any λ) |
| \(W+hW^{T}\) ≡ \(W-V\) | **FAIL** (χ non-Hermitian) — do not use |

**Verdict:** \(\Gamma^{\mathrm{II}}\) gold chain fully locked (m ≡ DIRECT ≡ Path B ≡ ethS).

### Re-run AMC

```bash
cd learn/amc_tts/factored_GII
amc -o output/chi_zeta.tex input/chi_zeta.txt
amc -o output/G2_Wbra_noperm.tex input/G2_Wbra_noperm.txt
amc -o output/G2_Wket_noperm.tex input/G2_Wket_noperm.txt
amc -o output/G2_Wdirect_noperm.tex input/G2_Wdirect_noperm.txt
amc -o output/G2_Vdirect_noperm.tex input/G2_Vdirect_noperm.txt
```
