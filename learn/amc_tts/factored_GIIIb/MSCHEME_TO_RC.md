# Γ^{III_b}: m-scheme → term1 Path A (straight path)

Trusted source only: [`learn/factorized_code_analyze.tex`](../../factorized_code_analyze.tex).
Ignore arxiv.

## 1. Unfactored

\[
\Gamma^{\mathrm{III}_b}_{ijkl}
=-\sum_{abcd}(\bar n_b n_c n_d+n_b\bar n_c\bar n_d)\,
(1-\hat P_{ij})(1-\hat P_{kl})
\bigl(
\Omega_{dcbk}\Omega_{biac}\Gamma_{jald}
+\Omega_{jcbd}\Omega_{balc}\Gamma_{diak}
\bigr).
\]

## 2. Factorized m-scheme

\[
\chi^\eta_{ijkl}
=\sum_{ab}(\bar n_a n_b\bar n_k+n_a\bar n_b n_k)\,
\Omega_{iabl}\Omega_{bjka}.
\]

(χ is scalar, packaged **reduced** because Ω is reduced tensor: \(\chi_{\mathrm{red}}=S/\hat J\).)

Bare factorized fold (P stripped), locked any λ in
`run/test_G3b_pathB_fold_mscheme.py`:

\[
W=-\sum_{ab}\bigl(
  \chi^\eta_{bkai}\,\Gamma_{jbla}
 +\chi^\eta_{lajb}\,\Gamma_{aibk}
\bigr),\qquad
Z=(1-\hat P_{ij})(1-\hat P_{kl})W.
\]

**Term1** (bare Inv locks): \(W_1=-\sum\chi_{bkai}\Gamma_{jbla}\).

**Term2** (needs \((1-P)^2\)): \(W_2=-\sum\chi_{lajb}\Gamma_{aibk}\).

## 3. Why RC on χ (Pandya on Γ)

| Object | Hermitian / AS? | Transform |
|---|---|---|
| \(\chi^\eta\) (occ-weighted) | no | **RC** (code JT primitive) |
| \(\Gamma\) | yes | Fac scalar **Pandya** |

Plain Pandya of \(\chi_{bkai}\) does **not** share the \((b,\bar a)\)
channel with Pandya of \(\Gamma_{jbla}\). Code RC puts χ on the shared
ph channel for DGEMM.

## 4. Locked term1 dual (any λ_Ω) — Path B DGEMM

χ̄ product is **DGEMM** (no Σ_ab orbital loop):

| λ | χ̄ construction |
|---|---|
| 0 | Fac Pandya Ω̄ @ (occ⊙Ω̄) equal-J |
| ≠0 | AMC Pandya mid-J: \(P=\sum_{J_2}(-1)^{J_2+\lambda}/\hat\lambda\, L@R\), then \(\mathrm{CHI}[a,b;c,d]=(-1)^J/\hat J^2\, P[(a,b),(d,c)]\) ≡ mid\((a,d,c,b)/\hat J\) |

Then for all λ:

3. \(\mathrm{RC}_1\) = code RC of \(\bar\chi_{bc,ad}\) only (term1; no pack)
4. \(\bar\Gamma\) = Fac scalar Pandya(Γ)
5. \(\bar W_1 = \bar\Gamma\cdot\mathrm{RC}_1\) (DGEMM)
6. \(Z_1\) = InvPandya_noperm (\(j\bar l,\,i\bar k\to ijkl\); no \((1-P)^2\))

\[
\mathrm{WE}(Z_1)\equiv W_{1,m}\qquad(\lambda_\Omega=0,\ldots,4).
\]

Bench: `run/test_G3b_term1_pathA.py`.

Do **not** use Fac Pandya of normal χ, nor AMC RC[χ]×RC[Γ]→InvRC (product weight still open).

## 5. Locked term2 dual (any λ_Ω) — pack split; bare is \(jilk\)

Same Path B χ̄ DGEMM as term1. Term2 is the other RC pack piece
(\(ad\leftrightarrow bc\)), **not** a Γ transpose:

3. \(\mathrm{RC}_2\) = code RC of \(\bar\chi_{ad,bc}\) only
4. \(\bar W_2 = \bar\Gamma\cdot\mathrm{RC}_2\) (same \(\bar\Gamma\) as term1)
5. \(Z_2\) = InvPandya_noperm (\(j\bar l,\,i\bar k\to ijkl\))

**Bare identity (checked λ=0 and λ=2):**
\[
\mathrm{WE}(Z_2)_{ijkl}\equiv W_{2,m}(jilk)
\qquad\bigl(\text{also }\equiv W_{1,m}(klij)\bigr).
\]
So Inv of the `adbc` leg writes the fold's second string at the
**fully exchanged** slot \(ij\leftrightarrow kl\) composed with \(i\leftrightarrow j\),
\(k\leftrightarrow l\) — i.e. \(ijkl\mapsto jilk\) — not at \(ijkl\).
That is why bare \(\mathrm{WE}(Z_2)\not\equiv W_{2,m}(ijkl)\) while
\((1-P)^2\) still locks:
\[
(1-P)^2 Z_2(ijkl)=(1-P)^2 W_2(ijkl)
\]
because \(Z_2(ijkl)=W_2(jilk)\) maps the four \(P\)-images onto each other.

Gold used in the bench:
\[
(1-P)^2\,\mathrm{WE}(Z_2)\equiv(1-P)^2 W_{2,m}
\qquad(\lambda_\Omega=0,\ldots,4).
\]

Sanity: \(\mathrm{RC}_1+\mathrm{RC}_2\equiv\mathrm{RC}(\mathrm{pack})\);
\((1-P)^2(Z_1+Z_2)\equiv(1-P)^2(W_1+W_2)\).

Bench: `run/test_G3b_term2_pathA.py`.

## 6. Full Path B pack (any λ_Ω)

Production one-GEMM path, locked against fold:

\[
\mathrm{RC}=\mathrm{RC}[\bar\chi_{bc,ad}+\bar\chi_{ad,bc}],\quad
\bar W=\bar\Gamma\cdot\mathrm{RC},\quad
Z=\mathrm{Inv}\to(1-P)^2,
\]
\[
\mathrm{WE}(Z)\equiv Z_{\mathrm{fold}}\equiv Z_{\text{4-index}}
\qquad(\lambda_\Omega=0,\ldots,4).
\]

Bench: `run/test_G3b_pathB_pack_mscheme.py`.
Fold≡4-index: `run/test_G3b_pathB_fold_mscheme.py`.

## 7. What we are not doing (deferred)

- AMC Pandya→CC print ≡ code RC (code RC is the JT primitive)
- AMC RC×RC→InvRC angular product (InvRC∘RC=id OK; product not locked)
- C++ ethS production wiring (bench dual is locked)

## 8. Why RC×RC was abandoned for this step

`InvRC∘RC=id` on χ and Γ PASSes at λ=2, but
\(\mathrm{InvRC}[D]_{lkji}\equiv\mathrm{extract}(W_1)\) with
SixJ\((0,0,0)\) / AMC 04c guesses does **not**. The production dual
(Pandya Γ × code RC χ̄) locks once χ̄ is Path B mid with adcb/\(\hat J\).
