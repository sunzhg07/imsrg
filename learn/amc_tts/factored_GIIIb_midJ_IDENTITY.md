# Mid-\(J\) identity: \([\bar\Omega\times\bar\Omega]^0\) as Factorized DGEMM

**Goal.** Turn AMC’s scalar projection
\(\hat\lambda^{-1}(-1)^{J+J'+\lambda}\Omega^{JJ'\lambda}\Omega^{J'J\lambda}\)
into a **Factorized-layout** Pandya product that reduces to
`barCHI_III = bar_Eta * nnnbar_Eta` at \(\lambda=0\).

**Oracle.** Factorized-continuous (match λ=0 Factorized bit-exact). May ≠ TTS at λ≠0;
`*_slow` keeps TTS.

---

## 1. Factorized Pandya (legs \(adcb\))

Scalar Factorized (`FactorizedDoubleCommutator.cc`):

\[
\bar\Omega^{J}_{ab,cd}
=
-\sum_{J'}(2J'+1)
\begin{Bmatrix} j_a & j_b & J \\ j_c & j_d & J' \end{Bmatrix}
\Omega^{J'}_{a\,d\,c\,b}.
\]

**Tensor continuous** (NineJ; same legs; recovers SixJ at \(\lambda=0\), \(J'=J\)):

\[
\begin{aligned}
\bar\Omega^{J J'\lambda}_{ab,cd}
&=
-\sum_{J_1 J_2}
\hat J_1\hat J_2\hat J\hat J'\,
\sqrt{\frac{2J_1+1}{2J+1}}\,
(-1)^{(j_b+j_d)/2+J'+J_2}
\begin{Bmatrix}
j_a & j_d & J_1 \\
j_b & j_c & J_2 \\
J & J' & \lambda
\end{Bmatrix}
\Omega^{J_1 J_2\lambda}_{a\,d\,c\,b}.
\end{aligned}
\]

Layout: Factorized **\(2\times n_{\mathrm{Kets}}\)** exchange blocks (never ph-only TensorCommutators mats).

Occupancy for \(\chi^\eta\): Factorized `occ_AbarBC` / partners on the **right** factor (`nnnbar_Eta`).

---

## 2. Scalar projection in Pandya space (the identity)

AMC (ordinary basis), e.g. \(\chi^\alpha\):

\[
\chi^{0}
\sim
\hat\lambda^{-1}
\sum_{JJ'}
(-1)^{J+J'+\lambda}
\Omega^{JJ'\lambda}\Omega^{J'J\lambda}.
\]

**Claim (Factorized continuous).** With the Pandya of §1,

\[
\boxed{
\bar\chi^{\eta\,0\,J}
=
\hat\lambda^{-1}
\sum_{J'}
(-1)^{J+J'+\lambda}\,
\bar\Omega^{J J'\lambda}
\cdot
\bigl(\mathrm{occ}\odot\bar\Omega^{J' J\lambda}\bigr)
}
\]

as a **matrix product** in the \(2n_J\) Factorized layout.

| \(\lambda\) | Content of the sum |
|---|---|
| \(0\) | Only \(J'=J\), weight \(1\) → `bar_Eta * nnnbar_Eta` |
| \(\neq 0\) | Sum mid-\(J'\) with \(\lvert J-J'\rvert\le\lambda\le J+J'\) |

**Why this isolates \(\Lambda=0\):** the weight \(\hat\lambda^{-1}(-1)^{J+J'+\lambda}\) is AMC’s reduced form of \(\sum_\mu T_{\lambda\mu}^* S_{\lambda\mu}\). A naked `Ā*B̄` without that weight mixes \(\Lambda=0,\ldots,2\lambda\).

**Inverse:** \(\bar\chi^{\eta\,0}\) is **scalar** (equal-\(J\) blocks only). Use Factorized **scalar** SixJ inverse Pandya — **no** extra three-6j “tensor reverse” factor \(\mathcal R\).

---

## 3. Why the old ethS mid-\(J\) attempt gave norm \(0\)

| Bug | Detail |
|---|---|
| Wrong Pandya legs | `barred_eta_tensor` used \(\Omega_{abjk}\), not Factorized \(\Omega_{adcb}\) |
| Empty `barCHI_III` | `zeros()` on default-constructed \(0\times 0\) mat (no allocate) |
| Fake “tensor inv” | Audit claimed \(\mathcal R\); code still used scalar inv (correct *if* \(\bar\chi\) is scalar — but \(\bar\chi\) was wrong) |
| GIIIb gate | λ≠0 always TTS; never exercised a fixed mid-\(J\) path |

---

## 4. Downstream (unchanged once \(\bar\chi^{0}\) is right)

\[
\bar\chi
\;\xrightarrow{\mathrm{RC\ (SixJ)}}\;
\bar\chi_{\mathrm{RC}}
\;\xrightarrow{\bar\Gamma\cdot(\,\cdot\,)}\;
\bar C
\;\xrightarrow{\mathrm{inv\ Pandya}}\;
Z.
\]

Same as scalar Factorized IIb / ethS GIIIb λ=0.

---

## 5. Checklist for implementation

1. ~~Rectangular Factorized NineJ Pandya of Ω (GIVa Path B helper style).~~
2. ~~`barCHI_III[J] = Σ_{J'} λ̂^{-1}(-1)^{J+J'+λ} Ω̄[J,J'] * nnnbar[J',J]`.~~
3. ~~Allocate `2n×2n` before fill.~~
4. ~~Scalar RC + inv only.~~
5. ~~Gate: `*_slow`→TTS; else Factorized (all λ).~~
6. ~~λ=0 regression vs current Fact; λ≠0 smoke (Fact≠TTS OK until optionally calibrated).~~

**Validated (emax=1 He4 RandomOp):** mid-J@λ=0 recovers equal-J Factorized bit-exact
(`TwoBodyNorm = 1.2957759368332740e+03`); λ=2 Fac nontrivial (`~6.27e+02`, ≠ TTS OK).
