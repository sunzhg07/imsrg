# EOM N-kernel: m-scheme → scalar AMC → tensor AMC

`ConstructNormMatrix` is a bilinear kernel: `ComputeNorm(Op1,Op2) = v1^T N v2`.
Each diagram is a Wick string of χ† ρ χ contracted to a **0-body**. There is no Pandya;
J-scheme hats/6j/9j come from recoupling only.

## How to read the code

1. Strip angular-momentum factors from `EOM::ConstructNormMatrix` → m-scheme
   (keep √2 identical-particle norms and ket-exchange `Phase`; those are not recoupling).
2. AMC that string with `scalar=true` → index structure should match the scalar
   loops. **Hat powers often do not** (ρ packaging: code 1b uses ĵ, 2b uses Ĵ;
   AMC identity is ĵ² / Ĵ²). See `learn/amc_tts/REDUCED_UNREDUCED.md`.
3. Same string with `X,Y` tensors and **Y index-reversed** so λ×λ→0
   (`X_ai Y_ia`, `X_abij Y_ijab`). `N` is 0-body `reduce=false`.
4. `ConstructNormMatrix_tensor` implements the tensor AMC formulas, folding
   reversed Y into **stored** MEs (`Flatten` reads `GetOneBody(p,q)` / `GetTwoBody`
   with `ch_bra ≤ ch_ket`).

Regenerate: `learn/amc_tts/eom_norm/regenerate.sh`

## Tensor AMC gotcha

`X_ai Y_ai` fails (“Inconsistent M projections”). Reverse Y: `X_ai Y_ia`.
Same pattern as `comm110tts` / `comm220tts`.

## Stored ME ↔ AMC Y (IMSRG tensor flip)

Tensors are stored as reduced RMEs. Hermitian flip (`Operator::SetOneBody`,
`GetTwoBody` when `ch_bra > ch_ket`):

- 1b: \(Y_{ia}=(-1)^{j_a-j_i} Y_{ai}\)
- 2b: \(Y^{J_1 J_0}_{ijab}=(-1)^{J_1-J_0} Y^{J_0 J_1}_{abij}\)

The kernel multiplies stored vector components, so these phases sit in `N`.

Overall-sign caveat (`learn/amc_tts/NOTES.md`): take 6j / phase / λ̂ from AMC;
λ=0 tensor N will **not** numerically match scalar `ConstructNormMatrix`
(reduced vs unreduced hats, and AMC’s extra minus from λ×λ→0).

## Diagrams

| Label | Block | m-scheme (hats off) | Scalar code hats | Tensor AMC (unreduced 0-body) |
|---|---|---|---|---|
| B4 | qv-qv | \(X_{ai}\rho_{ij}Y_{aj}\) same q | \(\rho\,\hat j_v\) | \(-\,(-1)^{j_a+j_i+\lambda}\hat\lambda^{-1} X_{ai}\rho_{ij}Y_{ja}\) |
| A1 | ph-ph id | \(X_{ai}Y_{ai}\) | \(\hat j_h^2\) | \(-\,(-1)^{j_a+j_i+\lambda}\hat\lambda^{-1} X_{ai}Y_{ia}\) |
| B5 | ph-ph | \(-X_{ai}\rho_{ba}Y_{bi}\) same h | \(-\rho\,\hat j\) | \(+\,(-1)^{j_a+j_i+\lambda}\hat\lambda^{-1} X_{ai}\rho_{ba}Y_{ib}\) |
| C1 | ppvv | \(X_{abij}\rho_{ijkl}Y_{abkl}\) | \(\hat J_{vv}^{-1}\hat\lambda^{-1}\) (overlap already has \(\rho\cdot\hat J_{vv}\)) | \((-1)^{J_0+J_1+\lambda}\hat J_{vv}^{-1}\hat\lambda^{-1} X^{J_0 J_1\lambda}\rho^{J_1}Y^{J_1 J_0\lambda}\) |
| B1 | pphv | \(X_{abei}\rho_{ij}Y_{abej}\) | \(\rho\,\hat J^2/\hat j_c\) | \(\delta_{j_j j_i}(-1)^{J_0+J_1+\lambda}\hat\lambda^{-1} X_{eiab}^{J_0 J_1} \rho Y_{abej}^{J_1 J_0}\) |
| C4 | pphv | \(X_{abcd}\rho_{dfae}Y_{bfce}\) + 3 perms | `Core_Diagram` 9j | \(-\,(-1)^{J_1+\lambda+J_3+J_4+J_5+j_b+j_e}\hat J_0\hat J_1\hat\lambda^{-1}\hat J_3^2\hat J_4\hat J_5\hat\jmath_0^2\) \(6j\{J_1\lambda J_0;j_c j_d j_0\}\) \(6j\{J_4\lambda J_5;j_c j_e j_0\}\) \(9j\{j_b j_a J_1;j_f J_3 j_d;J_4 j_e j_0\}\), \(X_{cdab}\), \(Y_{bfce}\) |
| A2 | pphh id | \(X_{abij}Y_{abij}\) | \(\hat J^2\) | \((-1)^{J_0+J_1+\lambda}\hat\lambda^{-1} X Y_{ijab}\) |
| C3 | pphh | \(X_{abij}\rho_{abcd}Y_{cdij}\) | \(\hat J_{pp}^{-1}\hat\lambda^{-1}\) (overlap already has \(\rho\cdot\hat J_{pp}\); same leftover packaging as C1) | \((-1)^{J_0+J_1+\lambda}\hat\lambda^{-1} X\rho^{J_0}Y_{ijcd}\) |
| B2 | pphh | \(-X_{abij}\rho_{ca}Y_{cbij}\) + 3 perms | \(\hat J^2/\hat j\) | \(-\delta_{j_c j_a}(-1)^{J_0+J_1+\lambda}\hat\lambda^{-1} X\rho Y_{ijcb}\) |
| B3 | pphv-ph | \(X_{cdab}\rho_{bd}Y_{ac}\) (X daggered) | \(\hat J^2/\hat j\) | \(\delta_{j_d j_b}(-1)^{J_{ab}+j_a+j_b}\hat J_{ab}\hat J_{hv}\hat\lambda^{-1}6j\{\lambda J_{ab} J_{hv};j_b j_c j_a\}\), times \((-1)^{J_{ab}-J_{hv}}\) |
| C5 | pphv-ph | \(-X_{cdab}\rho_{abed}Y_{ec}\) (X daggered) | \(\hat J_{hv}\hat\lambda^{-1}\) (overlap already has \(\rho\cdot\hat J_{ab}\)) | \(-(-1)^{J_{ab}+j_d+j_e}\hat J_{hv}\hat\lambda^{-1}6j\{J_{hv} J_{ab}\lambda;j_e j_c j_d\}\), times \((-1)^{J_{ab}-J_{hv}}\) |
| C2 | ppvv-qv | \(X_{abij}\rho_{akij}Y_{kb}\) | \(\hat J_{qv}\hat\lambda^{-1}\) (not \(\hat J_{qv}\hat J_{vv}\hat\lambda^{-2}\); overlap already has \(\rho\cdot\hat J_{vv}\)) | \((-1)^{J_0+j_a+j_b}\hat J_0\hat\lambda^{-1}6j\{J_0 J_1\lambda;j_k j_b j_a\}\), ρ at \(J_1\) |

B5 tensor **input** has a leading minus (same as scalar); AMC recoupling cancels it.
B2’s four terms are m-scheme particle exchanges, not Pandya.

C4 uses the same four particle-exchange perms as scalar `ConstructNormMatrix`.
`Core_Diagram_tensor` is the AMC of daggered \(X_{cdab}\) and stored \(Y_{bfce}\).
The caller multiplies \((-1)^{J_{ab}-J_{hv}}\) and the same `Ket::Phase` factors as the scalar loop.

## Scalar AMC vs code (hats)

Index structure matches. Prefactors that differ (do not treat as angular-momentum bugs):

- B4: AMC \(\hat j_a^2\), code \(\hat j_v\)
- C1/C3/C5/C2: AMC \(\hat J^2\), code \(\hat J\)
- B1/B2/B3: AMC \(\hat J^2\), code \(\hat J^2/\hat j\)

ρ in the kernel is the **scalar** valence RDM (`RdmOB` / `RdmTB_J`) even when χ is a tensor.
