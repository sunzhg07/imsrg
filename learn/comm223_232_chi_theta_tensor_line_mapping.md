# comm223_232 chi_theta tensor line mapping

Scope: chi^theta path only in comm223_232_chi2b.

## Equation reference

- chi^theta in diag2_compact:
  amc/examples/sample_output/diag2_compact.tex:84

Equation structure used:
- two Eta operators (Eta*Eta)
- occupation factor: (n_a n_b nbar_k + nbar_a nbar_b n_k + n_a n_b nbar_j + nbar_a nbar_b n_j)
- final chi^theta is scalar.

## Symbol -> code mapping

1. Branch gate
- Role: tensor chi^theta only when eta is tensor and Gamma scalar.
- Symbol: tensor_eta_case.
- Anchor: src/FactorizedDoubleCommutator_eths.cc:2260.

2. Scalar legacy chi^theta (unchanged)
- Symbol: CHI_IV[ch] = Eta_matrix * Eta_matrix_c + (Eta_matrix * Eta_matrix_d)^T.
- Anchors:
  - src/FactorizedDoubleCommutator_eths.cc:2261
  - src/FactorizedDoubleCommutator_eths.cc:2262

3. Tensor chi^theta construction (new)
- Role: explicit Eta*Eta sum, scalar output matrix CHI_IV[ch].
- Tensor block starts: src/FactorizedDoubleCommutator_eths.cc:2264.

4. External indices for CHI_IV element
- ibra_ext -> (i,j)
- iket_ext -> (k,l)
- Anchors:
  - src/FactorizedDoubleCommutator_eths.cc:2267-2276
  - src/FactorizedDoubleCommutator_eths.cc:2282-2291

5. Internal summation indices
- a,b over all_orbits
- J2 coupling loop with triangle(J0,J2,lambda)
- Anchors:
  - src/FactorizedDoubleCommutator_eths.cc:2298-2319

6. Occupation factor for chi^theta
- Code:
  occ = n_a*n_b*nbar_k + nbar_a*nbar_b*n_k + n_a*n_b*nbar_j + nbar_a*nbar_b*n_j
- Anchor: src/FactorizedDoubleCommutator_eths.cc:2307-2308.

7. Two-Eta product term
- eta1 = Eta.TwoBody.GetTBME_J(J0, J2, i, j, a, b)
- eta2 = Eta.TwoBody.GetTBME_J(J2, J0, a, b, k, l)
- Contribution:
  chi_theta += occ * phase(J2 + lambda) * hat_lambda_inv * eta1 * eta2
- Anchors:
  - src/FactorizedDoubleCommutator_eths.cc:2314-2320
  - src/FactorizedDoubleCommutator_eths.cc:2324-2325

8. Scalar output normalization
- CHI_IV[ch](ibra_ext, iket_ext) = phase(J0) * chi_theta / (2J0+1)
- Anchor: src/FactorizedDoubleCommutator_eths.cc:2335-2336.

## Reuse of scalar downstream contraction

After CHI_IV is built (scalar object), downstream code remains unchanged and reused:
- bar_CHI_IV recoupling and later chi*Gamma->Z path still use existing scalar machinery.
- Initial contraction block including Chi_III_Op path anchor:
  src/FactorizedDoubleCommutator_eths.cc:1803.

This matches your requirement: chi^theta computed from Eta*Eta in tensor-aware way, then treated as scalar for the rest.
