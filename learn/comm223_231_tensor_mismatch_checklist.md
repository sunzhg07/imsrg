# 223_231 Tensor Mismatch Checklist (Code vs diag2_compact)

## Scope
- Target implementation: `src/FactorizedDoubleCommutator_eths.cc`
- Tensor branch: `Eta` tensor, `Gamma` scalar, output scalar.
- Reviewed pieces: chi1b tensor (`I`, `II`) and chi2b tensor (`III_a`, `III_b`).

## Resolved Items
- Barred Omega index-order fix (cross-coupled coding convention):
  - `barred_tbme(Op,a,j,k,b,...)` now maps to underlying TBME access
    `Op.TwoBody.GetTBME_J(J2,J3,a,b,j,k)`.
  - This follows scalar Pandya coding convention where barred indices are
    interpreted in cross-coupled operator order rather than direct leg order.
- f(II) closure loop ranges decoupled:
  - Before: `J1` reused `J0`-constrained range.
  - Now: `J0` from `(b,p)` coupling, `J1` from `(a,q)` coupling, plus triangle filter.
  - Files updated:
    - `src/FactorizedDoubleCommutator_eths.cc`
    - `learn/comm223_231_chi1b_tensor_code_equations.tex`
- `chi_delta` normalization/prefactor consistency:
  - Implemented `chi_delta` carries `0.25 * (-1)^{J0} / (2J0+1)`.
  - The final `f^(III_b)` accumulation multiplies by `(2J0+1)`.
  - Net effect reproduces the expected `hat(J0)^{-2}` (inside `chi_delta`) with
    outer `hat(J0)^2` cancellation pattern from `diag2_compact` structure.
  - No code change required for this specific item.
- `barred_tbme` summation-domain tightening:
  - Added explicit triangle guards for external couplings:
    - `(j_a,j_b,J0)`, `(j_j,j_k,J1)`, `(J1,J0,lambda)`.
  - Tightened `j0` loop bounds to the overlap implied by the 6j chain, instead
    of a looser orbit-only range.
  - This aligns the coded summation domain more closely with the tensor Pandya equation.

## Open Review Items
- Confirm overall prefactors for `chi^alpha`, `chi^beta`, `f^(I)`, `f^(II)` against exact `diag2_compact.tex` hats/phases.
- Confirm `barred_tbme` recoupling chain matches the chosen Pandya convention term-by-term (ordering and phase conventions).
- Confirm `bar_chi_gamma` occupancy and phase prefactors against `\bar\chi^\gamma` in `diag2_compact.tex`.
- Add numerical validation against a trusted reference for tensor 223_231 (small-space regression test).

## Suggested Next Validation Pass
- Derive one explicit low-dimensional test case and compare:
  - code output from tensor 223_231 branch,
  - direct reference implementation (or AMC-reduced expression evaluation),
  - term-by-term decomposition (`I`, `II`, `III_a`, `III_b`).
