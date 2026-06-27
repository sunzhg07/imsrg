# Validation Status: 223_231 Tensor vs diag2_compact

## Verdict
- Implementation completeness (tensor branch pieces present): YES
- Strict equation-level validation against diag2_compact: PARTIAL
- Numerical validation against independent reference: NOT YET

## Scope Validated
- File: `src/FactorizedDoubleCommutator_eths.cc`
- Branch: `Eta` tensor, `Gamma` scalar, `Z` scalar in `comm223_231`

## Term-by-Term Status

| Block | In Code | Structural Match to diag2 | Prefactor/Phase Fully Verified | Numeric Regression |
|---|---|---|---|---|
| `chi^alpha` | Yes | Mostly yes | Partial | No |
| `f^(I)` | Yes | Mostly yes | Partial | No |
| `chi^beta` | Yes | Mostly yes | Partial | No |
| `f^(II)` | Yes | Improved after J-range fix | Partial | No |
| `bar{Omega}` recoupling helper | Yes | Improved after triangle-domain tightening | Partial | No |
| `bar{chi}^gamma` (`III_a` intermediate) | Yes | Partial | Partial | No |
| `chi^delta` (`III_b` intermediate) | Yes | Mostly yes | Mostly yes (normalization item resolved) | No |
| `f^(III_a)` | Yes | Partial | Partial | No |
| `f^(III_b)` | Yes | Mostly yes | Partial | No |

## What is already resolved
1. f(II) used independent `J0` and `J1` loop ranges (no longer tied).
2. `chi_delta` normalization cancellation (`hat(J0)^2` outer with `hat(J0)^{-2}` inside) is consistent in current coding convention.
3. Barred recoupling helper had stricter triangle-domain enforcement added.

## Remaining blockers for calling this “validated”
1. End-to-end term-by-term numerical comparison against trusted values from AMC/diag2 workflow.
2. Final phase/sign convention audit for all barred quantities in `III_a` path.
3. Small-space deterministic regression test in `UnitTest` for tensor `223_231`.

## Testing availability check
- Existing `UnitTest.cc` references `comm223_231` in scalar/factorized contexts.
- No dedicated tensor `FactorizedDoubleCommutator_eths` regression found yet.

## Bottom line
- We are **implementation-complete** for tensor `223_231` blocks in this branch,
  but **not fully validated** against `diag2_compact` until numerical and final sign checks are added.
