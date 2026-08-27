# \(f^{\mathrm{III}_a}\) DIRECT / AMC packaging audit

Status after step 1 (`run/test_tts_fIIIa_mscheme.py`): **m ≠ DIRECT**.
This note records what we checked in the AMC → DIRECT chain (no Path B yet).

## Locked facts

| Check | Result |
|---|---|
| AMC input | `learn/amc_tts/input/reduced/f3a.txt` (= analyze §unfact string) |
| AMC output (re-run) | identical to `output/case2_reduced/f3a.tex` |
| `tts_fIIIa` vs literal AMC tex | **bit-match** (\(\max\|\Delta\|\sim10^{-15}\)) |
| `tts_fIIIa` vs ethS slow (\(W_1W_2\)) | **bit-match** |
| Phase integerization variants | none match m (orbit-dependent ratios) |
| Denom \(\hat\jmath^{-1}\) vs \(\hat\jmath^{-2}\) vs none | none match m |
| WE convention `wigner` vs `sakurai` | formulas **differ**; not yet coded vs m |
| Known AMC caveat ([NOTES.md](../NOTES.md)) | tensor AMC can drop/flip overall rational prefactor |

**Conclusion so far:** the C++ DIRECT is a faithful transcription of AMC `case2_reduced/f3a`.
The mismatch with m is **not** a coding typo in sixj / hats / TBME index order.
Either (a) AMC’s reduced equation is wrong for this network, (b) our m-oracle CG packaging is wrong for this topology, or (c) `reduce=` / WE convention packaging for the *inputs* is wrong.

## AMC input packaging (case 2)

```text
declare Omega { mode=4, scalar=false }           # tensor → always reduced RME
declare Gamma { mode=4, scalar=true, reduce=true }
declare f3a   { mode=2, latex="f^{(III_a)}", scalar=true, reduce=true }

f3a_ij = sum_abcde((nbar_a*nbar_b*n_c*n_d - n_a*n_b*nbar_c*nbar_d)
    * (Omega_abcd*Omega_idae*Gamma_cejb - Omega_abcd*Omega_edaj*Gamma_cieb));
```

From AMC `docs/ug.rst`:
- nonscalar tensors **always** reduced;
- `reduce=true` on scalars ⇒ reduced MEs (default for scalars is *unreduced*);
- default `--wet-convention wigner` ≡ Edmonds (IMSRG).

TTS packaging matches that recipe:
- `Eta` reduced tensor, `Gamma.MakeReduced()`, body has AMC’s \(\hat\jmath_i^{-1}\), store unreduced via \(\div(2j+1)\equiv\div\hat\jmath\) after that.

Same recipe works for \(f^{\mathrm{I/II/III}_b}\) at \(\lambda\neq0\). So case-2 flags alone are not an obvious smoking gun.

## AMC output (Term 1 / Term 2) — what DIRECT implements

**Term 1** (overall −):
\[
-\delta_{j_i j_j}\,\hat\jmath_i^{-1}
\sum (-1)^{J_1+J_3+j_b+j_e+\lambda}\,
\hat J_0\cdots\hat J_4\,\hat J_5^{2}\,\hat\jmath_0^{2}\,\hat\lambda^{-1}
\times(5~\mathrm{sixj})\,
\Omega^{J_0 J_1\lambda}_{abcd}\,
\Omega^{J_2 J_3\lambda}_{idae}\,
\Gamma^{J_4}_{cejb}
\]

**Term 2** (overall \(+(-1)^{j_i}\)):
\[
+\delta_{j_i j_j}\,(-1)^{j_i}\,\hat\jmath_i^{-1}
\sum (-1)^{J_1+J_3+j_b+\lambda}\,
(\text{same hats})\,
\times(5~\mathrm{sixj})\,
\Omega^{J_0 J_1\lambda}_{abcd}\,
\Omega^{J_2 J_3\lambda}_{edaj}\,
\Gamma^{J_4}_{cieb}
\]

`tts_fIIIa` implements this literally (combined integer phases for half-integer \(j\)).

## What is *not* wrong

- Sixj argument order vs AMC tex  
- \(\Omega\) / \(\Gamma\) GetTBME_J index order vs AMC superscripts  
- \(\hat J_4\) once with `MakeReduced(Γ)` (≡ unreduced \(\hat J_4^{2}\Gamma\))  
- Double-counting an extra \((-1)^{j_i}\) on top of the combined phase (checked)  
- Path B (not in scope until DIRECT ≡ m)

## Suspects (ordered)

1. **Isolate \(\chi^\gamma\) first** (`chi_gamma_direct.txt` → `chi_gamma_direct_plain.tex`)  
   **DONE — FAIL** (`run/test_chi_gamma_mscheme.py`, emax=1 λ=2):  
   m-scheme \(\chi^\gamma\) with \([\Omega\times\Omega]^{(0)}\) ≠ AMC/ethS W1×W2 DIRECT.  
   Ratios orbit/m-dependent (not a global phase/λ̂).  
   Consistent with [OMEGA_TT_TO_SCALAR.md](OMEGA_TT_TO_SCALAR.md): **AMC printed TT→0 direct is already wrong**; use Neithan Path A / corrected Path B.

2. **Overall rational / sign from analyze** ([NOTES.md](../NOTES.md) tensor caveat)  
   Secondary once a Neithan-correct χ exists.

3. **`--wet-convention sakurai`**  
   Regenerated f3a differs in hat powers; IMSRG uses wigner/Edmonds — not the primary issue given (1).

4. **Γ ladder / f3a assembly**  
   Blocked until χ^γ DIRECT ≡ m (or Neithan χ ≡ m).

5. **λ=0 red herring**  
   Rank-0 `RandomOp` Ω is unreduced while AMC/TTS assume reduced rank-λ. Fair tests use \(\lambda\neq0\).

## Next action

Benchmark **m vs Neithan Path A** (and/or Path B with \(-\hat\lambda(-1)^{J+\lambda}\)) for \(\chi^\gamma\).  
Do **not** fix Path B or f3a ladder against AMC direct / `tts_fIIIa`.
