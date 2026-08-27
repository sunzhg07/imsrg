# AMC check: IMSRG `comm222_phst` Pandya equations

**Verdict: forward Pandya (scalar + tensor) follows AMC exactly.** Inverse uses the same 6j/9j structure; scalar inverse includes AMC’s overall minus, tensor inverse **intentionally omits** it (see code comment).

## What was checked

AMC I/O: `learn/amc_tts/comm222_phst/{input,output}/`

| Piece | AMC | IMSRG |
|---|---|---|
| Scalar forward | `scalar_pandya*.tex` | `DoPandyaTransformation` / ref `Xbar` |
| Tensor forward | `tensor_pandya*_ninej.tex` | `DoTensorPandyaTransformation_SingleChannel` / ref `Ybar` |
| Scalar inverse | `scalar_pandya_inv.tex` | `AddInversePandyaTransformation` (`commij -= …`) |
| Tensor inverse | `tensor_pandya_inv_ninej.tex` | `AddInverseTensorPandyaTransformation` (`commij += …`) |
| Product | `ph_product_pandya*.tex` | coupled DGEMM (see below) |

## Index map

AMC Pandya `scheme=((1,-4),(3,-2))` writes `_ijkl` but couples **(iℓ)(kj)**.

IMSRG stores \(\langle ab J_{\mathrm{bra}}|\,\bar O\,|cd J_{\mathrm{ket}}\rangle\) from normal ME \(O(a,d,c,b)\).

**Map:** evaluate AMC at \((i,j,k,l)=(a,d,c,b)\) ↔ IMSRG channel \(\langle ab|cd\rangle\).

## Numeric kernel match

### Forward (exact)

- **Scalar:** AMC \(-\hat{J}'^{2}\{j_d\,j_c\,J;\,j_b\,j_a\,J'\}\) at `adcb` = IMSRG \(-\hat{J}'^{2}\{j_a\,j_b\,J;\,j_c\,j_d\,J'\}\).
- **Tensor (ninej):** full coeff (overall minus + phase + hats + 9j) matches IMSRG for all sampled \((j,J,\lambda)\), including \(\lambda=0\).

### Inverse

- **Scalar:** AMC overall minus present in IMSRG via `commij -= (2J'+1) sixj me`.
- **Tensor:** same 9j/phase/hats as AMC `Y=-barY`, but IMSRG uses `+=` instead of AMC’s overall minus → **IMSRG = −AMC** on the isolated inverse kernel.

This is documented in `Commutator.cc`:

> “technically there's a minus sign missing in what is done here, but that's on purpose.”

So the *isolated* tensor inverse does not match naive AMC `Y_ijkl = - barY_ijkl`; the *ph commutator pipeline* (two forward Pandya minuses × product × this inverse) is what was validated against m-scheme (`ReferenceImplementations::comm222_phst`).

## Product / DGEMM

AMC expands \(\sum_{ab}\bar X_{ilab}\bar Y_{abkj}\) with extra 6js (m-scheme SP contraction). IMSRG multiplies already J-coupled Pandya blocks at shared intermediate \(J\) (DGEMM). Same physics, different packaging.

## Takeaway for ethS \(f^{\mathrm{III}_a}\)

IMSRG `comm222_phst` / tensor Pandya **does** follow AMC on the forward transform that builds \(\bar\Omega\). χ mismatches vs AMC \(\chi^\gamma\) are from Ω×Ω / occupancy / leg packaging, not from “IMSRG Pandya ≠ AMC Pandya.”
