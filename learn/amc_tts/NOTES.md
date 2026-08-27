# TTS unfactorized J-coupled equations (trusted source only)

**Trusted M-scheme source:** `amc/examples/sample_output/factorized_code_analyze.tex` §unfact only.

**Do not trust** for TTS physics: `amc/examples/sample_input/reference_tensor.txt`,
`amc/examples/sample_output/reference_tensor.tex`, or the existing
`comm223_*_tts_BruteForce` implementations derived from those files.

**Production code:** `comm223_231_tts` / `comm223_232_tts` in
`src/ReferenceImplementations.cc` (regenerator: `run/gen_tts_trusted.py`).

## Setup

- $\Omega$ / `Eta`: tensor (`scalar=false`), rank $\lambda$
- $\Gamma$ / `Gamma`: scalar (`scalar=true`)
- Final $Z$: scalar one- and two-body
- IMSRG defaults: $h_\Omega=-1$, $h_\Gamma=+1$ (kept general in formulas)
- Two-body AMC inputs omit $\hat P$; restore $(1-\hat P_{ij})$ / $(1-\hat P_{kl})$ as in analyze

## Tensor reduced-ME storage (IMSRG++)

Tensors (`rank_J>0`) are stored as **reduced** MEs by default (`is_reduced=true`).

**One-body** (`Operator::SetOneBody` / `Symmetrize`):

$$
O_{ji} = h\,(-1)^{j_i-j_j}\,O_{ij}
$$

(`phase((j2_i-j2_j)/2)` in code; $h=+1$ Hermitian, $h=-1$ anti-Hermitian).

**Two-body** (`TwoBodyME::GetTBME` when `ch_bra>ch_ket`):

$$
\langle J_b\|O\|J_k\rangle
=
h\,(-1)^{J_b-J_k}\,
\langle J_k\|O\|J_b\rangle
$$

Access tensor MEs with different bra/ket $J$ via `GetTBME_J(J_bra,J_ket,a,b,c,d)`.
Scalar $Z$ uses equal-$J$ channels only.

## AMC inputs / outputs

| Diagram | Input | Notes |
|---|---|---|
| $f^{\mathrm{I}}$ | `input/f1.txt` | two $\Gamma$ strings |
| $f^{\mathrm{II}}$ | `input/f2a.txt` + `f2b.txt` | see hermiticity note |
| $f^{\mathrm{III}_a}$ | `input/f3a.txt` | |
| $f^{\mathrm{III}_b}$ | `input/f3b.txt` | |
| $\Gamma^{\mathrm{I..IV}_c}$ | `input/G*.txt` | direct strings; restore $P$ |

Outputs: `output/*.tex`. Raw dump: `EQUATIONS_RAW.md`.

### $f^{\mathrm{II}}$ hermiticity note

AMC cannot orient $\Omega_{diej}$ (Yutsis). Use bra–ket swap
$\Omega_{diej}=h_\Omega\,\Omega_{ejdi}$ (m-scheme, real) so

$$
h_\Omega h_\Gamma\,\Omega_{diej}=h_\Gamma\,\Omega_{ejdi}.
$$

Assemble $f^{\mathrm{II}}=f_{2a}+h_\Gamma\,f_{2b}$ with AMC of $\Omega_{eidj}$ and $\Omega_{ejdi}$.

### Overall-sign caveat (tensor AMC)

For some tensor reductions AMC drops or flips the overall rational prefactor
relative to the m-scheme input (e.g.\ scalar $G2$ keeps $-\tfrac12$, tensor $G2$
prints $+\tfrac12$; $G3a$/$G3b$/$G4a$ lose the leading $-1$).
**Implementation rule:** take 6j / hat / phase structure from AMC, but set the
overall rational prefactor from analyze §unfact (and restore $(1-\hat P)$ by hand).

## Name map (analyze ↔ old reference_tensor labels)

| Analyze | Old TTS BruteForce label |
|---|---|
| $f^{\mathrm{I}}$ | I |
| $f^{\mathrm{II}}$ | IIIa, IIIb |
| $f^{\mathrm{III}_a}$ | IIa, IIc |
| $f^{\mathrm{III}_b}$ | IIb, IId |
| $\Gamma^{\mathrm{I}}$ | Ia, Ib |
| $\Gamma^{\mathrm{II}}$ | IVa, IVb |
| $\Gamma^{\mathrm{III}_a}$ | IIa, IIc |
| $\Gamma^{\mathrm{III}_b}$ | IIb, IId |
| $\Gamma^{\mathrm{III}_c}$ | IIe, IIf |
| $\Gamma^{\mathrm{IV}_a}$ | IIIa, IIIb |
| $\Gamma^{\mathrm{IV}_b}$ | IIIc, IIId |
| $\Gamma^{\mathrm{IV}_c}$ | IIIe, IIIf |

## Factorized ethS (implemented)

**Living doc (J-scheme equations + optimizations):**  
[FACTORIZED_TTS_IMPLEMENTED.md](FACTORIZED_TTS_IMPLEMENTED.md)

Update that file when adding diagrams or changing optimized prefactors/storage.

## Lessons learned (factorized vs unfactorized)

See also [LESSONS.md](LESSONS.md) for a short agent-oriented checklist.

Factorized \(f^{\mathrm{II}}\) / tensor \(\chi^\beta\): [factored_fII/NOTES.md](factored_fII/NOTES.md).

