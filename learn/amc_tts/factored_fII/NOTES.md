# Factorized TTS \(f^{\mathrm{II}}\) / \(\chi^\beta\)

AMC inputs/outputs under `input/` and `output/`.

## Declarations (critical)

| Object | `scalar` | Why |
|---|---|---|
| \(\Omega\) | `false` | tensor rank \(\lambda\) |
| \(\Gamma\) | `true` | scalar |
| \(\chi^\beta\) | **`false`** | tensor intermediate \(\Omega\Gamma\) |
| \(f^{\mathrm{II}}\) | **`true`** | final scalar (\(j_i=j_j\)) |

## Equations (AMC)

**\(\chi^\beta_{de}^\lambda\)** (`chi_beta.tex`) — \(\tfrac12\) kept here:

\[
\chi^\beta_{de}^{\lambda}
=
\tfrac12 (-1)^{j_e+\lambda}
\sum_{abc J_0 J_1}
(\,n_a n_b\bar n_c\bar n_e-\bar n_a\bar n_b n_c n_e\,)
(-1)^{J_0+j_c}\,\hat J_0\hat J_1
\begin{Bmatrix}\lambda & J_1 & J_0\\ j_c & j_d & j_e\end{Bmatrix}
\Gamma_{cdab}^{J_0}\,\Omega_{abce}^{J_0 J_1\lambda}.
\]

Occupation follows free index \(e\) (not \(d\)). No \(\hat j_d^{-2}\): store as reduced tensor 1b ME.

**\(f^{\mathrm{II}}_a\)** (`f2a_from_chi.tex`) — tensor\(\times\)tensor\(\to\)scalar:

\[
f_a_{ij}
=
\delta_{j_i j_j}(-1)^{j_i}\hat j_i^{-2}
\sum_{de J_0 J_1}
(-1)^{J_1+j_d}\,\hat J_0\hat J_1\hat\lambda^{-1}
\begin{Bmatrix}J_1 & \lambda & J_0\\ j_e & j_i & j_d\end{Bmatrix}
\chi^\beta_{de}^{\lambda}\,\Omega_{eidj}^{J_0 J_1\lambda}.
\]

(\(\lambda\) fixed to \(\mathrm{rank}(\Omega)\); AMC also writes a dummy sum over \(\lambda_0\).)

**\(f^{\mathrm{II}}_b\)**: same with \(\Omega_{ejdi}\). Assemble

\[
f^{\mathrm{II}}=f_a+h_\Gamma\,f_b
\]

(see parent `NOTES.md` hermiticity note: AMC cannot orient \(\Omega_{diej}\)).

## Implementation

- ethS: `FactorizedDoubleCommutator_eths::comm223_231_chi1b_tensor` (`use_TypeII_1b`)
- Reference extract: `ReferenceImplementations::comm223_231_tts_fII`
- Benchmarks: `run/test_tts_fII.py` (DIRECT↔PathB), `run/test_tts_f_mscheme.py` (m↔J)

**Phases (no tune):** AMC prints \((-1)^{j_i}(-1)^{J_1+j_d}\). ethS uses the
**combined** integer form \((-1)^{(j_{2i}+j_{2d})/2+J_1}\) only — do not multiply
an extra `phase((jp2+1)/2)` (that double-counts). DIRECT `tts_fII` applies
`phase((jp2+1)/2)` once because its unfactored loop phase has only
\(j_c{+}j_d{+}j_e\). See `learn/factorized_code_analyze.tex` §verify-chain /
§code-fII. Verified: m ≡ DIRECT ≡ Path B.
