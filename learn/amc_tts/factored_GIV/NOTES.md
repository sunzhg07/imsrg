# Γ^{IV_a,b,c} / χ^{ι,κ,λ} (Ω×Γ)

Cross-cutting playbook: [../LESSONS.md](../LESSONS.md) §Playbook.

## Product-rank rule

| × | scalar | tensor |
|---|---|---|
| **scalar** | scalar | tensor |
| **tensor** | tensor | **scalar** |

- III (η, θ): \([\Omega\times\Omega]^0\) → **scalar** χ
- **IV (ι, κ, λ):** \([\Omega\times\Gamma]^\lambda\) or \([\Gamma\times\Omega]^\lambda\) → **tensor** χ (rank λ)
- Ladder: χ(T)×Ω(T) → **scalar** Z (allowed). Never invent \(T\times T\to T\) for χ itself.

## Naming map (analyze ↔ Factorized)

| Analyze | Factorized object | Topology |
|---|---|---|
| Γ^{IV_a} / χ^κ ↔ `Chi_VI_II`; χ^{κ'} ↔ `Chi_VI` | `bar_CHI_VI_II` / `bar_CHI_VI` | Inv Pandya → −Ωχ_VI − χ_VI,II Ω |
| Γ^{IV_b} / χ^ι | `bar_CHI_V` → RC → `CHI_V_final` | Pandya → RC → ×Ω̄ → Inv Pandya |
| Γ^{IV_c} / χ^λ | `CHI_VII` → `bar_CHI_VII_CC` | Normal → CC → ×Ω → Inv CC (+½) |

## Status (ethS)

| Piece | λ=0 | λ≠0 | Gold chain (m≡DIRECT≡Path B) |
|---|---|---|---|
| IV_a / χ^κ | **one** Path B any λ | rectangular (= square at λ=0) | **DONE** m≡AMC≡Path B≡ethS (λ≠0); λ=0 Factorized fork = packaging only |
| IV_b / χ^ι | Pandya→RC→×Ω̄→inv (λ=0); χ→fold (any λ) | rectangular Path B χ + fold | **DONE** full \(Z\): m ≡ AMC_direct ≡ Path B |
| IV_c / χ^λ | CHI_VII→CC→DGEMM→inv (+½) | Path B Pandya | **DONE** |

## Γ^{IV_a} / χ^κ — gold chain

**Trusted equation** ([factorized_code_analyze.tex](../../factorized_code_analyze.tex) `\eqref{eq:chikappa}`):

\[
\chi^\kappa_{ijbd}
=\sum_{ac}
\bigl(\bar n_c\bar n_d n_a+n_c n_d\bar n_a\bigr)
\,\Omega_{aicd}\,\Gamma_{jcba}.
\]

The **old arxiv** print \(\Omega_{ajbl}\Gamma_{ibka}\) with occ on free \(k\) is **obsolete**
(explicitly marked “Fixes vs.\ older print” / summary table in analyze). Do not use
`chi_kappa.txt` / `arxiv_eq.tex` as gold.

AMC twin: `input/chi_kappa_analyze.txt` → `output/chi_kappa_analyze.tex`.

### Path B (Factorized) — χ^κ ↔ `Chi_VI_II`

From analyze §code-pandya-chi / §code-GIVa:

\[
\bar\chi^{\kappa,\mathrm{II}}
=h_\Omega\,(\mathtt{nnnbar\_Eta\_d})^{T}\cdot\bar\Gamma
\quad(\mathtt{bar\_CHI\_VI\_II}),
\qquad
\bar\chi^{\kappa'}
=\bar\Gamma\cdot\mathtt{nnnbar\_Eta\_d}
\quad(\mathtt{bar\_CHI\_VI}).
\]

Occ = `occ_ABbarD` family. **Bare χ^κ** = inv Pandya of `bar_CHI_VI_II` only
(direct term; never AS χ; never Pandya-exchange bake-in / √2). `bar_CHI_VI` is χ^{κ'}.

**Bare χ^κ Path B ≡ AMC direct** — `run/test_chi_kappa_pathB_vs_direct.py` — **PASS all λ**
(emax=1, tensor/WE Ω):

| Step | Formula |
|---|---|
| Pandya Γ | IMSRG scalar adcb, 1×6j |
| Pandya Ω | IMSRG tensor adcb, **rectangular** \(J_{\mathrm{bra}}\neq J_{\mathrm{ket}}\), 1×9j |
| DGEMM | \(\bar\chi^{J_0 J_1}(il;kj)=h_\Omega(-1)^{J_0+J_1}\sum_{ab}\mathrm{occ}_{ABbarD}(a,b,l)\,\bar\Omega^{J_1 J_0}(ab;il)\,\bar\Gamma^{J_1}(ab;kj)\) (`VI_II`) |
| Inv | IMSRG `AddInverseTensorPandya` kernel **or** AMC via Eq4 (no printed leading minus) |

\((-1)^{J_0+J_1}\) is required for λ≠0 (trivial +1 when \(J_0=J_1\)). Do **not** use AMC `tensor_pandya_inv` leading − (gives overall −χ).

**Bare χ^κ m ↔ AMC direct:** `run/test_chi_kappa_m_vs_amc.py` — **PASS** for
**all λ including 0** (emax=1). Ω is WE-reduced on the **tensor path for every λ**;
λ=0 does **not** mean unreduced scalar.

**Ladder twin (full Z):** `run/test_G4a_pathB_mscheme.py` — Path B ≡ m (**all λ**, one method).
ethS: `run/test_tts_GIVa_eths_vs_pathB.py` — **PASS** λ≠0. λ=0 ethS still forks to Factorized
scalar CHI_VI (unreduced-Ω habit); physics Path B is the same — unify when convenient.

### Physics / debugging lessons (IV_a lock, 2026-08)

**1. One Path B for any λ — no second “scalar method”.**
Rectangular Pandya + VI_II + inv + ladder is continuous in λ. At λ=0: \(J_0=J_1\),
\((-1)^{J_0+J_1}=+1\), 9j→6j. Do **not** invent a separate algorithm for λ=0; only
name Ω packaging (always **WE-reduced** on this tensor path — λ=0 ≠ unreduced scalar).

**2. Reduced tensor transpose → phase (the big one).**
Factorized VI_II is \(h_\Omega\,(\mathrm{occ}\odot\bar\Omega)^{T}\cdot\bar\Gamma\).
For a **rank-λ reduced** ME, Hermitian conjugation / channel flip is **not** a bare `.t()`:

\[
\bar\Omega^{J_0 J_1}(il;ab)
=
h_\Omega\,(-1)^{J_0-J_1}\,
\bar\Omega^{J_1 J_0}(ab;il)
\quad\text{(real MEs)}.
\]

For integer \(J\): \((-1)^{J_0-J_1}=(-1)^{J_0+J_1}\). Same phase as `comm222_phst`:

```cpp
flipphaseY = hY * phase(Jbra - Jket);  // ≡ phase(Jbra+Jket) for integer J
```

So the locked DGEMM is

\[
\bar\chi^{J_0 J_1}
=
h_\Omega\,(-1)^{J_0+J_1}
\sum_{ab}\mathrm{occ}\,
\bar\Omega^{J_1 J_0}(ab;il)\,
\bar\Gamma^{J_1}(ab;kj).
\]

At equal \(J\) the phase is trivial and only \(h_\Omega\) remains (Factorized λ=0).
**Never** raw `arma::mat::t()` of a rectangular tensor block without this phase.

**3. Same diagram as `comm222_phst`, stop before AS.**
Pandya(Γ)+Pandya(Ω)+DGEMM+inv tensor = phst pipeline. Differences: occ `ABbarD` not
\((n_a-n_b)\); χ stays **non-Hermitian** — no \((1-P)\) on χ. Apply \((1-P)\) only on
the ladder intermediate \(W=-\chi\Omega\) (and \(W_{klij}\)).

**4. Inv sign: IMSRG / via Eq4, not printed AMC leading −.**
Isolated AMC `tensor_pandya_inv` has overall −; IMSRG `AddInverseTensorPandya`
intentionally omits it (`Commutator.cc` comment). Path B gold = **invPlus** (no
printed −). With printed − → overall −χ (r=−1). Round-trip: invPlus∘fwd ≈ −Id.

**5. \(h_\Omega\) once, in χ̄ build — not again on \(W_{klij}\).**
Hermitian assemble \(Z=(1-P_{ij})W+(1-P_{kl})W_{klij}\) is \(A+A^{T}\). Multiplying
\(h_\Omega\) on \(W_{klij}\) when both braces reuse the same χ^κ gives \(A-A^{T}\)
(anti-Hermitian). \(h_\Omega\) already sits in VI_II.

**6. Packaging checklist (before any m↔J).**
| Object | Packaging |
|---|---|
| Ω | **reduced** all λ (WE) |
| Γ | unreduced scalar |
| χ^κ | tensor WE (AMC analyze as printed) |
| W, Z | \(X_{\mathrm{red}}=\sum\mathrm{CG}\,\mathrm{CG}\,X(m)/\hat J\); ethS store \(Z_{\mathrm{unred}}=Z_{\mathrm{red}}/\hat J\) |
| `GetTBME_J` | undoes √2 store norms — compare without extra √2 on the Python side |

Fake FAIL ratios \(\sim\hat J^{\pm1}\), \(\sqrt{2}^{\pm1}\) → packaging, not 6j.

**7. Phases: combined-integer only.**
AMC may print \((-1)^{J+j}(-1)^{J'+j'}\). Evaluate as **one** integer phase via
\((j2_a+j2_b+\ldots)/2 +\) integer \(J\)s (χ^β lesson). Per-orbit `j2//2` drops a sign.

**8. Gold equation is analyze, not old arxiv.**
\(\chi^\kappa=\sum(\mathrm{occ})\,\Omega_{aicd}\Gamma_{jcba}\). Old \(\Omega_{ajbl}\Gamma_{ibka}\)
with occ on free \(k\) is obsolete.

Benches: `test_chi_kappa_{m_vs_amc,pathB_vs_direct}.py`, `test_G4a_{Wbra,pathB}_mscheme.py`,
`test_tts_GIVa_eths_vs_pathB.py`. See also [../comm222_phst/AMC_CHECK.md](../comm222_phst/AMC_CHECK.md).

Ladder after χ (both partners):

\[
Z\;{+}{=}\;
-\Omega\,\chi_{\mathrm{VI}}
-\chi_{\mathrm{VI,II}}\,\Omega
\quad\text{with }(1-P).
\]

Forward/inv Pandya minuses cancel (analyze §code-GIVa).

## Γ^{IV_b} / χ^ι — bare χ m ≡ AMC (2026-08)

**Trusted equation** (η-analog; Factorized `bar_CHI_V`):

\[
\chi^\iota_{ijkl}
=\sum_{ab}
\bigl(\bar n_a n_b\bar n_k+n_a\bar n_b n_k\bigr)
\,\Omega_{bjka}\,\Gamma_{iabl}.
\]

**Arxiv / older print \(\Omega_{bika}\Gamma_{iabl}\) is wrong** — free \(j\) is a
spectator (AMC Yutsis fails). It only coincides on fold patterns
\(\chi_{aibk}\), \(\chi_{akbi}\). Do not use it as gold.

AMC: `input/chi_iota_analyze.txt` → `output/chi_iota_analyze.tex`
(Ω·Γ product order → κ-like mid-\(j_0\) print).

**Bare χ^ι m ↔ AMC:** `run/test_chi_iota_m_vs_amc.py` — **PASS all λ** (emax=1;
Ω WE-reduced for every λ; Γ unreduced; combined-integer phases).

**Bare χ^ι Path B ≡ AMC** — `run/test_chi_iota_pathB_vs_direct.py` — **PASS all λ**
(no RC; Pandya→DGEMM→inv only):

| Step | Formula |
|---|---|
| Pandya Γ | IMSRG scalar adcb, 1×6j |
| Pandya Ω | IMSRG tensor adcb, **rectangular** \(J_{\mathrm{bra}}\neq J_{\mathrm{ket}}\), 1×9j |
| DGEMM | \(\bar\chi^{J_0 J_1}(il;kj)=\sum_{ab}\mathrm{occ}_{AbarBC}(a,b,k)\,\bar\Gamma^{J_0}(il;ab)\,\bar\Omega^{J_0 J_1}(ab;kj)\) (`CHI_V`) |
| Inv | IMSRG / AMC invPlus (**no** printed leading −; no AS) |

Unlike χ^κ VI_II: **no** \(h_\Omega\), **no** \((-1)^{J_0+J_1}\), **no** transpose — Factorized
`bar_CHI_V = bar_Gamma * nnnbar_Eta` as written. Pandya-space occ is
`occ_AbarBC` (≠ normal analyze weight; remapped by Pandya).

### Procedure (follow **scalar code**, not arxiv)

Same discipline as [factored_GIIIb/NOTES.md](../factored_GIIIb/NOTES.md):

1. **Read Factorized** — which Pandya→RC, exact index order (ignore arxiv / stale comments).
2. **AMC** — strip pack/rewire; same-label scheme change only; put pack back by hand.
3. **Tensor** — `scalar=false` on the same AMC string; keep code pack/rewire.

#### Γ^{IV_b} full \(Z\) (term1+term2) — **WRAP-UP LOCKED** (2026-08)

**Gold chain (any λ):**

\[
Z_m
\equiv
Z[\chi^{\mathrm{AMC\,direct}}]
\equiv
Z[\chi^{\mathrm{Path\,B}}]
\]

with the same analyze fold + pack:

\[
W=\sum_{ab}\bigl(\chi_{aibk}\Omega_{jbla}-\chi_{akbi}\Omega_{jalb}\bigr)
=W_1-W_2,
\qquad
Z=(1-P_{ij})(1-P_{kl})W.
\]

| Link | Bench | Status |
|---|---|---|
| \(\chi_m\equiv\chi^{\mathrm{AMC}}\) | `test_chi_iota_m_vs_amc.py` | **PASS** all λ |
| \(\chi^{\mathrm{Path\,B}}\equiv\chi^{\mathrm{AMC}}\) | `test_chi_iota_pathB_vs_direct.py` | **PASS** all λ |
| \(Z[\chi^{\mathrm{Path\,B}}]\equiv Z_m\) | `test_G4b_pathB_fold_mscheme.py` | **PASS** λ=0,1,2 |

So **m-scheme = AMC direct = Path B** for the full diagram. No TTS.

**Production Path B (any λ)** — same discipline as GIVa:

| Step | What | \(P\)? | Status |
|---|---|---|---|
| 1 | Pandya \(\bar\chi^\iota=\bar\Gamma\cdot(\mathrm{occ}\odot\bar\Omega)\) | no | **LOCKED** |
| 2 | InvPlus → normal \(\chi^\iota\) | no | **LOCKED** (≡ m ≡ AMC) |
| 3 | Analyze fold \(W=W_1-W_2\) | no | **LOCKED** |
| 4 | \(Z=(1-P_{ij})(1-P_{kl})W\) | **yes — only here** | **LOCKED** |

**Optional duals (not required for gold):**
| Dual | Scope | Bench |
|---|---|---|
| Factorized Pandya→RC→×Ω̄→Inv | λ=0 full \(Z\) | `test_G4b_factorized_fullZ.py` **PASS** |
| Path A term1 RC×RC | bare \(W_1\) | `test_G4b_pathA_term1_correct.py` **PASS** |
| Path A term2 = klij remap of term1 | bare \(W_2\), all λ | `test_G4b_pathA_term2_klij.py` **PASS** |
| Factorized term2 pack \(\mathrm{RC}[-h_Z\bar\chi^T]\) | λ=0 \(Z_2\) | `test_G4b_term2_pandya_RC.py` **PASS** |

Rectangular Factorized RC×Ω̄×Inv at λ≠0 still open as a speed twin — **do not block**; fold is gold.

**Do not use TTS.** Gold = m / AMC bare χ / Path B fold.

`run/test_G4b_pathB_fullZ.py` (AMC 03b→04) remains WIP alternate dual.

#### Path A term1 (RC×RC) — **LOCKED** bare \(W_1\), λ=2 (2026-08)

Correctness dual of Path B for term1 only:
\(W_1=\sum_{ab}\chi_{aibk}\Omega_{jbla}\) (no \((1-P)^2\), no pack).

Flow: **RC[Ω], RC[χ]** (scheme `((1,-3),(4,-2))`) → product \(D\) → **InvRC**.

\[
D^{J_0 J_1\lambda_0}_{ijkl}
=(-1)^{J_0+J_1+j_i+j_j+j_k+j_l}\,
c(\lambda,\lambda_0)
\sum_{ab J_2}
\begin{Bmatrix}\lambda_0&\lambda&\lambda\\ J_2&J_0&J_1\end{Bmatrix}
\mathrm{RC}\Omega^{J_1 J_2\lambda}_{jbla}\,
\mathrm{RC}\chi^{J_2 J_0\lambda}_{aibk},
\quad
W=\mathrm{InvRC}[D].
\]

**Product weight:** use \(c(\lambda,\lambda_0)\), **not** AMC’s printed \(\hat\lambda_0\).
Derived from RC–WE conjugation; \(c=0\) for odd \(\lambda_0\);
special \(c(\lambda,0)=1/\hat\lambda\). For \(\lambda=2\):
\(c\approx\{0{:}0.447,\,1{:}0,\,2{:}1.542,\,3{:}0,\,4{:}4.597\}\).

**Gold:**
| Compare | Target |
|---|---|
| RME | \(\mathrm{InvRC}[D]\equiv\mathrm{extract}(W_{1,m})\) |
| magnetic | \(\mathrm{WE}(W)\equiv P(W_1):=\mathrm{WE}\circ\mathrm{extract}(W_{1,m})\) |

Do **not** gold-compare to raw crossed \(W_{1,m}(ijkl)\) (not fully in Edmonds span).

Bench: `run/test_G4b_pathA_term1_correct.py` — **PASS**.
RC/InvRC alone were fine; only the GEMM weight needed fixing.

#### Factorized term2 (Pandya Ω × RC) — **LOCKED** λ=0 (2026-08)

Bare: \(W_2=\sum_{ab}\chi_{akbi}\Omega_{jalb}\).

Pack split of Factorized RC (`bar_CHI_V_RC`):
\[
\mathrm{RC}=\mathrm{RC}[\bar\chi]
+\mathrm{RC}[-h_Z\bar\chi^{T}]
=\mathrm{RC}_1+\mathrm{RC}_2,
\]
where \(\bar\chi^{T}\) is the Pandya-matrix transpose \(\bar\chi_{bc,ad}\).
**Do not** replace \(\mathrm{RC}_2\) by \(-h_Z(\mathrm{RC}_1)^{T}\) — finished RC does not transpose that way.

Flow (same Pandya \(\bar\Omega\) as term1):
\[
\bar W_2=\bar\Omega\cdot\mathrm{RC}_2,
\quad
Z_2=\mathrm{InvPandya}_{(1-P)^2}[\bar W_2]
\equiv(1-P_{ij})(1-P_{kl})(-W_{2,m}).
\]

Also: \(\mathrm{RC}_1+\mathrm{RC}_2=\mathrm{RC}_{\mathrm{full}}\);
term1 alone \(\bar\Omega\cdot\mathrm{RC}_1\to(1-P)^2 W_1\) (λ=0).

Bench: `run/test_G4b_term2_pandya_RC.py` — **PASS**.

#### Path A term2 via \(ijkl\leftrightarrow klij\) — **LOCKED** all λ (2026-08)

Reuse term1 χ via \(W_1\) at \(klij\) (no second GEMM). M-scheme Ω hermiticity
gives \(W_2=h_\Omega W_1|_{klij}\) at λ=0. For λ≠0 the same map holds on RMEs with
the reduced-tensor conjugation factor from the RC–WE product algebra:

\[
W_2^{J_0J_1\lambda_0}_{ijkl}
=h_\Omega(-1)^{J_0-J_1}f(\lambda,\lambda_0)\,W_1^{J_1J_0\lambda_0}_{klij},
\]

\[
f(\lambda,\lambda_0)
=\frac{\displaystyle\sum_{M_0,M_1}
\langle\lambda\,M_1,\,\lambda_0\,(M_0-M_1)|\,\lambda\,M_0\rangle}
{\displaystyle\sum_{M_0,M_1}(-1)^{M_1}
\langle\lambda\,M_1,\,\lambda_0\,(M_0-M_1)|\,\lambda\,M_0\rangle}
\quad(\lambda_0\text{ even; odd}\Rightarrow c=0).
\]

Equivalent: \(f=(-1)^\lambda c^\dagger/c\) with Path A weight \(c(\lambda,\lambda_0)\) and
\(c^\dagger\) using Ω† WE CG \(\mathrm{CG}(J_1,-M_1;\lambda,-\mu;J_2)\,J_1/J_2\)
(reference \(J_0=J_1=\lambda\), \(J_2=0\)).

Specials: \(f(\lambda,0)=(-1)^\lambda(2\lambda+1)\); \(f(\lambda,\lambda)=1\) (even λ).

Path A: \(W_1=\mathrm{InvRC}[D_1]\) at \((k,l,i,j)\) with \(J_1\leftrightarrow J_0\).

Bench: `run/test_G4b_pathA_term2_klij.py` — **PASS** λ=0..4 all even λ₀.

**Still open (do not block):** rectangular Factorized \(\mathrm{RC}_2\!\times\!\bar\Omega\!\to\!\mathrm{Inv}\);
AMC term2 Yutsis crash. Fold remains alternate gold.

#### When to apply \((1-P)\) — scalar code (`FactorizedDoubleCommutator.cc`)

| Stage | Code object | \((1-P_{ij})(1-P_{kl})\)? |
|---|---|---|
| Build \(\bar\chi^\iota\) | `bar_CHI_V = Γ̄·(occ⊙Ω̄)` | **no** |
| Inv Pandya → normal χ (gold) | Path B inv of `bar_CHI_V` | **no** (bare χ) |
| RC / pack | Factorized `bar_CHI_V_RC` (or AMC 03b + hand pack) | **no** — pack ≠ external \(P\) |
| ×Ω | `CHI_V_final = Ω̄ · RC` | **no** |
| Back to \(Z\) | Inv Pandya of `CHI_V_final` (L1884–2013) | **yes — only here** |

The relative minus \(\chi_{aibk}\Omega_{jbla}-\chi_{akbi}\Omega_{jalb}\) (analyze) is the **topology pack**, identical in role to code \(\chī_{adbc}-h_Z\chī_{bcad}\). Do **not** confuse it with bra/ket antisymmetrization.

#### Factorized Pandya→RC (locked dual)

**Full equation flow (scalar + tensor continuous):**
[`output/steps/FACTORIZED_RC_FLOW.md`](output/steps/FACTORIZED_RC_FLOW.md).

Factorized does Pandya→RC in one shot (`bar_CHI_V_RC` L1833–1841):

\[
\overline{\overline\chi}{}^{\iota\,J}_{a\bar b\,c\bar d}
=
+\sum_{J'}(2J'+1)\,(-1)^{j_b+j_c+J'}
\begin{Bmatrix} j_a & j_b & J \\ j_c & j_d & J' \end{Bmatrix}
\Bigl(
  \bar\chi^{\iota\,J'}_{a\bar d\,b\bar c}
  -h_Z\,\bar\chi^{\iota\,J'}_{b\bar c\,a\bar d}
\Bigr).
\]

**Full chain ≡ m at λ=0** via Factorized RC dual (`test_G4b_factorized_fullZ.py`,
`test_G4b_factorized_tensor.py` λ=0).
**Any λ** via Path B χ → fold (`test_G4b_pathB_fold_mscheme.py`).
RC-only angular check: `test_chi_iota_rc.py`.
Tensor RC×Ω̄×Inv: λ→0 NineJ layout fixed in `FACTORIZED_RC_FLOW.md` §3′/§5′;
rectangular full \(Z\) still **FAIL** vs m (`test_G4b_factorized_tensor.py` λ≠0) — fold remains gold.

### Which arxiv transform for \(\overline{\overline\chi}{}^{\iota}_{\mathrm{RC}}\) / \(\Gamma^{\mathrm{IV}_b}\)?

Arxiv `\eqref{cross-coupled_Pandya}` has **two** normal→CC wirings:

| | LHS | Used by |
|---|---|---|
| **CC #1** | \(\overline{\overline{A}}^{J}_{j\bar l\,k\bar i}\) | \(\Gamma^{\mathrm{III}_b}\), \(\Gamma^{\mathrm{IV}_b}\) |
| **CC #2** | \(\overline{\overline{A}}^{J}_{i\bar k\,l\bar j}\) | \(\Gamma^{\mathrm{III}_c}\), \(\Gamma^{\mathrm{IV}_c}\) (single-bar products) |

**Choose CC #1** for IV_b: arxiv writes

\[
\overline{\overline\Gamma}^{\mathrm{IV}_b\,J}_{j\bar l\,k\bar i}
=(1-P)^2\sum_{ab}
\overline{\overline\Omega}^{J}_{j\bar l\,a\bar b}
\bigl(\bar\chi^{\iota\,J}_{a\bar b\,k\bar i}-\bar\chi^{\iota\,J}_{k\bar i\,a\bar b}\bigr).
\]

Note arxiv keeps \(\chi^\iota\) as **single-bar Pandya** and puts **Ω in CC #1**.
Factorized dualizes that: keep \(\bar\Omega\) (Pandya) and build `bar_CHI_V_RC` so

\[
\bar\Omega\cdot\mathrm{RC}[\bar\chi^\iota]
\quad\text{implements the same fold as}\quad
\overline{\overline\Omega}\cdot\bar\chi^\iota
\]

in the \(j\bar l\,k\bar i\) slot. So `bar_CHI_V_RC` is **not** “apply CC #1 or #2 to χ”;
it is the IIb leg-recouple \(X^{J}_{ab,cd}\leftarrow X^{J'}_{ad,bc}\) plus ι pack, aimed at **CC #1
contraction layout**. Do **not** pick CC #2 for this diagram.

Pipeline for Path B \(\Gamma^{\mathrm{IV}_b}\) (**LOCKED any λ**):
1. Pandya Γ, Pandya Ω → `bar_CHI_V = Γ̄·(occ⊙Ω̄)`
2. InvPlus → normal χ (no AS)
3. Analyze fold \(W=\chi\Omega\) pack (topology minus) — **or** λ=0 Factorized RC→×Ω̄
4. \((1-P)^2\) only on the way to \(Z\)
Benches: `test_G4b_pathB_fold_mscheme.py` (any λ); `test_G4b_factorized_fullZ.py` (λ=0 RC dual).

### Pandya → RC (AMC angular skeleton only)

AMC schemes (user): Pandya `((1,-4),(3,-2))` → RC `((1,-3),(4,-2))`.

| Variant | Input | Output |
|---|---|---|
| **tensor** | `input/chi_iota_pandya_to_rc.txt` | `output/chi_iota_pandya_to_rc.tex` (`--collect-ninejs`) |
| **scalar** | `input/chi_iota_pandya_to_rc_scalar.txt` | `output/chi_iota_pandya_to_rc_scalar.tex` |

Angular skeleton only (`dbarχ_ijkl = χ̄_ijkl`); **not** the Factorized pack
\((\bar\chi_{adbc}-h_Z\bar\chi_{bcad})\) — that pack fails Yutsis (same as η).

**Tensor (ninej):**
\[
\overline{\overline\chi}{}^{\iota\,J_0 J_1\lambda}_{ijkl}
=
(-1)^{J_0+J_1+j_i+j_j+j_k+j_l+\lambda}\,
\hat J_0\hat J_1
\sum_{J_2 J_3}
(-1)^{J_2+J_3}\,
\hat J_2\hat J_3
\begin{Bmatrix}
\lambda & J_0 & J_1 \\
J_3 & j_k & j_j \\
J_2 & j_i & j_l
\end{Bmatrix}
\bar\chi^{\iota\,J_2 J_3\lambda}_{ijkl}.
\]

**Scalar** (bit-identical to η `03a_RC_scheme_change`):
\[
\overline{\overline\chi}{}^{\iota\,J}_{ijkl}
=
(-1)^{J+j_k+j_l}
\sum_{J'}
(-1)^{J'}\,
\hat J'^2
\begin{Bmatrix} j_l & j_j & J \\ j_k & j_i & J' \end{Bmatrix}
\bar\chi^{\iota\,J'}_{ijkl}.
\]
Tensor \(\lambda\to 0\) collapses to this.

**Vs arxiv / analyze:**

| Source | What it is | Match? |
|---|---|---|
| Arxiv `\eqref{Pandya}` | normal → Pandya | AMC `chi_iota_normal_to_pandya_scalar` (with leading −) |
| Arxiv `\eqref{cross-coupled_Pandya}` #2 | normal → CC | **≡** AMC `chi_iota_normal_to_cc_scalar` (phase/6j spot-check r=1) |
| Arxiv note | “translation via Pandya and CC” | composition, not a printed Pandya→CC formula |
| Analyze / Factorized RC | 6j **plus** pack \((\chī_{adbc}-h_Z\chī_{bcad})\) | AMC skeleton = 6j only; pack is code-side |
| factorized_latex RC print | \(\sum (-1)^{j_b+j_c+J'}\hat J'^2\{ja jb J; jc jd J'\}\bar X_{ad,bc}\) | same family as analyze (includes index rewire) |

So: scalar AMC Pandya→RC is the right angular twin of the tensor formula and of η;
arxiv’s printed transforms are normal↔{Pandya,CC}. Full Factorized RC = AMC map on
each of \(\bar\chi_{adbc}\) and \(\bar\chi_{bcad}\), then relative \(-h_Z\).

Full \(Z\) locked any λ (`test_G4b_pathB_fold_mscheme.py`); λ=0 also via Factorized RC dual.
Next: ethS λ≠0 Path B wiring (χ invPlus → fold); rectangular RC speed twin optional.

**Phase convention (AMC analyze, λ≠0):** combined-integer phases as in
`test_chi_kappa_m_vs_amc.py` (χ^β pattern). Do **not** phase each half-int with
`j2//2` separately.

### AMC inputs (`factored_GIV/input/`)

| File | Role | Status |
|---|---|---|
| `chi_kappa.txt` | direct arxiv product | AMC OK → `output/chi_kappa.tex` |
| `chi_kappa_as.txt` | AS-realigned product | AMC OK (printed overall phase differs; m AS exact) |
| `chi_kappa_via_pandya.txt` | Pandya→product→inv | AMC OK |
| `G4a_Wbra_noperm.txt` | ladder bra, \(P\) stripped | AMC OK → `W=-\chi_{ijbd}\Omega_{dbkl}` |
| `G4a_Wklij_noperm.txt` | ket via \(W_{klij}\) (same kernel) | AMC OK → transpose with \(h_\Omega\) |
| `G4a_Wket_noperm.txt` | literal \(\Omega_{ijdb}\chi\) | AMC Yutsis fail — do not use |
| `G4a_from_chi.txt` | full ladder with \(P\) | AMC Yutsis fail — use Wbra + Wklij |

Unfactored expand (legacy): `../input/G4a.txt`.

### Ladder (after χ) — \(P\) only on the result, never on \(\chi^\kappa\)

\(\chi^\kappa\) is **non-Hermitian**. Do **not** antisymmetrize / permute \(\chi\) itself.
Do **not** apply hermiticity / antihermiticity fill to \(\chi\) or \(\bar\chi\)
(\(h_\chi\) does not exist). Only \(\Omega\) and \(\Gamma\) are Hermitian/anti-Hermitian
and fermionic-AS; use that solely when reading their TBMEs / Pandya of those ops.

\((1-\hat P)\) acts on the **ladder DGEMM intermediate** \(W=-\chi\Omega\)
(and ket twin \(W_{klij}\)), not on \(\chi\) and not as a late pass on stored \(Z\).
Full gold bench: `run/test_G4a_pathB_mscheme.py` (Path B χ → \(W\) → \((1-P)\) on \(W\) ≡ m, \(Z_{\mathrm{red}}\)).

Bra and ket are two separate contractions (analyze \(\chi^\kappa\) vs \(\chi^{\kappa'}\)):

\[
\begin{aligned}
W_{ijkl}&=-\sum_{bd}\chi^\kappa_{ijbd}\,\Omega_{dbkl}
&&\text{(bra: }\chi\Omega\text{)}\\
V_{ijkl}&=-\sum_{cb}\Omega_{ijcb}\,\chi^{\kappa'}_{klcb}
&&\text{(ket: }\Omega\chi'\text{; }\chi'\neq P(\chi)\text{)}\\
\Gamma^{\mathrm{IV}_a}
&=(1-\hat P_{ij})\,W
+(1-\hat P_{kl})\,V.
\end{aligned}
\]

\(\chi^{\kappa'}\) is the \(\Gamma\Omega\) partner (Factorized \(\mathtt{Chi\_VI}\)), not \((1-P)\) or \(P\) applied to \(\chi^\kappa\).

**Ket without a second AMC string:** same bra kernel with legs swapped, then
**Hermitian** completion onto \(Z\) (final \(\Gamma^{\mathrm{IV}_a}\) is Hermitian):

\[
\begin{aligned}
W_{ijkl}&=-\sum_{bd}\chi^\kappa_{ijbd}\,\Omega_{dbkl},\\
W_{klij}&=-\sum_{bd}\chi^\kappa_{klbd}\,\Omega_{dbij},\\
\Gamma^{\mathrm{IV}_a}
&=(1-\hat P_{ij})\,W
+(1-\hat P_{kl})\,W_{\bullet\,\mathrm{as\ }klij}.
\end{aligned}
\]

That is \(Z=A+A^{T}\) in equal-\(J\) blocks. Do **not** assemble
\(V=h_\Omega W_{klij}\) into \(Z\): with \(h_\Omega=-1\) that is \(A-A^{T}\)
(anti-Hermitian). \(h_\Omega\) belongs in the Pandya build of \(\chi^{\kappa'}\)
(\(\mathtt{Chi\_VI\_II}=h_\Omega(\mathrm{occ}\odot\bar\Omega)^{T}\bar\Gamma\)),
not as a second factor when both braces reuse the same \(\chi^\kappa\).

| Step | AMC / bench | Status |
|---|---|---|
| Bra \(W=-\chi\Omega\) (no \(P\) on \(\chi\)) | `G4a_Wbra_noperm.tex` | **PASS** `test_G4a_Wbra_mscheme.py` (A) |
| \((1-P_{ij})\) on **result** \(W\) | same bench (B) | **PASS** |
| Ket \(W_{klij}\) (same kernel) | `G4a_Wklij_noperm.tex` | **PASS** `test_G4a_Wket_mscheme.py` (A) |
| Hermitian assemble \(Z=(1-P_{ij})W+(1-P_{kl})W_{klij}\) | equal-\(J\): \(A+A^{T}\) | **PASS** vs m |
| Full \(\Gamma^{\mathrm{IV}_a}\) Path B | χ Path B → \(W\) → \((1-P)\) on \(W\) | **PASS** `test_G4a_pathB_mscheme.py` |

AMC bra / \(W_{klij}\) RME (`reduce=true` → \(W_{\mathrm{red}}\)):
\[
\begin{aligned}
W^{J_0}_{ijkl}
&=
-(-1)^{J_0}\hat J_0^{-1}
\sum_{bd J_2}
(-1)^{j_b+j_d+\lambda}\hat\lambda^{-1}
\,\chi^{J_0 J_2\lambda}_{ijbd}\,\Omega^{J_2 J_0\lambda}_{dbkl},\\
W^{J_0}_{klij}
&=
-(-1)^{J_0}\hat J_0^{-1}
\sum_{bd J_2}
(-1)^{j_b+j_d+\lambda}\hat\lambda^{-1}
\,\chi^{J_0 J_2\lambda}_{klbd}\,\Omega^{J_2 J_0\lambda}_{dbij},\\
Z^{J_0}_{ijkl}
&=(1-P_{ij})W^{J_0}_{ijkl}+(1-P_{kl})W^{J_0}_{klij}.
\end{aligned}
\]
Phase: one integer `phase((j2_b+j2_d)/2 + λ)` (χ^β convention). m uses \([\chi^\lambda\times\Omega^\lambda]^{(0)}\) CG.

### Status (IV_a)

| Piece | Status |
|---|---|
| AS identity (arxiv ≡ AS) | **PASS** `test_chi_kappa_mscheme.py` |
| χ^κ AMC print | OK (`chi_kappa_analyze.tex`) |
| χ^κ m ≡ AMC J (bare) | **PASS** all λ `test_chi_kappa_m_vs_amc.py` — tensor/WE Ω |
| χ^κ Path B ≡ AMC direct | **PASS all λ** `test_chi_kappa_pathB_vs_direct.py` (rect. Pandya + `VI_II` + \((-1)^{J_0+J_1}\)) |
| bra \(W=-\chi\Omega\); \((1-P)\) on \(W\) only | **PASS** `test_G4a_Wbra_mscheme.py` |
| ket \(W_{klij}\); Hermitian \(Z=W+W_{klij}\) (not \(h_\Omega W_{klij}\)) | **PASS** (full Z vs m) |
| Full \(\Gamma^{\mathrm{IV}_a}\): Path B χ → ladder → \((1-P)\) on \(W\) ≡ m | **PASS** `test_G4a_pathB_mscheme.py` (λ=0..; \(Z_{\mathrm{red}}\)) |
| AMC via_pandya ≡ analyze | open (input still old arxiv product) |
| full m ≡ Path B ≡ ethS (λ≠0) | **PASS** `test_tts_GIVa_eths_vs_pathB.py` (rect. Path B in `comm223_232_GIVa`) |
| Do not retune to TTS Path A | Path A removed from ethS |

### Re-run AMC (IV_a)

```bash
cd learn/amc_tts/factored_GIV
amc -o output/chi_kappa.tex input/chi_kappa.txt
amc -o output/chi_kappa_as.tex input/chi_kappa_as.txt
amc -o output/chi_kappa_via_pandya.tex input/chi_kappa_via_pandya.txt
amc -o output/G4a_Wbra_noperm.tex input/G4a_Wbra_noperm.txt
amc -o output/G4a_Wklij_noperm.tex input/G4a_Wklij_noperm.txt
```

---

## Note on Factorized vs TTS

At λ=0, Factorized DGEMM ≠ AMC TTS G4 (same class of mismatch as \(f^{\mathrm{III}_a}\)). Default = Factorized; `*_slow` = TTS oracle.

## λ≠0 — what was tried / do not repeat

| Attempt | Result |
|---|---|
| Equal-\(J\) NineJ Pandya of Ω (same \(J_{\mathrm{bra}}=J_{\mathrm{ket}}=J_{\mathrm{cc}}\)) then scalar CHI_VI layout | **Wrong vs TTS** (dead stubs remain in `comm223_232_GIVa` after early return) |
| Reuse III `barred_eta_tensor` / `barCHI_III` (Ω̄×Ω̄) | **Wrong product** (T×T→S); also norm 0 for GIIIb |
| Stuff Path-A strips into `GetMatrix` | Same failure mode as GIIIa |

## λ≠0 — two viable paths (pick one and calibrate)

### Path A (TTS-faithful strips — recommended first)

G4a Term1 6j **cannot** eliminate both Ω strings into a 4-index χ without keeping mid indices \((J_6,j_0)\):

\[
\chi_{pgbd}^{J_6 j_0}
=\sum_{ac J_2 J_3 J_5}
(\mathrm{occ})\,(\text{6j without }J_4)\,
\Omega^{J_2J_3\lambda}_{apcd}\,\Gamma^{J_5}_{gcba},
\]
then contract with \(\Omega^{J_4 J_1}_{dbqh}\) and the remaining 6j that carry \(J_4\).

- Mirror `FillChiEtaG3a` / `ChiTab`, but strides include \(J_6,j_0\).
- Copy phases/hats **verbatim** from `tts_GIVa` (Term1 then Term2).
- Term2 has free \(q/h\) inside one Ω — needs its own strip layout (do not force Term1 indices).
- Validate strip vs on-the-fly TTS **before** any Pandya.

### Path B (Factorized rectangular CC — production speed)

Same topology as scalar CHI_VI, but Ω̄ is **rectangular**:

\[
\bar\Gamma^{J},\quad
\bar\Omega^{J J'\lambda}
\;\xrightarrow{\mathrm{DGEMM}}\;
\bar\chi^{\kappa\,J J'\lambda}
=\bar\Gamma^{J}\cdot(\mathrm{occ}\odot\bar\Omega^{J J'}),
\]
then **tensor** inv Pandya (NineJ; cf. `TensorCommutators::AddInverseTensorPandyaTransformation`) into rank-λ `Chi`, then mid-\(J\) ladder

\[
Z^{J}\;{+}{=}\;
-\sum_{J_m}\Omega^{J J_m}\chi^{J_m J}
-\sum_{J_m}\chi'^{J J_m}\Omega^{J_m J}.
\]

Use Factorized **2×nKets** layout consistently, **or** full `TensorCommutators` ph Pandya — do not mix. Occ weight = Factorized `occ_ABbarD` / partner for χ^{κ'}. Materialize `.t()` for χ^{κ'}.

`comm222_phst` is the closest existing DGEMM pattern (scalar × tensor in Pandya space).

## Path A notes (GIVa — calibrated)

- **χ_ex ≠ swap(p,g) of χ:** TTS exchange keeps \(\Gamma_{gcba}\) and puts \(j_p\) in the last 6j.
- **Term2 χ_T2_ex ≠ swap(q,h):** TTS keeps \(\Gamma_{bahd}\); only η leg + \(j_q\leftrightarrow j_h\) in 6j change. Layout \((q,h,b,c,J_0,J_6,j_{02})\) — **not** Term1’s \((p,g,b,d,J_6,j_{02})\).
- **Phase T1:** χ absorbs \((j_{2p}+j_{2b}+j_{2c}+j_{2d})/2 + J_2+J_3+J_5+\lambda\); ladder uses \(J_1+J_4\).
- **Phase T2:** χ absorbs \((j_{2b}+j_{2c}+j_{2d}+j_{02})/2 + J_2+\lambda\); ladder uses \(J_4\) (J0-dependent 6js stay inside χ).
- **TTS bug fixed:** Term1/Term2 must not `continue` past exchange when direct ME is tiny (same class as GIIIa LESSONS). Fixed in `tts_GIVa`.

## Path A notes (GIVb — calibrated)

- **(1-P_ij)(1-P_kl) ≠ index swaps of one χ:** four sectors pick \((i,k,\mathrm{spec})\) explicitly.
- **Spectator index:** bake \((j_{2\mathrm{spec}}+\ldots)/2\) into χ (do not split half-integers across χ/ladder).
- **Topo-1 vs topo-2:** same id layout, different ME/6j; topo-1 phase uses \(j_b\), topo-2 uses \(j_a\).
- Gate: `*_slow`→TTS; λ≠0→Path B χ→fold (locked); λ=0→Factorized CHI_V RC dual also OK.

## Path A notes (GIVc — calibrated)

- Overall \(1/2\); two topologies with distinct χ layouts (not shared with GIVb).
- Topo-1: \(\chi_1[i,k,b_p,c,d,J_2,J_3,j_{02}]\) — γ first \(i\), η₁ ket \(k\), phase bra-partner \(b_p\).
- Topo-2: \(\chi_2[i,k,b_g,b_k,d,J_2,J_3,j_{02}]\) — η₁ bra \(i\), γ ket \(k\), phase legs \((b_g,b_k)\).
- Four (1-P) sectors = four index picks; never invent by swapping one strip.

## Γ^{IV_c} / χ^λ — m / AMC / Path A (locked 2026-07-29)

χ^λ is **T×S + S×T → tensor**; fold is **T×T→S** ring.

\[
\chi^\lambda_{ijkl}
=\sum_{ab}\Bigl[
  w_l\,\Gamma_{ijab}\Omega_{abkl}
 +w_j\,\Omega_{ijab}\Gamma_{abkl}
\Bigr],
\quad
\Gamma^{\mathrm{IV}_c}
=\tfrac12(1-P_{ij})(1-P_{kl})
\sum_{ab}\chi^\lambda_{ialb}\Omega_{bjak}.
\]

### AMC inputs (`factored_GIV/`)

| File | Role |
|---|---|
| `input/chi_lambda.txt` | χ^λ unreduced Γ (no 6j; Term1+Term2) |
| `input/chi_lambda_gamma_reduced.txt` | + \(\hat J_0^{-1}/\hat J_1^{-1}\) |
| `input/G4c_from_chi.txt` | fold unreduced |
| `input/G4c_from_chi_reduced.txt` | fold Case-2 (×\(\hat J_0\)) |
| `input/G4c_direct{,_reduced}.txt` | unfactored ΩΩΓ with `P(i/j)P(k/l)` |

### Packaging (same class as χ^θ)

| Object | Formula |
|---|---|
| \(Z(m)\) | \(\tfrac12(R_1-R_2-R_3+R_4)\), \(R=\sum CG(\lambda\mu;\lambda-\mu;00)\,\chi(m)\Omega(m)\) |
| \(S\) | \(\sum CG\,CG\,Z(m)\) |
| \(Z_{\mathrm{red}}\) | \(S/\hat J\) ≡ AMC `reduce=true` fold |
| \(Z_{\mathrm{unred}}\) | \(S/\hat J^2\) ≡ AMC default fold |

**Do not** compare bare \(S\) to AMC unreduced (fake \(\hat J^{\pm2}\) ratios). χ^λ itself is always WE-reduced (tensor); AMC χ with unreduced Γ has no extra 6j.

### Status

| Piece | Status |
|---|---|
| χ^λ m ≡ AMC | **PASS** `run/test_chi_lambda_mscheme.py` |
| m ≡ ring fold (red/unred) | **PASS** exhaustive `run/test_tts_GIVc_mscheme.py` |
| Path A ≡ TTS Case-2 DIRECT | **PASS** (~1e-14) — ethS Path A strips |
| Path A / TTS ≟ χ-fold | **OPEN** — TTS DIRECT ≠ χ-fold; do not retune ring to Path A |
| ethS Path B (Pandya→DGEMM→inv) | **PASS** `use_TypeGIVc_factorized` → `run/test_tts_GIVc_pathB.py` |

### Fold angular gold (locked 2026-07-29)

Printed AMC `G4c_from_chi*_ninej.tex` is **wrong** vs m (~1–3% on e.g. `2345`; non-Hermitian under bra↔ket).

Correct kernel = locked **tts_ring** Path A (`z_pbar_aqsb_direct`) with χ≡Ω₁, Ω≡Ω₂:

\[
X_{pqsr}=\sum_{ab}\chi_{pbar}\,\Omega_{aqsb},
\qquad
Z_{ijkl}=\tfrac12\bigl(
  X_{ijkl}
  -(-1)^{j_k+j_l-J}X_{ijlk}
  -(-1)^{j_i+j_j-J}X_{jikl}
  +(-1)^{j_i+j_j+j_k+j_l}X_{jilk}
\bigr).
\]

(m-products \(\chi_{ialb}\Omega_{bjak}\equiv\chi_{ibal}\Omega_{ajkb}\).)

**Path B speed form (λ≠0):** IMSRG tensor Pandya (`adcb`) with map
AMC \(\bar\chi(p,b,a,r)=\) IMSRG\((p,r,a,b)\), mid-J DGEMM
\(\bar X^{J}\mathrel{+}=\lambdâ^{-1}(-1)^{J'+\lambda}\bar\chi^{JJ'}\bar\Omega^{J'J}\),
corrected inv (drop AMC-sample minus), then fermionic AS. λ=0 stays Factorized CHI_VII.

### Path B AMC (starting point for ethS)

| File | Role |
|---|---|
| `input/G4c_chi_omega_pandya.txt` | fwd Pandya(χ,Ω) → product → inv; **Z reduce=true** |
| `input/G4c_chi_omega_pandya_unred.txt` | same, unreduced Z |
| `input/G4c_chi_omega_pandya_invP.txt` | product without P; AS on inv |
| `output/G4c_chi_omega_pandya*_ninej.tex` | printed RME |

**Reduce flags:** χ/Ω tensors always reduced; `barG4c`/`G4c` use `reduce=true` for \(Z_{\mathrm{red}}=S/\hat J\).

**Locked composition rules** (same as `tts_ring` / `omega_cross`):

1. **Inv:** `G4c = barG4c` — drop AMC-sample overall minus so Path B ≡ Path A ≡ m.
2. **Fwd Pandya (printed AMC):** evaluate RHS \(O_{ijkl}\) with **same-label** ME. Do **not** silently apply IMSRG `adcb` remap when coding the printed ninej. When wiring ethS via `DoTensorPandyaTransformation`, map indices explicitly (see `comm222_phst/AMC_CHECK.md`).
3. Gold: Path B ≡ Path A (`G4c_from_chi_reduced`) ≡ m packaging.

```bash
amc --collect-ninejs -o output/G4c_chi_omega_pandya_ninej.tex input/G4c_chi_omega_pandya.txt
amc --collect-ninejs -o output/G4c_chi_omega_pandya_unred_ninej.tex input/G4c_chi_omega_pandya_unred.txt
```

### Re-run AMC

```bash
cd learn/amc_tts/factored_GIV
amc -o output/chi_lambda.tex input/chi_lambda.txt
amc -o output/chi_lambda_gamma_reduced.tex input/chi_lambda_gamma_reduced.txt
amc --collect-ninejs -o output/G4c_from_chi_ninej.tex input/G4c_from_chi.txt
amc --collect-ninejs -o output/G4c_from_chi_reduced_ninej.tex input/G4c_from_chi_reduced.txt
amc -o output/G4c_direct.tex input/G4c_direct.txt
amc -o output/G4c_direct_reduced.tex input/G4c_direct_reduced.txt
amc --collect-ninejs -o output/G4c_chi_omega_pandya_ninej.tex input/G4c_chi_omega_pandya.txt
```

## Speed TODO (λ≠0 only)

1. ~~GIVa / GIVb / GIVc Path A PASS.~~
2. ~~**GIVc Path B**~~ — λ≠0 Pandya→DGEMM→inv (`use_TypeGIVc_factorized`); λ=0 = Factorized CHI_VII.
3. ~~III λ≠0 mid-J~~ — done for III_a/b/c.
4. GIVc: close Path A/TTS ≟ AMC from_chi (DIRECT vs χ-fold).
