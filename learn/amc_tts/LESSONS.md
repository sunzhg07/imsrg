# TTS / ethS lessons (agent reference)

**Read this first** before implementing or optimizing any new TTS / ethS diagram.
Findings from debugging factorized vs unfactored TTS
(\(f^{\mathrm{I\ldots III}}\), \(\Gamma^{\mathrm{I\ldots IV}}\)), χ storage, and Pandya→DGEMM ports.

| Doc | Role |
|---|---|
| **This file** | Cross-cutting rules + playbook (recall every time) |
| **[tensor_pro_final/](tensor_pro_final/)** | **Curated gold:** m / AMC direct / Path B / benches per diagram |
| **[F_1B_PROCEDURE.md](F_1B_PROCEDURE.md)** | **\(f^{\mathrm{I\ldots III}}\): m → AMC direct → Path B → DGEMM** |
| **[REDUCED_UNREDUCED.md](REDUCED_UNREDUCED.md)** | **CRITICAL — packaging before any m↔J compare** |
| [factored_DGEMM_STATUS.md](factored_DGEMM_STATUS.md) | Which pieces are Factorized DGEMM vs Path A / TTS |
| [FACTORIZED_TTS_IMPLEMENTED.md](FACTORIZED_TTS_IMPLEMENTED.md) | Piece checklist + flags + tests |
| `factored_*/NOTES.md` | Per-diagram detail |

---

## CRITICAL: reduced vs unreduced (do not skip)

**User has emphasized this repeatedly. Agents still get it wrong. Lock packaging before numerics.**

Full note: **[REDUCED_UNREDUCED.md](REDUCED_UNREDUCED.md)**.

- AMC `reduce=true` / `IsReduced()=true` → **reduced** \(X_{\mathrm{red}}\) (often overall \(\hat J\) in the printout).
- m-scheme MEs are **always physical** (unreduced). There is no reduced m.
- Compare: \(X^{J}_{\mathrm{red}}=\hat J^{-1}\sum_m\mathrm{CG}\,\mathrm{CG}\,X(m)\) vs AMC reduced formula.
- If ratios look like \(1/\hat J\) → **packaging bug**, not a bad 6j / “AMC is wrong”.
- Gold benches: `run/test_z_ring_mscheme_sign.py`, `run/test_omega_cross_mscheme.py`.

---

## Playbook: implementing / optimizing a new diagram

Follow this order. Skipping steps is how we burned days re-solving the same bugs.

### 0. Classify the product (before any code)

| Product | χ rank | Downstream |
|---|---|---|
| \(T\times T\to S\) (III: η, θ) | **scalar** χ | Scalar Factorized path (Pandya / RC / inv / DGEMM) |
| \(T\times S\to T\) or \(S\times T\to T\) (IV: ι, κ, λ) | **tensor** χ rank λ | Need mid-\(J\) tensor χ; λ=0 can copy scalar Factorized |
| Never | \(T\times T\to T\) | Do not invent \([\Omega\times\Omega]^{\lambda\neq 0}\) |

Map AMC name → Factorized label (`NOTES.md` table): e.g. \(\Gamma^{\mathrm{III}_a}\)=IIa/IIc, \(\Gamma^{\mathrm{III}_b}\)=IIb/IId, \(\Gamma^{\mathrm{III}_c}\)=IIe/IIf, IV↔IIIc…IIIf.

### 1. Get a trusted slow oracle first

1. Unfactored TTS / Path A / AMC dense table that **PASS**es vs `ReferenceImplementations.comm223_*_tts_*`.
2. If χ ≠ TTS, **audit the reference** before rewriting factorization (exchange phase, independent exchange, wrong \(j\) in 6j) — see §\(\Gamma^{\mathrm{III}_a}\).
3. Wire `use_Type*_slow` (or Path A gate) so tests can force the oracle.

### 2. Port λ=0 as a literal Factorized copy

- Copy from `FactorizedDoubleCommutator.cc` (same Pandya fill, DGEMM, inv, phases, √2).
- Template: standalone `comm223_232_GIIIb` / `GIIIa` / `GIVa` in ethS.
- **Do not** edit scalar `FactorizedDoubleCommutator.cc` for TTS experiments.
- χ from \(T\times T\) must live in **scalar** storage (`Chi_III_Op = Z.TwoBody`, never copy Eta’s tensor rank).

### 3. Dual-oracle pattern (production vs TTS)

At **λ=0**, Factorized DGEMM often **≠** AMC TTS even when both are “correct” for their partition (seen: \(f^{\mathrm{III}_a}\), GIV, GIIIa, GIIIc).

| Role | Default | Flag |
|---|---|---|
| Match **scalar Factorized** (production speed) | Factorized Pandya→DGEMM | `use_Type*_2b` |
| Match **TTS / AMC** (tests) | Path A or TTS | `SetUse_Type*_slow(True)` |

- Tests that claim TTS PASS **must** set `*_slow=True` (or use λ≠0 Path A).
- Do **not** “fix” Factorized to equal TTS by random pref factors — decide which oracle is production truth first.

### 4. λ≠0 only after λ=0 Factorized is wired

| χ type | λ≠0 strategy |
|---|---|
| Scalar χ (\(T\times T\to S\)) | Mid-\(J\) \([\bar\Omega\times\bar\Omega]^0\) → **same** RC/inv/DGEMM as λ=0 |
| Tensor χ (\(T\times S\to T\)) | Mid-\(J\) / rectangular CC for rank-λ χ → Factorized-style DGEMM |

**GIIIb mid-\(J\) (done):** Factorized continuous
\(\bar\chi^{\eta\,0\,J}=\hat\lambda^{-1}\sum_{J'}(-1)^{J+J'+\lambda}\bar\Omega^{JJ'}\cdot(\mathrm{occ}\odot\bar\Omega^{J'J})\)
in Factorized \(adcb\) layout — see [factored_GIIIb_midJ_IDENTITY.md](factored_GIIIb_midJ_IDENTITY.md).
λ=0 recovers `bar_Eta * nnnbar_Eta` bit-exact; λ≠0 default is Factorized (may ≠ TTS);
`SetUse_TypeGIIIb_slow(True)` → TTS. The old `chi2b` `barred_eta_tensor` path remains broken (wrong legs / empty mats) — do not re-enable it.

### 5. Path A vs Path B (χ storage)

| Path | When | Rule |
|---|---|---|
| **A** | TTS match, open-index AMC strips | Dense `ChiTab` (AMC index order); explicit \((1-P)\) ladder |
| **B** | Factorized speed | Pandya → inv into **normalized** `TwoBodyME` → DGEMM |

**Never** stuff Path-A Term strips into `TwoBodyME::GetMatrix` and expect DGEMM to match (√2 / ket order / exchange). Failed repeatedly for GIIIa.

### 6. Engineering checklist (every port)

| Item | Rule |
|---|---|
| SixJ | `PreCalculateSixJ()` **before** OMP; never fill SixJ cache inside parallel (`DANGER!!!!!!! Updating SixJList…` → segfault). Mid-\(J\) nests: prefer `AngMom::SixJ` or serial build |
| Arma `.t()` | Materialize `arma::mat T = M.t();` then `T * N` — chained `(M.t())*N` was unstable in Factorized CHI_VI_II |
| Hermiticity | Occupancy-weighted χ: `SetNonHermitian` + full matrix; no upper-triangle+conjugate |
| Tensor Ω | Never raw `.t()` of tensor Ω (bra/ket \(J\) flip) — need \(h_\Omega(-1)^{J_0\pm J_1}\) (see §Reduced tensor transpose) |
| χ_ex ≠ orbit swap | IV exchanges often keep Γ on one leg and only retarget η / 6j — build a **separate** strip |
| Phase + half-integers | Do not split \(j_a+j_b+\ldots\) across χ/ladder with per-orbit `j2/2` ints; use one \(( \sum j_2 )/2\) or keep the spectator \(j\) inside χ (GIVb) |
| TTS `continue` | Never `continue` past an exchange when the direct ME is tiny — fix the reference first |
| Shared \(J\) windows | Two swapped Γ/Ω strings → **union** of \(J\) windows (or loop all \(J\)); never one window from the first string only |
| Isolate pieces | Standalone `comm223_232_G*` + `run/test_tts_G*.py`; don’t only bench monolithic `chi2b` |
| Bench | Always report λ=0 Fact vs TTS **and** λ≠0 PathA/TTS; emax=1 He4 is the smoke size |

### 7. Topology → DGEMM target

| Topology | Target |
|---|---|
| Ladder (\(\chi\Gamma\), \(\Gamma\chi\)) | Ordinary-channel **DGEMM** |
| Not a ladder | **Pandya/CC** so it becomes a ladder in \(\bar{}\) space → **DGEMM** → inverse if needed |

Maximize DGEMM; Pandya only to *make* a ladder.

### 7b. Do NOT inverse-Pandya χ if the fold is also a CC contraction

**Read the scalar template before porting.** For a **one-body** final, the
scalar Factorized code often keeps χ **inside** the Pandya representation and
contracts it with \(\bar\Gamma\) there — the inverse Pandya never happens:

```
IntermediateTwobody[ch_cc] = (2J_cc+1) * Xbar * Xbar_nnnn * Gamma_bar
f_pq = ĵ_p^{-2} Σ_{e,J_cc} Intermediate[(pe),(qe)] − partner
```

(`FactorizedDoubleCommutator.cc` `comm223_231_chi2b` §II_a/II_c, Eq. B5c.)

Inserting an inverse Pandya + `MakeNotReduced` + dense `chi_tab` + ordinary
ladder is **not** a harmless detour:

- it manufactures a reduced/unreduced question (`chi = ±barChi`, "drop the
  sample minus") that the CC route never has to answer;
- it throws away every DGEMM (\(f^{\mathrm{III}_a}\): ~9× slower at emax=4).

**Rule:** before writing any \(\bar\chi\to\chi\) step, check whether the fold
partner (\(\Gamma\) here) can be Pandya-transformed instead. If yes, stay in CC.
Burned a long time on \(f^{\mathrm{III}_a}\) this way — see
[factored_fIIIa/NOTES.md](factored_fIIIa/NOTES.md) §Path B.

**Tensor Ω inside a CC χ:** use the plain IMSRG tensor Pandya (`hat` only, no
per-element `scale`) on the **reduced** Ω and put the packaging in the mid-\(J\)
prefactor:
\[
\bar\chi^{J}=\frac{\hat\lambda^{-1}}{2J+1}\sum_{J'}(-1)^{J+J'+\lambda}
\bar\Omega^{JJ'}\bigl(\mathrm{occ}\odot\bar\Omega^{J'J}\bigr).
\]
\(\lambda=0\) then reproduces the scalar `Xbar * Xbar_nnnn` **bit-exact**
(because \(\bar\Omega^{JJ}(\Omega_{\mathrm{red}})=\hat J\,\bar\Omega_{\mathrm{scalar}}\)),
so one code path covers all λ and Ω never needs `MakeNotReduced`.
**Anchor every new tensor CC χ on that λ=0 identity** instead of tuning prefactors.

### Naming map (AMC ↔ Factorized)

| AMC / ethS | Factorized | χ product |
|---|---|---|
| \(\Gamma^{\mathrm{III}_a}\) | IIa, IIc | \(T\times T\to S\) → inv → ladder DGEMM |
| \(\Gamma^{\mathrm{III}_b}\) | IIb, IId | same χ̄ → RC → DGEMM → inv |
| \(\Gamma^{\mathrm{III}_c}\) | IIe, IIf | ordinary χ^θ → Pandya → ×Γ̄ → inv (+½) |
| \(\Gamma^{\mathrm{IV}_a}\) | IIIc, IIId (CHI_VI) | \(T\times S\to T\) |
| \(\Gamma^{\mathrm{IV}_b}\) | IIIa, IIIb (CHI_V) | \(T\times S\to T\) |
| \(\Gamma^{\mathrm{IV}_c}\) | IIIe, IIIf (CHI_VII) | \(T\times S\to T\) |

---

## Default coupling rule (TTS / ethS)

Product ranks (operator \(J\)-rank):

| × | scalar | tensor |
|---|---|---|
| **scalar** | scalar | tensor |
| **tensor** | tensor | **scalar** |

- **Never** \(\mathrm{tensor}\times\mathrm{tensor}\to\mathrm{tensor}\) (no \([\Omega\times\Omega]^{\lambda\neq 0}\)).
- \(\mathrm{tensor}\times\mathrm{tensor}\to\mathrm{scalar}\): χ from Ω×Ω (\(\chi^\eta,\chi^\theta,\ldots\)) has \(\lambda_\chi=0\); λ of Ω enters only while building that scalar χ (Pandya / mid-\(J\)). Downstream: **scalar** inv-Pandya, DGEMM, RC — same as `FactorizedDoubleCommutator.cc`.
- \(\mathrm{tensor}\times\mathrm{scalar}\to\mathrm{tensor}\) (and \(S\times T\to T\)): IV family χ^{ι,κ,λ} inherit rank λ from Ω; do **not** treat them as scalar CHI_III-style intermediates.

## Trust boundary

- **Trusted M-scheme:** `amc/examples/sample_output/factorized_code_analyze.tex` §unfact only.
- **Do not trust:** `amc/.../reference_tensor.*`, `comm223_*_tts_BruteForce` derived from them.
- **AMC:** declare \(\Omega\) `scalar=false`, \(\Gamma\) and final \(f\)/`Z` `scalar=true` (forces \(j_i=j_j\)).
- Tensor AMC overall rationals (\(\pm\tfrac12\), \(\pm1\)) can be wrong; keep analyze overall signs; take 6j/hat/phase from AMC.

## Scalar selection rules (TTS)

| Object | Constraint | Comment |
|---|---|---|
| \(\chi^\alpha\) (scalar intermediate) | \(j_d = j_e\) | AMC \(\delta_{j_e,j_d}\) only — **not** full `OneBodyChannel(l,j,tz)` |
| \(f^{\mathrm{I}}\) (scalar final) | \(j_p = j_q\) | AMC \(\delta_{j_j,j_i}\) |
| Same-\(j\), different-\(l\) | allowed for \(\chi\) pairs | Restricting to `GetOneBodyChannel` drops terms for tensor \(\Omega\) |

## Critical bug: shared \(J\) window for two \(\Gamma\) strings (\(f^{\mathrm{I}}\))

**Symptom:** factorized \(f^{\mathrm{I}}\) disagreed with unfactorized `comm223_231_tts_fI` by a large factor (~3× too small). Exact copy of the unfactored loops inside ethS matched → algebra OK; factorization loop was wrong.

**Cause:** when contracting
\[
\chi_{de}\,(\Gamma_{epdq}+\Gamma_{depq}),
\]
code used **one** triangular \(J\) bound taken from the first string only:
\[
J\in\big[|j_e-j_p|,j_e+j_p\big]\cap\big[|j_d-j_q|,j_d+j_q\big].
\]
The second string \(\Gamma_{depq}\) needs the \(d\leftrightarrow e\) window, which can differ. Shared bounds **drop valid \(J\)** for the second string.

**Fix (match unfactorized reference):**
- Loop \(J=0\ldots J_{\max}\) (or take the **union** of both strings’ windows).
- Let `GetTBME_J` return 0 when forbidden.
- Do **not** reuse a single `Jmin`/`Jmax` for both \(\Gamma\) orderings.

**Check:** `run/test_tts_fI.py` — factorized ethS `comm223_231_chi1b_tensor` / `comm223_231_st` (1b only) vs `ReferenceImplementations.comm223_231_tts_fI`.

**Where fixed:** `src/FactorizedDoubleCommutator_eths.cc` → `comm223_231_chi1b_tensor`.

## Factorization bookkeeping for \(f^{\mathrm{I}}\)

Unfactored (AMC `f1.tex`):
\[
Z_{pq}\leftarrow \tfrac12\hat\jmath_p^{-2}\sum
(\ldots)\,
\hat J^2\,\hat\jmath_d^{-2}\,\hat\lambda^{-1}\,
\Omega\,\Omega\,(\Gamma_{epdq}+\Gamma_{depq}).
\]

Factorized:
\[
\chi_{de}=\sum_{abcJ_0J_1}
(\mathrm{occ})\,(-1)^{J_0+J_1+\lambda}\,\hat\lambda^{-1}\,
\Omega^{J_0J_1\lambda}_{cdab}\Omega^{J_1J_0\lambda}_{abce}
\big/\hat\jmath_d^{2},
\]
\[
Z_{pq}\leftarrow \tfrac12\hat\jmath_p^{-2}\sum_{deJ}(2J+1)\,
\chi_{de}\,(\Gamma^{J}_{epdq}+\Gamma^{J}_{depq}).
\]

Put the overall \(\tfrac12\) in the \(\chi\Gamma\) contraction (same as unfactored), not twice.

## Tensor \(\chi^\beta\) / \(f^{\mathrm{II}}\) (tensor\(\times\)tensor\(\to\)scalar)

- Declare \(\chi^\beta\) `scalar=false`, final \(f\) `scalar=true`.
- No \(j_d=j_e\) on \(\chi^\beta\); triangle with \(\lambda\).
- Put overall \(\tfrac12\) in \(\chi^\beta\) (matches unfactored); contraction has \(\hat j_i^{-2}\) and \(\hat\lambda^{-1}\) only.
- Assemble \(f=f_a+h_\Gamma f_b\) with \(\Omega_{eidj}\) / \(\Omega_{ejdi}\) (not \(\Omega_{diej}\)).
- Loop \(J_3,J_4\) independently for both \(\Omega\) strings (same shared-window lesson as \(f^{\mathrm{I}}\)).
- Checks: `run/test_tts_fII.py` vs `comm223_231_tts_fII`; docs in `factored_fII/NOTES.md`.

## Reduced tensor transpose / conjugation phase (Γ^{IV_a} lesson)

**Burned here:** Factorized CHI_VI_II is \(h_\Omega(\mathrm{occ}\odot\bar\Omega)^{T}\bar\Gamma\).
At equal \(J\) a bare transpose + \(h_\Omega\) is enough. For **rectangular** rank-λ
reduced MEs, conjugation under \(J_0\leftrightarrow J_1\) costs an extra phase.

\[
\langle J_0\|O^\lambda\|J_1\rangle
\;\longleftrightarrow\;
h_O\,(-1)^{J_0-J_1}\,
\langle J_1\|O^\lambda\|J_0\rangle
\quad\text{(real)}.
\]

Integer \(J\): \((-1)^{J_0-J_1}=(-1)^{J_0+J_1}\). Documented in `comm222_phst` as
`flipphaseY = hY * phase(Jbra - Jket)`.

Locked χ^κ VI_II DGEMM:

\[
\bar\chi^{J_0 J_1}
=
h_\Omega\,(-1)^{J_0+J_1}
\sum_{ab}\mathrm{occ}\,
\bar\Omega^{J_1 J_0}(ab;il)\,
\bar\Gamma^{J_1}(ab;kj).
\]

| Case | Extra phase |
|---|---|
| λ=0 (\(J_0=J_1\)) | \(+1\) — only \(h_\Omega\) |
| λ≠0 | \((-1)^{J_0+J_1}\) **required** |

**One Path B for any λ** — rectangular formula; λ=0 is the equal-\(J\) limit.
Do not maintain a separate “scalar” algorithm. Ω stays **WE-reduced for all λ**
(λ=0 ≠ unreduced scalar on this path). Full write-up: [factored_GIV/NOTES.md](factored_GIV/NOTES.md)
§Physics / debugging lessons.

Also: \((1-P)\) on ladder intermediate \(W=-\chi\Omega\), **never** on non-Hermitian χ;
inv without AMC printed leading − (IMSRG intentional); Hermitian \(Z=W+W_{klij}\)
without a second \(h_\Omega\) on \(W_{klij}\).

## Non-Hermitian intermediates / storage (future reference)

Applies to **any** occupancy-weighted \(\chi\) built as \(O\,(\mathrm{occ})\,O\) (or Pandya
analogue), including \(\chi^\delta\) (\(f^{\mathrm{III}_b}\)) and \(\chi^\gamma\) (\(f^{\mathrm{III}_a}\)).

| Rule | Why |
|---|---|
| `SetNonHermitian()` + fill **full** `(ibra,iket)` | Occ flips sign under hermiticity; upper-triangle + conjugate is **wrong** |
| Prefer \(\chi(ch,ch)\) squares only when topology is pp | Off-diagonal \(\chi(ch_b,ch_k)\) not needed for same-\(J\) ladders |
| Never raw `.t()` of **tensor** \(\Omega\) | Bra/ket \(J\) flip; use \(h_\Omega(-1)^{J_0-J_1}\) (= \((-1)^{J_0+J_1}\) for integer \(J\)); see §Reduced tensor transpose / Γ^{IV_a} |
| Scalar products: `.t()` only on **square** same-channel \(M\) | Unequal-channel: explicit \(\Gamma\chi\) (or \(\chi\Gamma\)) into that `MatEl` key |
| Physical → normalized `/SQRT2` when writing `GetMatrix` | Matches `GetTBME` conventions |

### \(f^{\mathrm{III}_a}\) / AMC \(\chi^\gamma\) vs `TwoBodyME`

AMC \(\chi^{\gamma}_{ijkl}\) is the **unsymmetrized J-coupled** formula. Round-tripping through
`TwoBodyME::GetTBME` (antisymmetrized ket order + exchange phases) **does not** reproduce
on-the-fly AMC values for arbitrary index order — benchmark failed until we stored a dense
table keyed by AMC indices \((i,j,k,l,J)\).

- \(\chi^\delta\): `run/test_tts_fIIIb.py`, `factored_fIIIb/NOTES.md` (pp ladder; `GetMatrix` OK)
- \(\chi^\gamma\): `run/test_tts_fIIIa.py`, `factored_fIIIa/NOTES.md` (scalar χ; AMC via \(W_1/W_2\) factorization ~10⁴× vs unfactored at emax=2; Factorized B4c at λ=0)

## \(\Gamma^{\mathrm{III}_a}\) / \(\chi^\eta\) (232) — hard-won lessons

**Status:**
- **Path A** (λ≠0 or `*_slow`): Term strips + ladder **PASS** vs `tts_GIIIa` (`run/test_tts_GIIIa.py`).
- **Path B** (λ=0 default): Factorized IIa/IIc Pandya→inv→DGEMM (~1 ms). Matches scalar Factorized; **can ≠ TTS** (dual-oracle).
**Code:** `comm223_232_GIIIa`. **Naming:** AMC \(\Gamma^{\mathrm{III}_a}\) = Factorized IIa/IIc.

### What worked (Path A — TTS oracle)

1. Build **two** scalar open-index strips from G3a AMC:
   - Term1 / brace1 → \(\chi_L(i,j,c,b)\)
   - Term2 / brace2 → \(\chi_R(b,d,k,l)\)
2. Contract with \(\Gamma\) by an explicit \((1-P_{ij})(1-P_{kl})\) ladder (not \(f^{\mathrm{III}_b}\)’s `M*=2` / `0.25` / \(A+A^{T}\)).
3. Store strips in a dense `ChiTab` (AMC index order), same rationale as \(\chi^\gamma\).

Do **not** merge Term1 and Term2 into one generic AMC \(\chi^\eta_{ijkl}\): Ω routing differs; earlier single-χ attempts failed vs TTS.

### What worked (Path B — λ=0 Factorized)

Literal copy of Factorized IIa/IIc: Pandya \(\bar\Omega\) → \(\bar\chi^\eta=\bar\Omega\cdot(n\bar\Omega)\) → **scalar** inv-Pandya into `Chi_III_Op` (=`Z.TwoBody`) → `Chi*Γ + hZ*Γ*Chi.t()` DGEMM.

### What failed (do not retry blindly)

| Attempt | Why it broke |
|---|---|
| Load Term-strip χ into `TwoBodyME` + channel DGEMM | Open-index AMC order ≠ normalized ket basis / √2 / exchange phases |
| Single shared \(\chi^\eta\) for both braces | Wrong Ω legs vs G3a Term1/Term2 |
| Tensor reverse-Pandya of χ | χ is scalar (\([\Omega\times\Omega]^0\)); λ only in building χ |
| Copy \(f^{\mathrm{III}_b}\) contraction prefs | Different diagram; wrong overall factors |
| Force Fact = TTS at λ=0 with ad-hoc prefs | Dual-oracle; use `*_slow` instead |

### TTS reference bugs uncovered while matching χ

Unfactored `comm223_232_tts_GIIIa` / `tts_GIVa` / `gen_tts_trusted.py` had bugs that made “χ vs TTS” look like a factorization error. **Fix the reference first** when strip χ disagrees with on-the-fly TTS:

1. **Exchange phase:** recompute \((-1)^{j_p}\) / \((-1)^{j_q}\) on the exchanged string — do **not** reuse the direct string’s phase.
2. **Independent exchange:** if the direct ME is tiny, do **not** `continue` past the exchange; the exchange is a separate contribution.
3. **Term2 exchange 6j:** real \(q\leftrightarrow h\) swap uses \(j_q\), not \(j_b\).

After these fixes, Path-A χ + ladder matched TTS to \(\sim10^{-14}\).

### Compaction that is safe

- One `FillChiEtaG3a(side, …)` for both Term strips (shared `ChiTab` strides).
- Further speed: triangle-bound multipole loops; fill only \(J\) that appear in `Z.MatEl`.
- **Not** safe: dropping Term2, or collapsing to one TwoBodyME for both sides (left/right χ are different operators).

### Sibling status (232 III / IV)

| Piece | λ=0 default | λ≠0 | Notes |
|---|---|---|---|
| \(\Gamma^{\mathrm{III}_a}\) | Factorized IIa/IIc DGEMM | Path A | Fact ≠ TTS @ λ=0; `*_slow`→Path A |
| \(\Gamma^{\mathrm{III}_b}\) | Factorized IIb RC→DGEMM | TTS | mid-\(J\) barCHI (chi2b) still **norm 0** |
| \(\Gamma^{\mathrm{III}_c}\) | Factorized IIe/IIf DGEMM | Path A | Fact ≠ TTS @ λ=0; two occ tables (k/j) |
| \(\Gamma^{\mathrm{IV}_a}\) | Factorized CHI_VI | Path A Term1+Term2 | χ_ex ≠ orbit swap; TTS continue fixed |
| \(\Gamma^{\mathrm{IV}_b}\) | Factorized CHI_V | Path A topo-1+2 | spectator index for phase; ~1.3 s vs 74 s |
| \(\Gamma^{\mathrm{IV}_c}\) | Factorized CHI_VII | Path A topo-1+2 | overall ½; ~0.08 s vs 17 s |

Detail: [factored_GIIIa/NOTES.md](factored_GIIIa/NOTES.md), [factored_DGEMM_STATUS.md](factored_DGEMM_STATUS.md), [factored_GIV/NOTES.md](factored_GIV/NOTES.md).

## \(\Gamma^{\mathrm{III}_c}\) / \(\chi^\theta\) (232)

**Reuse from prior diagrams:**
- Mid-\(J\) ladder product like \(\chi^\delta\): \((-1)^{J+J_2+\lambda}/\hat\lambda\,\Omega^{JJ_2}\Omega^{J_2J}\) (physical `GetTBME`).
- Path A like \(\Gamma^{\mathrm{III}_a}\): dense `ChiTab`, **do not** force Factorized IIe until Path A trusted.
- Split Term1 / Term2 (G3c): occ on \(k\) vs occ on \(j\) — one combined analyze \(\chi^\theta\) fed into IIe disagreed; antiherm remapping \(\chi_j\leftrightarrow\chi_k\) also failed (bra-exchange phases).

**Status now:** default = Factorized IIe/IIf at **all λ** (λ≠0: mid-J ordinary χ^θ DGEMM → scalar Pandya → ×Γ̄ → inv +½). `*_slow` = Path A (= TTS). Same dual-oracle as III_a/b.

**What failed historically:** feeding Path-A `chi_k`/`chi_j` into IIe; ethS CHI_IV + scalar Pandya with ad-hoc \(\hat J^{-2}\) / `MakeNotReduced` prefs before Path A was trusted.

## \(\Gamma^{\mathrm{I}}\) / \(\chi^\varepsilon\) (232)

Same scalar-selection lesson as \(\chi^\alpha\): only \(j_p=j_q\).

**Critical bug — swapped \(J\) windows in \(\chi^\varepsilon\) (slow path):**
- \(o_1=\Omega^{J_0 J_1}(c p; a b)\), \(o_2=\Omega^{J_1 J_0}(a b; c q)\)
- \(J_0\) must couple \((c,p)\)/\((c,q)\); \(J_1\) must couple \((a,b)\)

**DGEMM (default):** ordinary-channel normalized `GetMatrix` (not Pandya; not a manual \(a\neq b\) weight).  
Same convention as scalar χ^α: leading **`2 * ang * Left * W * Right`**, then outer-leg \(N_{cp}N_{cq}=\sqrt{1+\delta}\) on extract.  
MatEl only has allowed \((J,\pi,T_z)\) with \(ch_bra\le ch_ket\); missing blocks must be skipped (`GetMatrix.at` throws).

**Write path:** tensor Γ^I uses `GetTBME_J` + √2; `GetTBME_norm`/`GetLocalIndex` misses same-\(j\) different-\((l,t_z)\).

**Γ^II @ λ=0:** same DGEMM with \(\Gamma W \Omega\); do not bench against `tts_GII`.

**Driver:** `comm223_232` must **not** `MakeNotReduced` tensor \(\Omega\).

**Flags:** `use_TypeGI_slow` forces orbit/`J` loops (debug).

## \(\Omega\times\Gamma\to\) tensor χ (IV family)

**Rank rule:** \(T\times S\to T\) (and \(S\times T\to T\)). So χ^{ι,κ,λ} are **rank-λ**, unlike III’s \(T\times T\to S\).

Do **not** force scalar Factorized CHI_VI / CHI_V / CHI_VII layouts when \(\lambda_\Omega\neq 0\) without mid-\(J\) tensor χ.

| Piece | Intermediate | ethS | Bench |
|---|---|---|---|
| \(\Gamma^{\mathrm{IV}_a}\) | \(\chi^\kappa\) λ=0 Factorized DGEMM; λ≠0 → TTS | `comm223_232_GIVa` | `test_tts_GIVa.py` |
| \(\Gamma^{\mathrm{IV}_b}\) | \(\chi^\iota\) λ=0 Pandya→RC→DGEMM; λ≠0 → TTS | `comm223_232_GIVb` | `test_tts_GIVb.py` |
| \(\Gamma^{\mathrm{IV}_c}\) | \(\chi^\lambda\) λ=0 CHI_VII DGEMM; λ≠0 → TTS | `comm223_232_GIVc` | `test_tts_GIVc.py` |

Detail: [factored_GIV/NOTES.md](factored_GIV/NOTES.md).

## Debugging recipe

1. Paste unfactored loops into the factorized function → must match reference (validates wiring).
2. Re-enable \(\chi\) factorization with **identical** \(J\) / orbit loops as the unfactored reference.
3. Only then tighten \(J\) bounds (union of windows, never a single shared window for swapped legs).
4. For DGEMM: compare `*_slow` / Path A vs default Factorized; expect possible λ=0 Fact≠TTS (dual-oracle).
5. For \(\Gamma^{\mathrm{III}_a}\): if strip χ ≠ TTS, **audit the reference** (exchange phase / independent exchange / Term2 \(j_q\)) before rewriting factorization.
6. Before mid-\(J\) λ≠0: confirm λ=0 Factorized is wired; compare mid-\(J\) \(\bar\chi\) to Path-A χ **elementwise** — do not enable if norm is 0.
7. SixJ: if you see `DANGER!!!!!!! Updating SixJList inside a parellel loop`, you filled the cache under OMP — fix with `PreCalculateSixJ` or `AngMom::SixJ` / serial.

**New diagram?** Start at §Playbook at the top of this file.
