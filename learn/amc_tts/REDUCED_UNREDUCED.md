# CRITICAL: reduced vs unreduced packaging

**Read before any m-scheme ↔ AMC / J-scheme numeric compare.**
This mistake has burned days more than once. Do not “remember later” — lock packaging first.

---

## Hard rule

**Never compare numbers until you have named both sides’ packaging.**

For every scalar (or intermediate) \(X\), write explicitly:

| Side | What it is |
|---|---|
| AMC `reduce=true` / IMSRG `IsReduced()=true` | **reduced** RME \(X_{\mathrm{red}}\) |
| AMC default scalar / `IsReduced()=false` | **unreduced** \(X_{\mathrm{unred}}\) |
| `GetMschemeMatrixElement_*` | **always physical / unreduced m** (no “reduced m”) |

Relation (scalar channel \(J\)):

\[
X^{J}_{\mathrm{red}}
=
\hat J\, X^{J}_{\mathrm{unred}}
\qquad\text{i.e.}\qquad
X_{\mathrm{unred}} = X_{\mathrm{red}} / \hat J.
\]

AMC with `reduce=true` prints **\(X_{\mathrm{red}}\)** (often an overall \(\hat J_0\) in the formula).
Treating that printout as unreduced (or unpacking with only 2 CGs as if unreduced) fakes FAIL — typical fake ratios \(\sim 1/\hat J\) (e.g. \(0.577\) for \(J=1\)).

---

## AMC is correct — mismatch is almost always packaging

**Do not blame AMC first.** When m ↔ AMC / code disagrees:

1. Name both sides (reduced vs unreduced) — see hard rule above.
2. Open the AMC **input** `reduce=` flags and the **printed** \(\hat J\) powers in the `.tex`.
3. If ratios cluster near \(\hat J^{\pm 1}\) or \(1/\sqrt{2J+1}\), **stop** — packaging bug, not a wrong 6j/9j.

AMC’s rule (scalar LHS):

| AMC input | Print is | Use against m |
|---|---|---|
| default (`reduce` off) | **unreduced** \(X_{\mathrm{unred}}\) | 2-CG unpack of \(X_{\mathrm{unred}}\), or \(\sum\mathrm{CG}\,\mathrm{CG}\,X(m)\) |
| `reduce=true` | **reduced** \(X_{\mathrm{red}}\) | \(X_{\mathrm{red}}=\sum\mathrm{CG}\,\mathrm{CG}\,X(m)/\hat J\) |

Same bare sum \(S\) in T×T→S intermediates often appears with \(\hat J^{-2}\) (unreduced) vs \(\hat J^{-1}\) (reduced). That is intentional — not two different physics formulas.

**Fix code / compare to AMC as printed. Do not retune AMC equations or invent extra hats to match a buggy store.**

Locked counterexamples:

- \(\chi^\theta\) / \(\Gamma^{\mathrm{III}_c}\) (2026-07-29): code stored bare \(S\) as “reduced”; AMC reduced is \(S/\hat J\). After fixing the store, m ≡ AMC DIRECT ≡ Path B with **no** sign/hat fudges. See [factored_GIIIc/NOTES.md](factored_GIIIc/NOTES.md).
- \(\chi^\lambda\) / \(\Gamma^{\mathrm{IV}_c}\) (2026-07-29): fold \(Z(m)\) projects as \(S=\sum CG\,CG\,Z(m)\), \(Z_{\mathrm{red}}=S/\hat J\), \(Z_{\mathrm{unred}}=S/\hat J^2\) ≡ AMC `G4c_from_chi{,_reduced}`. Comparing bare \(S\) to AMC unreduced fakes \(\hat J^{\pm2}\). See [factored_GIV/NOTES.md](factored_GIV/NOTES.md).
- \(\chi^\kappa\) / \(\Gamma^{\mathrm{IV}_a}\) (2026-08-06): Path B `inv_direct(bar_CHI_VI_II)` stores \(-\chi_{\mathrm{unred}}\). Compare \(\chi_{\mathrm{red}}\equiv-\hat J\cdot\mathrm{inv}\). Gold equation is analyze \(\Omega_{aicd}\Gamma_{jcba}\) — **not** old arxiv \(\Omega_{ajbl}\Gamma_{ibka}\). See [factored_GIV/NOTES.md](factored_GIV/NOTES.md).

---

## Correct m ↔ reduced J compare

m-scheme gold (always physical):

\[
X(m)=\sum\ldots
\quad\text{(include \([\Omega^\lambda\times\Omega^\lambda]^{(0)}\) CG when \(T\times T\to S\))}
\]

Project to **reduced** scalar (matches AMC `reduce=true`):

\[
X^{J}_{\mathrm{red}}
=
\frac{1}{\hat J}
\sum_{m}
\mathrm{CG}_{\mathrm{bra}}\,\mathrm{CG}_{\mathrm{ket}}\,X(m).
\]

Then compare \(X^{J}_{\mathrm{red}}\) to the AMC/J formula **as printed with `reduce=true`**.

Wrong (do not do this):

- Unpack AMC reduced ME with **unreduced** 2-CG formula and compare to \(X(m)\).
- Compare \(X(m)\) directly to a reduced TBME.
- Mix Path A (`reduce=true`) with Path B / code that stores unreduced without a \(\hat J\) conversion.
- Pauli-skip \(m_i=m_j\) when projecting a **non-AS** AMC object \(\chi_{ijkl}\) (that changes \(X_{\mathrm{red}}\)).

---

## Checklist (every bench)

1. Open the AMC input: is `reduce=true` on the LHS? On intermediates?
2. Open the `.tex`: is there an overall \(\hat J_0\) (or missing one)?
3. State the compare: `m → X_red` vs AMC, or `m → X_unred` vs AMC.
4. For tensors: `IsReduced()` is normal; unpack uses **WE** (3 CGs / \(\hat J_1\)).
5. If ratios cluster near \(\hat J^{\pm 1}\) or \(1/\sqrt{2J+1}\), **stop** — packaging bug, not angular-momentum bug.
6. Ring gold pattern: `run/test_z_ring_mscheme_sign.py`. Cross: `run/test_omega_cross_mscheme.py`.

---

## Locked examples (do not re-litigate)

| Object | AMC | Truth vs m |
|---|---|---|
| Ring \(Z=\sum\Omega_{p\bar a rb}\Omega_{aqsb}\) | `reduce=true` → \(Z_{\mathrm{red}}\) | m → \(Z_{\mathrm{red}}\) ≡ Path A |
| Cross \(\chi=\sum\Omega_{ajkb}\Omega_{ibal}\) | `reduce=true` → \(\chi_{\mathrm{red}}\) | m → \(\chi_{\mathrm{red}}\) ≡ Path A |
| \(\chi^\theta\) (T×T→S ladder) | red \(S/\hat J\), unred \(S/\hat J^2\) | m ≡ AMC; bare-\(S\)-as-reduced was **false FAIL** |
| Earlier “AMC direct FAIL” on cross | — | **false FAIL** from unreduced unpack |

See also: [tts_ring.md](factored_fIIIa/tts_ring.md), [LESSONS.md](LESSONS.md).
