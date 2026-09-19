# GIVc Path B — Pandya parity channels & DGEMM matching

Goal: keep **Pandya → DGEMM → inv**. Odd-π T1 fail is a **channel / topology**
mismatch (gold vs ring), not a mid-factor scratch.

---

## 1. Ordinary vs cross-coupled parity

Tensor \(\Omega^\lambda\) / \(\chi^\lambda\) with parity \(\pi_\eta\):

| \(\pi_\eta\) | Ordinary MatEl \((J_{\mathrm{bra}},J_{\mathrm{ket}})\) | CC pair \((\mathrm{ch}_b,\mathrm{ch}_k)\) |
|---|---|---|
| even \(0\) | \(\pi_{\mathrm{bra}}=\pi_{\mathrm{ket}}\) | **same-π** only |
| odd \(1\) | \(\pi_{\mathrm{bra}}\neq\pi_{\mathrm{ket}}\) | **opp-π** only |

Selection already in ethS:

```text
(tb.parity + tk.parity) % 2 == parity_eta
```

He4 emax=1 count: even → all same-π CC pairs; odd → all opp-π. That filter is
correct for **both** gold and ring.

---

## 2. Two topologies → same CC slot, different natural feed

Gold (m / NOTES):

\[
K=\sum_{ab}\chi^\lambda_{ialb}\,\Omega_{bjak}.
\]

Ring (production adcb / tts_ring):

\[
K_{\mathrm{ring}}=\sum_{ab}\chi^\lambda_{ibal}\,\Omega_{ajkb}.
\]

Both Pandya into the **same** CC layout \(\bar\chi^{J_{il}J_{ab}}_{il,ab}\) /
\(\bar\Omega^{J_{ab}J_{kj}}_{ab,kj}\), and both use the same DGEMM skeleton

\[
\bar G^{J_{il}}
\;+\!=\;
w(J_{il},J_{ab},\lambda)\,
\bar\chi^{J_{il}J_{ab}}\,
\bar\Omega^{J_{ab}J_{il}},
\qquad
(\pi_{il}+\pi_{ab})\bmod 2=\pi_\eta.
\]

Parity of the **target** CC slot is identical:

\[
\pi_{il}+\pi_{ab}
=\pi_i+\pi_l+\pi_a+\pi_b
=\begin{cases}
\pi_{ia}+\pi_{lb} & \text{gold natural}\\
\pi_{ib}+\pi_{al} & \text{ring natural}
\end{cases}
=\pi_\eta.
\]

So opp-π DGEMM matching does **not** distinguish gold from ring. What
differs is which **ordinary** (bra,ket) parity pair feeds the bar.

### Natural feed (χ)

| Topology | Natural coupling | \(\pi\) on the two pairs (odd \(\pi_\eta\)) |
|---|---|---|
| gold \(\chi_{ialb}\) | \((ia)\|(lb)\) | \((\pi_{ia},\pi_{lb})\in\{(0,1),(1,0)\}\) |
| ring \(\chi_{ibal}\) | \((ib)\|(al)\) | for the same \((i,l,a,b)\): **half exact swap**, half same pair |

He4-like emax=1 index scan (odd \(\pi_\eta\)): 320 tuples with
\((\pi_{ia},\pi_{lb})=(\pi_{ib},\pi_{al})\), 320 with exact swap, 0 other.
Aggregate counts into a fixed CC slot \((0,1)\) look identical for gold and
ring — that is a red herring. Per ME, the two topologies still read
\(\chi_{ialb}\) vs \(\chi_{ibal}\) (different ordinary legs / MatEl).

Even \(\pi_\eta\): both naturals are same-π \((0,0)\) or \((1,1)\), T1 weight
forces \(n_l=n_b\), χ looks AS under \((a\leftrightarrow b)\), and
\(K\equiv K_{\mathrm{ring}}\) — the even-π lock that hid the issue.

Odd \(\pi_\eta\): \(n_l\neq n_b\) on nnz T1 → gold and ring put **different**
content into the **same** \((\mathrm{ch}_{il},\mathrm{ch}_{ab})\) matrix.
Filling that matrix with `PandyaChiAdcb` implements ring, not gold.
The DGEMM channel *list* does not need to change; the **Pandya fill** of
those channels does.

### Natural feed (Ω)

| Topology | Natural | After Pandya |
|---|---|---|
| gold \(\Omega_{bjak}\) | \((bj)\|(ak)\) | \(\bar\Omega_{ab,kj}\) |
| ring \(\Omega_{ajkb}\) | \((aj)\|(kb)\) | \(\bar\Omega_{ab,kj}\) via adcb |

Ω **is** fermionic AS, so gold and ring bars are related by a known phase /
index swap once the χ side is right. Prefer locking χ first; then map Ω via AS
identity from adcb rather than a second untested 9j.

---

## 3. Equations to implement (gold Path B, both π)

\[
\begin{aligned}
\bar\chi^{J_{il}J_{ab}}_{il,ab}
&=\mathrm{Pandya}_{(ia)(lb)\to(il)(ab)}
  \bigl[\chi^{J_{ia}J_{lb}}_{ialb}\bigr],
\\
\bar\Omega^{J_{ab}J_{kj}}_{ab,kj}
&=\mathrm{Pandya}_{(bj)(ak)\to(ab)(kj)}
  \bigl[\Omega^{J_{bj}J_{ak}}_{bjak}\bigr],
\\
\bar G^{J_{il}}
&=\sum_{J_{ab}}
  \frac{(-1)^{J_{il}}}{\hat J_{il}}\,
  \hat\lambda^{-1}(-1)^{J_{ab}+\lambda}\,
  \bar\chi^{J_{il}J_{ab}}\,
  \bar\Omega^{J_{ab}J_{il}},
\\
&\quad\text{CC filter: }
  (\pi_{il}+\pi_{ab})\bmod 2=\pi_\eta
  \text{ (already in code)},
\\
Z_{\mathrm{red}}
&=\hat J_0\sum_{J_p}\hat J_p\,
  \{j_l\,j_k\,J_0;\,j_j\,j_i\,J_p\}\,
  \bar G^{J_p}
  \quad\text{(no AMC-sample overall minus)}.
\end{aligned}
\]

**Mid \(w\)** above is the locked **ring** mid. It is valid for gold **only after**
the bars are the gold Pandya maps — do not retune \(w\) to compensate a wrong
`PandyaChiIalb`. If a locked odd-π bar ME still needs a different \(w\), derive
it from AMC/tts after the fwd Pandya lock, not by fitting leftover ratios.

---

## 4. What is *not* the bug

1. Missing opp-π CC pairs — filter already keeps them when \(\pi_\eta=1\).
2. Inv Pandya / √2 packaging — even π and odd T2 already pass.
3. Antisymmetrizing χ to force ring≡gold — occupancy-AS only; would change physics.
4. Dropping `PandyaChiIalb` into the adcb slots with the ring mid without an
   odd-first bar lock — already tried; corr\(\sim 0.3\)–\(0.5\).

---

## 5. Lock order (odd-first)

1. **One odd-π bar ME:** for fixed \((i,l,a,b,J_{il},J_{ab})\) with
   \(\pi_{il}\neq\pi_{ab}\), compare
   `PandyaChiIalb` (or corrected 9j) to an independent m/CG construction of
   \(\bar\chi_{il,ab}\) from \(\chi_{ialb}\). Even π is degenerate with adcb —
   do not use it as the lock.
2. **Same for \(\bar\Omega_{ab,kj}\) from \(\Omega_{bjak}\)** (or prove AS map
   from `PandyaOmegaAdcb`).
3. **Single mid product** \(\bar\chi\cdot\bar\Omega\) vs CG bare \(K\) at one
   odd leftover ME (FoldAS off).
4. Scale to 2n×2n DGEMM + inv; keep production Path B speed path.

---

## 6. Code pointers

- CC filter / DGEMM: `comm223_232_GIVc_pathB` in
  `src/FactorizedDoubleCommutator_eths.cc`
- Ring fill (legacy even-π only): `PandyaChiAdcb` / `PandyaOmegaAdcb`
- Gold fill (production, all π): `PandyaChiIalb` / `PandyaOmegaBjak`
  — AMC `G4c_gold_il_ab_pandya.txt` schemes `((1,-3),(2,-4))` /
  `((3,-1),(4,-2))`
- Even-π mid lock: `run/test_givc_midj_bare_vs_cg.py`
