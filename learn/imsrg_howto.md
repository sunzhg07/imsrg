# imsrg++ — AI-Friendly How-To Guide (§20–24: Implementation & Flow)

> **Split document.** This file covers implementing commutator terms, the BCH/Magnus flow, generators, and operator transforms (§20–24, ~8K tokens).  
> For data structures, conventions, and APIs, see [imsrg_reference.md](imsrg_reference.md) (§1–19, ~18K tokens).

## 20. Implementing New Commutator Terms — Code Generation Guide

This section is the primary reference for converting a mathematical expression (diagram, algebraic formula) into working C++ code. Read this before writing any new commutator term.

### 20.1 Occupation Number Convention

```cpp
Orbit& oa = modelspace->GetOrbit(a);
double na    = oa.occ;        // n_a: hole=1, particle=0 (fractional for ensemble refs)
double nbara = 1.0 - oa.occ;  // n̄_a = 1 - n_a: hole=0, particle=1
```

So in a formula $n_a \bar{n}_b$ (`na * nbara`): a must be a **hole** ($n_a=1$), b must be a **particle** ($n_b=0, \bar{n}_b=1$). A factor $n_a(1-n_b) + (1-n_a)n_b = n_a + n_b - 2n_a n_b = $ ... is a general ph factor.

Common combinations in IMSRG(2):
| Factor | Meaning |
|---|---|
| `na` | a must be hole |
| `1-na` = `nbara` | a must be particle |
| `na * nb` | both holes (hh) |
| `nbara * nbarb` = `(1-na)*(1-nb)` | both particles (pp) |
| `na - nb` | ph, antisymmetric |
| `na*nb*(1-nc) + (1-na)*(1-nb)*nc` | two holes + one particle OR two particles + one hole (221 topology) |
| `(1-na)*(1-nb) - na*nb` | pp minus hh (222 pp/hh topology) |

Skip negligible terms (important for performance):
```cpp
if (std::abs(na) < 1e-9) continue;      // OCC_CUT = 1e-9
if (std::abs(nbara) < 1e-9) continue;
```

### 20.2 Angular Momentum: SixJ Call Convention

**CRITICAL**: `AngMom::SixJ` and `AngMom::CG` take **actual j values** (doubles, e.g. 0.5, 1.0, 1.5), **NOT doubled integers**.

```cpp
// orbit.j2 = 2j (stored as int). Must divide by 2 before passing to AngMom functions!
double sixj = AngMom::SixJ(oi.j2 * 0.5,   // j_i (actual)
                           oj.j2 * 0.5,   // j_j
                           (double)J,      // J  (already integer for integers)
                           ok.j2 * 0.5,   // j_k
                           ol.j2 * 0.5,   // j_l
                           (double)Jp);   // J'

double cg = AngMom::CG(oi.j2*0.5, mi, oj.j2*0.5, mj, (double)J, M);
double phase = AngMom::phase(x);  // returns (-1)^x, x must be integer
```

For the Pandya transformation formula $\bar{X}^J_{\bar{i}j\bar{k}l} = -\sum_{J'}(2J'+1)\{^{j_i\;j_j\;J}_{j_k\;j_l\;J'}\} X^{J'}_{ijkl}$:
```cpp
double sixj = AngMom::SixJ(oi.j2*0.5, oj.j2*0.5, (double)J,
                           ok.j2*0.5, ol.j2*0.5, (double)Jp);
xbar_ijkl -= (2*Jp + 1) * sixj * X2.GetTBME_J(Jp, Jp, i, j, k, l);
```

Use the cache instead for the repeated J loop in 2-body commutators:
```cpp
double sixj = Z.modelspace->GetSixJ(j1, j2, j3, J1, J2, J3); // cached version
```

### 20.3 Matrix Element Convention: Tilde vs Normalized

**Always be explicit about which convention you use.**

```cpp
// Returns TILDE (unnormalized): sqrt(1+dij)*(1+dkl) * Gamma  [USE IN FORMULAS]
double xijkl = X.TwoBody.GetTBME_J(J, J, i, j, k, l);

// Returns NORMALIZED: Gamma directly  [USE FOR STORING, or when you know what you want]
double xijkl = X.TwoBody.GetTBME_J_norm(J, J, i, j, k, l);
```

**Rule**: The IMSRG flow equations are written in the **tilde** convention. When you see $\tilde{\Gamma}^J_{ijkl}$ in a paper, use `GetTBME_J`. When computing $Z^J_{ijkl}$ to store, your accumulated `zijkl` is in tilde convention — you must convert to normalized before storing:

```cpp
// After computing zijkl in tilde convention (from tilde inputs):
if (i == j) zijkl /= PhysConst::SQRT2;   // un-tilde the bra
if (k == l) zijkl /= PhysConst::SQRT2;   // un-tilde the ket
// Now store normalized value:
Z.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);
```

`PhysConst::SQRT2 = sqrt(2)` is from `PhysicalConstants.hh`.

### 20.4 One-Body Channel Restriction

One-body matrix elements $Z_{ij}$ are non-zero only if $i$ and $j$ have the same $(l, j, t_z)$. Use `GetOneBodyChannel` to enforce this:

```cpp
for (size_t i : Z.modelspace->all_orbits) {
    Orbit& oi = Z.modelspace->GetOrbit(i);
    // Only loop over j with same l, j2, tz2 as i:
    for (size_t j : Z.GetOneBodyChannel(oi.l, oi.j2, oi.tz2)) {
        // ... compute Z.OneBody(i, j) ...
    }
}
```

For a **Hermitian** result, only compute $i \leq j$ and fill the other half:
```cpp
for (size_t j : Z.GetOneBodyChannel(oi.l, oi.j2, oi.tz2)) {
    if (j < i) continue;   // skip below diagonal (Hermitian)
    // ... compute zij ...
    Z.OneBody(i, j) += zij;
    if (i != j)
        Z.OneBody(j, i) += hZ * zij;  // hZ = +1 hermitian, -1 anti-hermitian
}
```

### 20.5 Template: 2-body → 1-body with Full J Sum

**Formula** (e.g. `comm121ss`):
$$Z_{ij} = \sum_{ab} \sum_J \frac{2J+1}{2j_i+1} (n_a - n_b) X_{ia} \tilde{Y}^J_{biaj}  - (X \leftrightarrow Y)$$

```cpp
void my_comm_121(const Operator& X, const Operator& Y, Operator& Z)
{
    int hZ = Z.IsHermitian() ? 1 : -1;
    // build all_orbits as vector for indexing
    for (size_t i : Z.modelspace->all_orbits) {
        Orbit& oi = Z.modelspace->GetOrbit(i);
        for (size_t j : Z.GetOneBodyChannel(oi.l, oi.j2, oi.tz2)) {
            if (j < i) continue;      // Hermitian: upper triangle only
            Orbit& oj = Z.modelspace->GetOrbit(j);

            double zij = 0;
            for (size_t a : Z.modelspace->all_orbits) {
                Orbit& oa = Z.modelspace->GetOrbit(a);
                for (size_t b : Z.modelspace->all_orbits) {
                    Orbit& ob = Z.modelspace->GetOrbit(b);
                    double na = oa.occ, nb = ob.occ;
                    // J range from triangle inequalities (doubled integers -> divide by 2)
                    int Jmin = std::max(std::abs(ob.j2-oi.j2), std::abs(oa.j2-oj.j2)) / 2;
                    int Jmax = std::min(ob.j2+oi.j2, oa.j2+oj.j2) / 2;
                    for (int J = Jmin; J <= Jmax; J++) {
                        // GetTBME_J returns tilde ME
                        double ybiaj = Y.TwoBody.GetTBME_J(J, J, b, i, a, j);
                        double xbiaj = X.TwoBody.GetTBME_J(J, J, b, i, a, j);
                        zij += (2*J+1) / (oi.j2+1.0) * (na - nb) * (X.OneBody(i,a)*ybiaj - Y.OneBody(i,a)*xbiaj);
                    }
                }
            }
            Z.OneBody(i, j) += zij;
            if (i != j) Z.OneBody(j, i) += hZ * zij;
        }
    }
}
```

### 20.6 Template: 2-body → 2-body pp/hh (Reference-level, Formula-direct)

**Formula** (`comm222_pp_hhss`):
$$Z^J_{ijkl} = \frac{1}{2}\sum_{ab}\left[(1-n_a)(1-n_b) - n_a n_b\right]\left(\tilde{X}^J_{ijab}\tilde{Y}^J_{abkl} - \tilde{Y}^J_{ijab}\tilde{X}^J_{abkl}\right)$$

```cpp
void my_comm_222_pphhss(const Operator& X, const Operator& Y, Operator& Z)
{
    // Collect all {ch_bra, ch_ket} pairs that exist in Z
    std::vector<size_t> ch_bra_list, ch_ket_list;
    for (auto& iter : Z.TwoBody.MatEl) {
        ch_bra_list.push_back(iter.first[0]);
        ch_ket_list.push_back(iter.first[1]);
    }
    int nch = ch_bra_list.size();

    #pragma omp parallel for schedule(dynamic, 1)
    for (int ich = 0; ich < nch; ich++) {
        size_t ch_bra = ch_bra_list[ich];
        size_t ch_ket = ch_ket_list[ich];
        TwoBodyChannel& tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
        TwoBodyChannel& tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
        int J = tbc_bra.J;

        int nbras = tbc_bra.GetNumberKets();
        int nkets = tbc_ket.GetNumberKets();

        for (int ibra = 0; ibra < nbras; ibra++) {
            Ket& bra = tbc_bra.GetKet(ibra);
            size_t i = bra.p, j = bra.q;

            // For Hermitian Z and diagonal channel: only upper triangle
            size_t iket_min = (ch_bra == ch_ket) ? ibra : 0;
            for (int iket = iket_min; iket < nkets; iket++) {
                Ket& ket = tbc_ket.GetKet(iket);
                size_t k = ket.p, l = ket.q;

                double zijkl = 0;
                for (size_t a : Z.modelspace->all_orbits) {
                    for (size_t b : Z.modelspace->all_orbits) {
                        if (b < a) continue;   // avoid double-counting; factor of 2 below
                        Orbit& oa = Z.modelspace->GetOrbit(a);
                        Orbit& ob = Z.modelspace->GetOrbit(b);
                        double na = oa.occ, nb = ob.occ;
                        double flip = (a == b) ? 1.0 : 2.0;   // b<a contribute equally
                        double fac = flip * 0.5 * ((1-na)*(1-nb) - na*nb);
                        if (std::abs(fac) < 1e-9) continue;
                        double xijab = X.TwoBody.GetTBME_J(J, J, i, j, a, b);
                        double yijab = Y.TwoBody.GetTBME_J(J, J, i, j, a, b);
                        double xabkl = X.TwoBody.GetTBME_J(J, J, a, b, k, l);
                        double yabkl = Y.TwoBody.GetTBME_J(J, J, a, b, k, l);
                        zijkl += fac * (xijab * yabkl - yijab * xabkl);
                    }
                }
                // Convert tilde z to normalized before storing
                if (i == j) zijkl /= PhysConst::SQRT2;
                if (k == l) zijkl /= PhysConst::SQRT2;
                Z.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);
            }
        }
    }
}
```

### 20.7 Template: 2-body → 2-body pp/hh (Production: Matrix Multiply)

Replace the O(N^6) reference loop with O(N^4) matrix multiplies using precomputed `KetIndex_pp`/`KetIndex_hh` subsets as `arma::uvec`:

```cpp
// Get per-channel subsets (arma::uvec for fast sub-matrix extraction):
auto& kets_pp    = tbc.GetKetIndex_pp();    // indices of pp kets
auto& kets_hh    = tbc.GetKetIndex_hh();    // indices of hh kets
auto& nanb_hh    = tbc.Ket_occ_hh;         // n_a * n_b for each hh ket (arma::vec)
auto& nbarnbar_hh= tbc.Ket_unocc_hh;       // (1-na)*(1-nb) for each hh ket

// Reference full matrices in the channel:
arma::mat& Xmat = X.TwoBody.GetMatrix(ch, ch);
arma::mat& Ymat = Y.TwoBody.GetMatrix(ch, ch);
arma::mat& Zmat = Z.TwoBody.GetMatrix(ch, ch);

// pp contribution: sum_abpp (1-na)(1-nb) Xijab Yabkl
// = X.cols(pp) * Y.rows(pp)  [pp kets have na=nb=0, so factor=1]
if (kets_pp.size() > 0)
    Zmat += X.cols(kets_pp) * Ymat.rows(kets_pp);

// hh contribution with sign and occ factors:
if (kets_hh.size() > 0) {
    Zmat += Xmat.cols(kets_hh) * arma::diagmat(nbarnbar_hh) * Ymat.rows(kets_hh);  // (1-na)(1-nb) part
    Zmat -= Xmat.cols(kets_hh) * arma::diagmat(nanb_hh)     * Ymat.rows(kets_hh);  // -na*nb part
}

// Symmetrize for Hermitian Z in diagonal channel:
if (Z.IsHermitian() && ch_bra == ch_ket)
    Zmat += Zmat.t();
// Then subtract Y*X (for [X,Y]):
// ... (handled in ConstructScalarMpp_Mhh with an XY - YX structure)
```

**Note**: `Xmat.cols(kets_pp)` extracts the columns of X indexed by `kets_pp` — this is a single BLAS call and the dominant operation in the pp/hh ladder. The matrix then has shape `[N_kets x N_pp]`.

### 20.8 Template: 2-body → 2-body ph (Pandya, Production)

The ph term requires the Pandya recoupling. In production code:

```cpp
void my_comm_222_ph(const Operator& X, const Operator& Y, Operator& Z)
{
    int hY = Y.IsHermitian() ? 1 : -1;
    size_t nch_cc = Z.modelspace->GetNumberTwoBodyChannels_CC();

    // Allocate intermediate Zbar in cross-coupled (CC) channels
    std::deque<arma::mat> Zbar(nch_cc);
    for (size_t ch = 0; ch < nch_cc; ch++) {
        size_t n = Z.modelspace->GetTwoBodyChannel_CC(ch).GetNumberKets();
        Zbar[ch].zeros(n, 2*n);   // 2 columns: direct and exchange terms
    }

    // Step 1: Pandya-transform X and Y into CC channels, contract to get Zbar
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t ch = 0; ch < nch_cc; ch++) {
        arma::mat Ybar_ph, Xt_bar_ph;
        // Returns X in CC chan (with na(1-nb) factors) and Y in CC chan:
        DoPandyaTransformation_SingleChannel_XandY(X, Y, Xt_bar_ph, Ybar_ph, ch);
        if (Ybar_ph.size() < 1 or Xt_bar_ph.size() < 1) continue;
        // Matrix multiply in CC space: Zbar = Xt * Ybar
        Zbar[ch] += Xt_bar_ph * Ybar_ph;  // (schematically)
    }

    // Step 2: Inverse Pandya transform Zbar back to normal channels
    AddInversePandyaTransformation(Zbar, Z);
}
```

For the **reference** (slow, maximally clear) ph implementation, the Pandya formula is:
```cpp
// Pandya transformation of X:
// X̄^{J'}_{ibjα} = -Σ_J (2J+1) {ji jj J; jk jl J'} X^J_{ijkl}
double xbar = 0;
int Jmin = std::max(std::abs(oi.j2-ob.j2), std::abs(oa.j2-ol.j2)) / 2;
int Jmax = std::min(oi.j2+ob.j2, oa.j2+ol.j2) / 2;
for (int Jpp = Jmin; Jpp <= Jmax; Jpp++) {
    double sixj = AngMom::SixJ(oi.j2*0.5, ol.j2*0.5, (double)Jp,
                               oa.j2*0.5, ob.j2*0.5, (double)Jpp);
    xbar -= (2*Jpp+1) * sixj * X.TwoBody.GetTBME_J(Jpp, Jpp, i, b, a, l);
}
```

### 20.9 Template: 2-body → 1-body monopole (comm221 via ConstructScalarMpp_Mhh)

The production `comm221ss` first builds $M_{pp}$ and $M_{hh}$ (the pp and hh intermediate matrices), then contracts them with a 6j sum:

```cpp
TwoBodyME Mpp(Z.modelspace, Z.GetJRank(), Z.GetTRank(), Z.GetParity());
TwoBodyME Mhh(Z.modelspace, Z.GetJRank(), Z.GetTRank(), Z.GetParity());
ConstructScalarMpp_Mhh(X, Y, Z, Mpp, Mhh);  // fills both intermediates

// Then contract to 1-body:
for (size_t i : Z.modelspace->all_orbits) {
    Orbit& oi = Z.modelspace->GetOrbit(i);
    for (size_t j : Z.GetOneBodyChannel(oi.l, oi.j2, oi.tz2)) {
        if (j < i) continue;
        double zij = 0;
        for (auto& c : Z.modelspace->all_orbits) {
            Orbit& oc = Z.modelspace->GetOrbit(c);
            double nc    = oc.occ;
            double nbarc = 1.0 - nc;
            int Jmin = std::max(std::abs(oc.j2-oi.j2), std::abs(oc.j2-oj.j2)) / 2;
            int Jmax = (oc.j2 + std::min(oi.j2, oj.j2)) / 2;
            for (int J = Jmin; J <= Jmax; J++) {
                if (std::abs(nc)    > 1e-9)
                    zij += (2*J+1) * nc    * Mpp.GetTBME_J(J, J, c, i, c, j);
                if (std::abs(nbarc) > 1e-9)
                    zij += (2*J+1) * nbarc * Mhh.GetTBME_J(J, J, c, i, c, j);
            }
        }
        Z.OneBody(i, j) += zij / (oi.j2 + 1.0);
        if (i != j) Z.OneBody(j, i) += hZ * zij / (oi.j2 + 1.0);
    }
}
```

### 20.10 Function Boilerplate

Every commutator function in this codebase follows this pattern:

```cpp
void comm_XYZss(const Operator& X, const Operator& Y, Operator& Z)
{
    double t_start = omp_get_wtime();    // timing start

    // ... body ...

    X.profiler.timer[__func__] += omp_get_wtime() - t_start;  // record time
}
```

And it must be:
1. Declared in `Commutator.hh` (inside `namespace Commutator`)
2. Added to `comm_term_on` map in `Commutator.cc` (default `true` for IMSRG(2), `false` for IMSRG(3))
3. Called from `CommutatorScalarScalar()` (or the appropriate dispatcher) with a `if (comm_term_on["..."])` guard

### 20.11 Hermitian/Anti-Hermitian Convention in Output

The generator `Eta` is **anti-Hermitian** (`IsAntiHermitian()`), $H$ is **Hermitian** (`IsHermitian()`). The commutator $[H, \eta]$ is Hermitian, $[\eta, H]$ is anti-Hermitian.

```cpp
int hZ = Z.IsHermitian()     ? +1 : -1;   // +1 for Hermitian, -1 for anti-Hermitian
// When filling upper triangle: Z(j,i) = hZ * Z(i,j)
Z.OneBody(j, i) = hZ * Z.OneBody(i, j);

// For 2-body: AddToTBME handles the hermitian update automatically at Level 1/2.
// At Level 3 (local indices): must do it manually if using AddToTBMENonHerm,
// or just use AddToTBME which does it for you.
```

### 20.12 J Range from Triangle Inequalities (Doubled Integer Arithmetic)

Given orbits with doubled-integer j2 values, the allowed two-body J range is:

```cpp
// For ket |ab J>: triangle inequality on a, b
int Jmin_ab = std::abs(oa.j2 - ob.j2) / 2;   // integer division (j2 are even for integer j)
int Jmax_ab = (oa.j2 + ob.j2) / 2;

// For a term involving <ij J| ... |kl J> with intermediate <cd J'>:
// where cd must couple to J (e.g. for loop over J'):
int J_prime_min = std::max(std::abs(oc.j2 - oi.j2), std::abs(od.j2 - oj.j2)) / 2;
int J_prime_max = std::min(oc.j2 + oi.j2, od.j2 + oj.j2) / 2;
for (int Jp = J_prime_min; Jp <= J_prime_max; Jp++) { ... }

// Also check parity and Tz conservation:
if ((oi.l + oj.l) % 2 != (ok.l + ol.l) % 2) continue;  // parity
if (oi.tz2 + oj.tz2 != ok.tz2 + ol.tz2) continue;       // isospin
```

### 20.13 Exchange Phase When Swapping Indices

When swapping two indices in a two-body ket:

```cpp
// |pq J> = (-1)^{jp + jq - J + 1} * |qp J>
int exchange_phase = Z.modelspace->phase((oi.j2 + oj.j2)/2 - tbc.J + 1);
// or equivalently using ket.phase_prefactor * (-1)^J (stored precomputed):
int ph = bra.Phase(J);  // = bra.phase_prefactor * (-1)^J
```

`Z.modelspace->phase(x)` returns $(-1)^x$.

### 20.14 Complete Example: New comm122-type Term

To implement a new 1-body × 2-body → 2-body term from a formula such as:
$$Z^J_{ijkl} = \sum_a A_{ia} \tilde{Y}^J_{ajkl} - (X \leftrightarrow Y)$$

```cpp
#include "Commutator.hh"
#include "AngMom.hh"
#include "PhysicalConstants.hh"

namespace Commutator {

void comm_my_122ss(const Operator& X, const Operator& Y, Operator& Z)
{
    double t_start = omp_get_wtime();

    // Collect channel pairs
    std::vector<size_t> ch_bra_list, ch_ket_list;
    for (auto& iter : Z.TwoBody.MatEl) {
        ch_bra_list.push_back(iter.first[0]);
        ch_ket_list.push_back(iter.first[1]);
    }
    int nch = ch_bra_list.size();

    // Use SortedTwoBodyChannels would be preferred, but MatEl loop above is general
    #pragma omp parallel for schedule(dynamic, 1)
    for (int ich = 0; ich < nch; ich++) {
        size_t ch_bra = ch_bra_list[ich];
        size_t ch_ket = ch_ket_list[ich];
        TwoBodyChannel& tbc_bra = Z.modelspace->GetTwoBodyChannel(ch_bra);
        TwoBodyChannel& tbc_ket = Z.modelspace->GetTwoBodyChannel(ch_ket);
        int J = tbc_bra.J;

        for (int ibra = 0; ibra < (int)tbc_bra.GetNumberKets(); ibra++) {
            Ket& bra = tbc_bra.GetKet(ibra);
            size_t i = bra.p, j = bra.q;

            size_t iket_start = (ch_bra == ch_ket) ? ibra : 0;
            for (int iket = iket_start; iket < (int)tbc_ket.GetNumberKets(); iket++) {
                Ket& ket = tbc_ket.GetKet(iket);
                size_t k = ket.p, l = ket.q;
                double zijkl = 0;

                for (size_t a : Z.modelspace->all_orbits) {
                    // GetTBME_J(J_bra, J_ket, i, j, k, l) — for scalar: J_bra==J_ket
                    double yajkl = Y.TwoBody.GetTBME_J(J, J, a, j, k, l);
                    double xiakl = X.TwoBody.GetTBME_J(J, J, i, a, k, l);
                    // ... apply formula ...
                    zijkl += X.OneBody(i, a) * yajkl - Y.OneBody(i, a) * xiakl;
                }

                // Normalize: divide out tilde normalization of bra and ket
                if (i == j) zijkl /= PhysConst::SQRT2;
                if (k == l) zijkl /= PhysConst::SQRT2;

                // AddToTBME with local indices: stores normalized value,
                // and handles the hermitian conjugate update automatically
                Z.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);
            }
        }
    }
    X.profiler.timer[__func__] += omp_get_wtime() - t_start;
}

} // namespace Commutator
```

### 20.15 Critical API Details: GetTBME_J, GetMatrix, and Loop Strategies

#### `GetTBME_J` — Two Overloads, Very Different Meanings

```cpp
// 1-argument form (SCALAR only): J_bra == J_ket, infers channel from orbital quantum numbers
double v = X.TwoBody.GetTBME_J(J, a, b, c, d);

// 2-argument form (GENERAL): explicit bra and ket J values — REQUIRED for tensor operators
// Also universally used in reference implementations even for scalars (J_bra == J_ket):
double v = X.TwoBody.GetTBME_J(J_bra, J_ket, a, b, c, d);
```

**Always use the 2-argument form** to avoid silent bugs when adapting code between scalar and tensor contexts. In all reference implementations you will see `X2.GetTBME_J(J, J, b, i, a, j)` — the repeated `J` is deliberate.

Norm variants:
```cpp
GetTBME_J(J_bra, J_ket, a,b,c,d)       // returns TILDE (multiplied by sqrt((1+dab)(1+dcd)))
GetTBME_J_norm(J_bra, J_ket, a,b,c,d)  // returns NORMALIZED (stored value)
GetTBME_J_norm_twoOps(OtherTBME, J_bra, J_ket, a,b,c,d, out_this, out_other)
// ↑ fetch from two operators simultaneously, one recoupling cost
```

#### `GetMatrix` — Zero-Copy Reference to the Underlying `arma::mat`

```cpp
// Returns reference: NO copy, NO allocation — required for large matrix operations
arma::mat& Xmat = X.TwoBody.GetMatrix(ch_bra, ch_ket);  // off-diagonal or diagonal block
arma::mat& Xmat = X.TwoBody.GetMatrix(ch);              // shorthand for (ch, ch)

// Const version (on const Operator):
const arma::mat& Xmat = X.TwoBody.GetMatrix(ch_bra, ch_ket);
```

**Critical**: `GetMatrix(ch_bra, ch_ket)` will throw `std::out_of_range` if the block doesn't exist (i.e., symmetry selection rules forbid it). Always guard with the channel pairing logic from `Z.TwoBody.MatEl` first.

For missing half (ch_bra > ch_ket, only ch_bra ≤ ch_ket is stored), use Hermiticity:
```cpp
// For Hermitian X, ch_bra > ch_ket block:
auto& X_block = (ch_bra <= ch_ket)
    ? X.TwoBody.GetMatrix(ch_bra, ch_ket)
    : X.TwoBody.GetMatrix(ch_ket, ch_bra).t() * hX;  // hX = +1 Hermitian, -1 AntiHerm
```

#### When to Loop Over `SortedTwoBodyChannels` vs `MatEl`

| Situation | Loop pattern | Reason |
|---|---|---|
| Scalar-scalar, 2b→2b | `for (ich in SortedTwoBodyChannels)` → `ch = SortedTwoBodyChannels[ich]` | Only diagonal {ch,ch} blocks exist; sorted for OpenMP load balance |
| Tensor or parity-changing, 2b→2b | `for (auto& iter : Z.TwoBody.MatEl)` | Off-diagonal {ch_bra,ch_ket} blocks may exist |
| 1b→2b or 2b→1b (any) | `for (ich in SortedTwoBodyChannels)` | Only J-diagonal 2-body blocks needed |
| Collecting ch_bra/ch_ket lists first | `for (auto& iter : Z.TwoBody.MatEl) { ch_bra_list.push_back(...) }` | Then parallelize over the list with OpenMP |

The `MatEl` key `{ch_bra, ch_ket}` has `ch_bra <= ch_ket` always.

---

### 20.16 Antisymmetry, Hermiticity, and CC Channel Symmetries

#### Operator Antisymmetry (Pauli Exclusion)

All two-body operators in this code act on **antisymmetrized** kets. The exchange relation is:

$$|abJ\rangle = (-1)^{j_a + j_b - J + 1} |baJ\rangle$$

In code:
```cpp
// From comm122ss: to find ket |aj> given |ja>:
int flipphaseij = -Z.modelspace->phase((oi.j2 + oj.j2)/2 - tbc.J);
// = (-1)^{(j2_i + j2_j)/2 - J + 1}  = (-1)^{j_i + j_j - J + 1}

// When looking up ket (a,j) where a > j, canonical order is (j,a):
size_t ind = tbc.GetLocalIndex(std::min(a,j), std::max(a,j));
double factor = (a > j) ? flipphaseij : 1.0;
// factor = (-1)^{j_a+j_j-J+1} because |aj J> = flipphaseij * |ja J>
```

`GetTBME_J` and `GetTBME` handle all of this automatically. You only need to apply the phase manually when doing index lookups at Level 3 (local ket index).

**Identical-particle factor (`SQRT2`)**: When `i == j` in a ket `|ii J>`, the normalized ME differs by $\sqrt{2}$:
```cpp
// The tilde norm:   tilde = sqrt(1 + d_{ab}) * normalized
// If i == j:        tilde = sqrt(2) * normalized
// So when accumulating zijkl in tilde convention and i==j or k==l:
if (i == j) zijkl /= PhysConst::SQRT2;  // before AddToTBME
if (k == l) zijkl /= PhysConst::SQRT2;
```

#### Hermitian/Anti-Hermitian Reduction (Computation Savings)

For a Hermitian output $Z$ and a **diagonal channel** ($ch\_bra == ch\_ket$):
- Only compute `ibra <= iket` (upper triangle)
- The lower triangle is $Z_{ijkl}(iket,ibra) = +Z_{ijkl}(ibra,iket)$ (Hermitian) or $-$ (anti-Hermitian)
- `AddToTBME` at Level 3 handles this automatically for diagonal channels

For **off-diagonal channels** ($ch\_bra \neq ch\_ket$):
- Only `ch_bra < ch_ket` block is computed and stored
- The `ch_bra > ch_ket` block is obtained on demand via `.t()` × `hX`
- Compute all `(ibra, iket)` pairs — no triangle restriction

```cpp
// In reference implementations, the pattern is explicit:
size_t iket_min = (ch_bra == ch_ket) ? ibra : 0;   // triangle for diagonal, full for off-diagonal
for (int iket = iket_min; iket < nkets; iket++) {
    // ... compute zijkl ...
    Z.TwoBody.AddToTBME(ch_bra, ch_ket, ibra, iket, zijkl);
    // AddToTBME automatically updates (iket, ibra) for diagonal channels
}
```

For 1-body, exploit Hermiticity explicitly:
```cpp
Z.OneBody(i, j) += zij;
if (i != j) Z.OneBody(j, i) += hZ * zij;   // hZ = +1 Hermitian, -1 AntiHermitian
```

#### Cross-Coupled (CC) Channel — Pandya ph Symmetry

The CC channel differs from the normal channel in three ways:

| Property | Normal `TwoBodyChannel` | Cross-Coupled `TwoBodyChannel_CC` |
|---|---|---|
| Ket ordering | `p <= q` (canonical: smaller index first) | No ordering restriction; both orderings of same quantum numbers appear |
| Tz conservation | `tz_p + tz_q = 2 Tz` (sum conserved) | `|tz_p - tz_q| = 2 Tz` (difference conserved) |
| Pauli exclusion | `p != q` if `J` odd | No Pauli restriction |
| Ket subsets | pp, hh, ph, cc, vv, ... | hh and ph (only these appear in ph diagrams) |

A CC ket is labeled `|i j_bar J_cc>` where the bar denotes the Pandya-transformed partner. The `KetIndex_hh` and `KetIndex_ph` of a CC channel refer to the **first** orbit being a hole or particle respectively.

**CC channel symmetry used in `comm222_phss`**: The Pandya-transformed $Z$-bar satisfies:

$$\bar{Z}^{J}_{i\bar{l}, k\bar{j}} = (-1)^{j_i + j_l + j_k + j_j} \cdot h_Z \cdot \bar{Z}^{J}_{l\bar{i}, j\bar{k}}$$

This means only computing the "direct" half of $\bar{Z}$ and recovering the "exchange" half by:
```cpp
// PhaseMatZ: (-1)^{(j_p + j_q)/2} per CC ket (stored as matrix for vectorized apply)
arma::mat PhaseMatZ(nKets_cc, nKets_cc, arma::fill::ones);
for (index_t iket = 0; iket < nKets_cc; iket++) {
    const Ket& ket = tbc_cc.GetKet(iket);
    if (Z.modelspace->phase((ket.op->j2 + ket.oq->j2)/2) < 0) {
        PhaseMatZ.col(iket) *= -1;
        PhaseMatZ.row(iket) *= -1;
    }
}
// Exchange = transpose * phase * hZ:
Zbar_ch.tail_cols(nKets_cc) += Zbar_ch.tail_cols(nKets_cc).t() % PhaseMatZ;
```

**Hermitian/anti-Hermitian simplification in CC space**: Because $[X,Y]$ with $X$ Hermitian and $Y$ anti-Hermitian gives a Hermitian $Z$ (or vice versa for $[\eta, H]$):

```cpp
// If Z is Hermitian (X hermitian != Y hermitian):
Zbar_ch.head_cols(nKets_cc) += Zbar_ch.head_cols(nKets_cc).t();   // ADD transpose
// If Z is anti-Hermitian (both same hermiticity):
Zbar_ch.head_cols(nKets_cc) -= Zbar_ch.head_cols(nKets_cc).t();   // SUBTRACT transpose
```

This halves the matrix multiply cost in CC space — only the "direct" contraction $\bar{X} \cdot \bar{Y}$ is computed, and the "exchange" $\bar{Y} \cdot \bar{X}$ is obtained from the transpose relation.

---

### 20.17 Checklist for a New Commutator Term

1. **Formula**: write out the full algebraic expression with correct tilde convention
2. **Topology**: identify N×M→P body structure — determines which template to follow
3. **Occupation factors**: write out $n_a$, $\bar{n}_a$ factors and verify ph topology
4. **J coupling**: identify which J runs are internal sums vs. fixed by channel
5. **SixJ**: if 6j appears, use `AngMom::SixJ` with actual j values (`j2/2` for each orbit)
6. **Tilde convention**: fetch with `GetTBME_J` (tilde), then normalize before `AddToTBME`
7. **Hermitian symmetry**: only compute upper triangle for diagonal channels; fill lower with `hZ * ...`
8. **Profiler**: add `X.profiler.timer[__func__] += ...` at start/end
9. **Register**: add to `comm_term_on` map in `Commutator.cc`
10. **Call site**: add call in `CommutatorScalarScalar()` gated by `comm_term_on["..."]` flag
11. **Test**: compare against corresponding function in `ReferenceImplementations.cc` (slow but guaranteed correct)


---

## 21. IMSRG Flow Equation and the Magnus Integrator

### 21.1 Physics: The Flow Equation

The IMSRG flow equation is:

$$\frac{dH(s)}{ds} = [\eta(s), H(s)]$$

where $s$ is the flow parameter and $\eta(s)$ is the **generator** (anti-Hermitian: $\eta^\dagger = -\eta$). The flow drives $H$ toward block-diagonal form. The **Magnus expansion** parametrizes the solution as $H(s) = e^{\Omega(s)} H_0 e^{-\Omega(s)}$ where $\Omega(s)$ is the **Magnus operator** (also anti-Hermitian). In practice the code takes an Euler Magnus step $\Omega(s+ds) \approx \Omega(s) + \eta(s)\,ds$, then recovers $H(s)$ via BCH.

### 21.2 The Solve Loop (`Solve_magnus_euler` — Default Method)

```
for each step:
  1. norm_eta = Eta.Norm()
  2. if norm_eta < eta_criterion (1e-6): converged, break
  3. if ||Omega.back()|| > omega_norm_max (2.0): split (NewOmega or GatherOmega)
  4. Adaptive:  ds = min(norm_domega/||eta||,
                          norm_domega/||eta||/(||omega||+eps),
                          omega_norm_max/||eta||,
                          ds_max)
  5. Euler step:     Eta *= ds
  6. BCH product:    Omega.back() = BCH_Product(Eta, Omega.back())
                     [exp(Omega_new) = exp(ds*eta) * exp(Omega_old)]
  7. BCH transform:  H(s) = BCH_Transform(H_saved, Omega.back())
  8. Update:         Eta = eta_func(H(s))
```

Key defaults (`IMSRGSolver` constructor):
- `ds = 0.1`, `ds_max = 0.5`, `smax = 2.0`
- `norm_domega = 0.1` (target per-step `||dOmega||`)
- `omega_norm_max = 2.0` (triggers Omega splitting)
- `eta_criterion = 1e-6` (convergence)
- `magnus_adaptive = true`

Other methods: `"magnus_backoff"` (adaptive with back-off on eta growth), `"magnus_modified_euler"` (predictor-corrector), `"flow_RK4"` (direct 4th-order RK on `dH/ds`), `"flow_adaptive"` / `"magnus_adaptive"` (boost::odeint adaptive steppers).

### 21.3 BCH_Transform: $e^\Omega X e^{-\Omega}$

Implements the nested commutator series (truncated at IMSRG(2) or IMSRG(3)):

```cpp
// Standard_BCH_Transform (BCH.cc):
double epsilon = nx * exp(-2*ny) * bch_transform_threshold / (2*ny);
Operator OpNested = OpOut;
for (int i = 1; i <= 40; ++i) {
    OpNested = Commutator::Commutator(Omega, OpNested);  // ith nested commutator
    factorial_denom /= i;
    OpOut += factorial_denom * OpNested;
    epsilon *= (i + 1);
    if (OpNested.Norm() < epsilon) break;
}
```

`bch_transform_threshold = 1e-9`. Warns at 12 iterations, gives up at 40.

### 21.4 BCH_Product: $\exp(Z) = \exp(X)\exp(Y)$

BCH series using Bernoulli numbers (up to 9th order):
```cpp
// Z = X + Y + (1/2)[Y,X] + (1/12)([Y,[Y,X]] - [X,[X,Y]]) - ...
std::vector<double> bernoulli = {1, -0.5, 1./6, 0, -1./30, 0, 1./42, 0, -1./30};
Z = X + Y;
Nested = Commutator(Y, X);
// special: also adds (1/12)*[Nested, X] if ||Nested||*||X|| > threshold
while (nxy > bch_product_threshold) {   // threshold = 1e-4
    Z += (bernoulli[k] / factorial[k]) * Nested;
    k++;
    if (k >= 9 or 2*ny*nxy < threshold) break;
    Nested = Commutator(Y, Nested);
}
```

Used at every Magnus step: `Omega.back() = BCH_Product(ds*eta, Omega.back())`.

### 21.5 Omega Splitting and H_saved

When `||Omega.back()|| > omega_norm_max`, the Magnus operator is split to keep BCH inexpensive:

**Standard `NewOmega()`:**
```
H_saved = FlowingOps[0]       // snapshot of H(s) now
Omega.push_back(zero)          // start a fresh Omega
// future: H(s) = BCH_Transform(H_saved, Omega.back())
```

**BCH_Transform call logic:**
```cpp
if (Omega.size() + n_omega_written < 2)
    H(s) = BCH_Transform(*H_0, Omega.back());    // very first Omega segment
else
    H(s) = BCH_Transform(H_saved, Omega.back()); // nth segment
```

Omega segments can be flushed to binary disk files via `FlushOmegaToScratch()`. `Transform_Partial` reads them back in sequence. `n_omega_written` tracks how many are on disk.

### 21.6 Hunter-Gatherer Mode

Enabled by `SetHunterGatherer(true)`. Always keeps exactly **two** Omega entries:
- `Omega.back()` = **hunter**: small, updated each step
- `Omega[size-2]` = **gatherer**: cumulative total $\Omega$

When `||hunter|| > omega_norm_max`, `GatherOmega()`:
```cpp
gatherer = BCH_Product(hunter, gatherer);   // exp(g_new) = exp(hunter)*exp(gatherer)
hunter.Erase();
H_saved = *H_0;
for (i = 0; i < Omega.size()-1; i++)        // all non-hunter Omegas
    H_saved = BCH_Transform(H_saved, Omega[i]);
```

**Advantage**: hunter stays small → fewer BCH nested commutators per step. At convergence a single gatherer holds the complete $\Omega_{\text{total}}$.

After `Solve()`, `GatherOmega()` is called once more to flush any remaining hunter into the gatherer.

---

## 22. Generators: Building $\eta(s)$

### 22.1 Generator Architecture

`Generator::Update(H_s, Eta_s)` calls `AddToEta`, which dispatches by `generator_type` string:

```cpp
void Generator::AddToEta(Operator& H_s, Operator& Eta_s) {
    if      (generator_type == "white")            ConstructGenerator_SingleRef(white_func);
    else if (generator_type == "atan")             ConstructGenerator_SingleRef(atan_func);
    else if (generator_type == "shell-model")      ConstructGenerator_ShellModel(white_func);
    else if (generator_type == "shell-model-atan") ConstructGenerator_ShellModel(atan_func);
    // ...
}
```

All generators:
1. Identify **off-diagonal** matrix elements of $H$
2. Compute **energy denominator** $\Delta$
3. Apply: $\eta_{ia} = f(H_{ia}, \Delta)$
4. Enforce anti-Hermitian: `Eta(a,i) = -Eta(i,a)`

### 22.2 Functional Forms

All are `std::function<double(double,double)>` of `(Hod, denom)`:

| `generator_type` | Formula | Notes |
|---|---|---|
| `"white"` | $H_{od}/\Delta$ | Default. Fast; can diverge near degeneracies |
| `"atan"` | $\frac{1}{2}\arctan(2H_{od}/\Delta)$ | Bounded $\pm\pi/4$; robust |
| `"wegner"` | $H_{od}\cdot\Delta$ | Original Wegner; slow convergence |
| `"imaginary-time"` | $H_{od}\cdot\operatorname{sgn}(\Delta)$ | Gradient descent-like |
| `"qtransfer-atan_N"` | $|\Delta M_N/\hbar^2|^{N/2}\cdot\text{atan}(\ldots)$ | Soft high-$q$ regulator |
| `"shell-model"` | white on shell-model blocks | Decouples core AND valence |
| `"shell-model-atan"` | atan on shell-model blocks | Robust shell-model version |
| `"shell-model-atan-npnh"` | atan + NpNh core excitations | For multi-reference truncations |
| `"hartree-fock"` | 1B only: $H_{ij}/\Delta_{ij}$ all pairs | Rarely used |

### 22.3 Energy Denominators

**1-body** ($\eta_{ia}$, $i$ = particle, $a$ = hole):
```
Moller-Plesset:  Δ = f_ii - f_aa
Epstein-Nesbet:  Δ += (n_i - n_a) * Γ_monopole(i,a,i,a)
MP_isospin:      Δ = (f_ii - f_aa + f_{ī,ī} - f_{ā,ā}) / 2   [ī = isospin partner orbit]
```

**2-body** ($\eta_{ij,ab}$, $ij$ = particle pair, $ab$ = hole pair):
```
Moller-Plesset:  Δ = f_ii + f_jj - f_aa - f_bb
Epstein-Nesbet:  Δ += (1-n_i-n_j)*Γ_mono(ij,ij) - (1-n_a-n_b)*Γ_mono(ab,ab)
                    + (n_i-n_a)*Γ_mono(ia,ia) + (n_i-n_b)*Γ_mono(ib,ib)
                    + (n_j-n_a)*Γ_mono(ja,ja) + (n_j-n_b)*Γ_mono(jb,jb)
```

**3-body** (Möller-Plesset only): $\Delta = f_{ii}+f_{jj}+f_{kk} - f_{aa}-f_{bb}-f_{cc}$

All denominators clipped: `if (|Δ| < denominator_cutoff) Δ = denominator_cutoff` (default 1e-6).

`GetTBMEmonopole(a,b,c,d)` = angle-averaged $\bar\Gamma_{abcd} = \frac{\sum_J(2J+1)\Gamma^J_{abcd}}{\sum_J(2J+1)}$.

### 22.4 Single-Reference Generator (white/atan/imaginary-time)

Decouples core from everything: zeros $\langle\text{particle}\ldots|H|\text{hole}\ldots\rangle$.

**1-body**: `i ∈ (valence ∪ qspace)`, `a ∈ core`:
```cpp
Eta.OneBody(i,a) = etafunc(H.OneBody(i,a), Get1bDenominator(i,a));
Eta.OneBody(a,i) = -Eta.OneBody(i,a);
```

**2-body**: `ibra ∈ (qq ∪ vv ∪ qv)` pairs, `iket ∈ cc` pairs:
```cpp
for (iket in tbc.GetKetIndex_cc())
  for (ibra in GetKetIndex_qq() ∪ GetKetIndex_vv() ∪ GetKetIndex_qv())
    ETA2(ibra, iket) = etafunc(H2(ibra, iket), denominator);
    ETA2(iket, ibra) = -ETA2(ibra, iket);   // enforces anti-symmetry
```

Ket subsets: `cc` = core-core (holes), `qq` = qspace-qspace (particles), `vv` = valence-valence, `qv`= qspace-valence.

### 22.5 Shell-Model Generator

Extends to decouple both core and the valence space (for VS-IMSRG):
- **1-body**: zero `(i,a)` for `a ∈ (core ∪ valence)`, `i ∈ (valence ∪ qspace)`
- **2-body, core**: η for `(vv ∪ qv ∪ qq)` vs `(cc ∪ vc)` kets
- **2-body, valence**: η for `(qv ∪ qq)` vs `vv` kets

`"shell-model-atan-npnh"` additionally includes NpNh (core-core off-diagonal) elements.

### 22.6 3-Body Generator

When `H.GetParticleRank() > 2 and not only_2b_eta`:
- Off-diagonal: bra = all v or q (particle) orbits, ket = all c (core/hole) orbits
- Denominator: Möller-Plesset $\Delta = \sum_\text{bra} f_{ii} - \sum_\text{ket} f_{aa}$
- dE3max and `occnat` cuts applied identically to 3B commutator loops
- Parallelized with `#pragma omp parallel for schedule(dynamic,1)` over 3B channels

---

## 23. Transforming Operators After the Flow

### 23.1 Three-Step Protocol

```cpp
// Step 1: build bare operator in the same modelspace
Operator O_bare = AllowedGamowTeller_Op(modelspace);

// Step 2: normal-order
Operator O_NO = O_bare.DoNormalOrdering();

// Step 3: apply stored unitary transformation e^Omega ... e^{-Omega}
Operator O_transformed = imsrgsolver.Transform(O_NO);
// Internally: sequential BCH_Transform(BCH_Transform(..., Omega[0])..., Omega[k])
```

`Transform(Op)` = `Transform_Partial(Op, 0)`: walks all stored Omega segments in order (from disk, then memory), applying `BCH_Transform` at each step.

`InverseTransform(Op)`: applies `−Omega_k, ..., −Omega_0` in reverse order.

### 23.2 Available Transition Operators

| String key | Factory | $J^\pi$ | Description |
|---|---|---|---|
| `"Fermi"` | `AllowedFermi_Op` | $0^+$, T=1 | $\hat\tau_\pm$ allowed Fermi |
| `"GamowTeller"` | `AllowedGamowTeller_Op` | $1^+$, T=1 | $\hat\sigma\hat\tau_\pm$ GT |
| `"E1"` – `"E3"` | `ElectricMultipoleOp(ms,L)` | $L^{(-1)^L}$ | Electric multipole B(EL) |
| `"M1"`, `"M2"` | `MagneticMultipoleOp(ms,L)` | $L^{(-1)^{L+1}}$ | Magnetic multipole B(ML) |
| `"ISQ"` / `"IVQ"` | `MultipoleResponseOp(ms,2,2,0/1)` | $2^+$ | IS/IV quadrupole response |
| `"ISM"` / `"IVM"` | `MultipoleResponseOp(ms,2,0,0/1)` | $0^+$ | IS/IV monopole response |
| `"Rp2"`, `"Rn2"`, `"Rm2"` | `Rp2_corrected_Op`, etc. | $0^+$ | Point-proton/neutron/matter radii |
| `"DGT"` | `M0nu::DGT_Op` | $0^+$ | Double GT operator |
| `"0vbb"` / `"0vbb_GT"` | `M0nu::GamowTeller(ms, Eclosure, src, hGT_AA)` | $0^+$ | $0\nu\beta\beta$ GT |
| `"0vbb_F"` | `M0nu::Fermi(ms, Eclosure, src, hF_VV)` | $0^+$ | $0\nu\beta\beta$ Fermi |

All are `Operator` instances with `rank_J`, `rank_T`, `parity` set to match the physics.

### 23.3 Tensor Operators in BCH_Transform

Tensor operators (`rank_J > 0`, e.g. E2, GT) use the same `Standard_BCH_Transform`:
- Auto-upgrade: `SetParticleRank(2)` if `nlegs%2==0 and rank<2`
- Auto-`MakeReduced()` if `rank_J==0 and (rank_T≠0 or parity≠0) and not IsReduced()`
- `Commutator::Commutator(Omega, OpNested)` dispatches to tensor commutators (`comm122st`, `comm222_pp_hhst`, etc. in `Commutator.cc`)
- Per-step diagnostics: prints 1B, 2B, 3B, total norms of `OpNested`

### 23.4 0νββ Operators (`M0nu.cc`)

Two-body scalar ($J=0$, $T=0$) operators violating lepton number. Key variants:

| Function | Channel | Term |
|---|---|---|
| `GamowTeller(ms, Eclosure, src, hGT_AA)` | GT | Axial-Axial (dominant) |
| `GamowTeller(ms, Eclosure, src, hGT_AP)` | GT | Axial-Pseudoscalar |
| `GamowTeller(ms, Eclosure, src, hGT_PP)` | GT | Pseudoscalar-Pseudoscalar |
| `GamowTeller(ms, Eclosure, src, hGT_MM)` | GT | Magnetic-Magnetic |
| `Fermi(ms, Eclosure, src, hF_VV)` | Fermi | Vector-Vector |
| `GamowTellerHeavy(ms, src, hGT_AA)` | GT | Heavy Majorana neutrino |

`Eclosure` = closure neutrino propagator energy. `src` = short-range correlation function (`"miller"`, `"AV18"`, `"CD-Bonn"`).

### 23.5 Co-Flowing Multiple Operators (Flow Method)

For ODE methods (`"flow_RK4"`, `"flow_adaptive"`):
```cpp
imsrgsolver.AddOperator(GT_op);        // FlowingOps[1]
imsrgsolver.AddOperator(E2_op);        // FlowingOps[2]
imsrgsolver.Solve();
// All operators get dO_i/ds = [eta(s), O_i(s)] at each RK sub-step
Operator GT_final = imsrgsolver.GetOperator(1);
```

For the **Magnus method** (`"magnus_euler"`): `FlowingOps` holds only $H(s)$. Observables are post-transformed by calling `solver.Transform(Op)` after `Solve()` finishes.

### 23.6 Parity-Violating (PV) Operators

For PC+PV mixed Hamiltonians, two operators flow with cross-coupled commutators:
$$\frac{d\mathcal{O}_{PV}}{ds} = [\eta_{PV}, \mathcal{O}_{PC}] + [\eta_{PC}, \mathcal{O}_{PV}]$$

```cpp
auto [O_pc, O_pv] = BCH::BCH_TransformPV(Op, Op_PV, Omega, Omega_PV);
```

`BCH_TransformHPV`: $H$ evolves under PC-only Omega; $H_{PV}$ gets both cross-contributions.

---

## 24. Goose-Tank and Factorized Corrections

### 24.1 The IMSRG(2) Truncation Error

At each BCH step, `[Omega_2, X_2]` generates 3- and 4-body pieces discarded at IMSRG(2). The leading missed contribution is a 4th-order quadruple where the discarded 3B term feeds back into the next nested commutator. Two corrections recover parts of this.

### 24.2 Goose-Tank Correction (`use_goose_tank_correction`)

Enable: `BCH::SetUseGooseTank(true)`.

Tracks an **auxiliary 1-body operator** $\chi$ representing the 2B→1B projected missing contribution:

```cpp
// Standard_BCH_Transform, at iteration i:
auto chi_last = chi.OneBody;          // chi from previous iteration
chi = GooseTankUpdate(Omega, OpNested);
OpNested.OneBody += chi_last;         // feed previous chi into current nested commutator
```

`GooseTankUpdate`:
```cpp
comm221ss(Omega, OpNested, chi);          // 2B+1B → 1B commutator
chi(i,j) *= n_i*n_j + (1-n_i)*(1-n_j);  // keep only pp and hh blocks
```

### 24.3 Factorized Double Commutator (`use_factorized_correction`)

Enable: `BCH::SetUseFactorizedCorrection(true)`.

At BCH step $i \geq 2$, augments the current nested commutator with trident-diagram contributions:
```cpp
comm223_231(Omega, chi2, Op_correction);  // 2+2+3→2+3 factorized
comm223_232(Omega, chi2, Op_correction);
OpNested += Op_correction;
```

`chi2` = nested commutator from **two steps back**. Recovers leading 3B corrections without explicitly storing 3-body operators. Also applied in `BCH_Product` (controlled by `use_factorized_correction_BCH_Product`).

Flags: `TurnOn_comm223_231`, `TurnOn_comm223_232` (both `true` by default).

### 24.4 Brueckner BCH Transform (`use_brueckner_bch`)

Enable: `BCH::SetUseBruecknerBCH(true)`.

Splits $\Omega = \Omega_1 + \Omega_2$ (1B + 2B parts) and applies two sequential BCH transforms:
```cpp
OpOut = BCH_Transform(OpIn, Omega1);   // 1-body rotation first
OpOut = BCH_Transform(OpOut, Omega2);  // then 2-body rotation
```

Can improve convergence for certain Hamiltonians (T.D. Morris variant).

---

## 25. Tensor-Tensor Commutators — Current Gap and Implementation Roadmap

### 25.1 What Exists

Three top-level dispatchers are declared in `Commutator.hh` (lines 59–61):
```cpp
Operator CommutatorScalarScalar(const Operator &X, const Operator &Y);
Operator CommutatorScalarTensor(const Operator &X, const Operator &Y);
Operator CommutatorScalarDagger(const Operator &X, const Operator &Y);
```

`CommutatorTensorTensor` is **NOT declared anywhere**. The dispatcher at `Commutator.cc` line ~278 explicitly hard-exits:
```cpp
else  // Both are tensor. Uh oh.
{
    std::cout << " Tensor-Tensor commutator not yet implemented." << std::endl;
    std::exit(EXIT_FAILURE);
}
```

The `[T,S]` case is handled by symmetry — the dispatcher calls `-CommutatorScalarTensor(Y, X)`.

### 25.2 Full Scalar-Tensor (`st`) Function Inventory

All `st` functions live in `TensorCommutators.cc` (production) and `ReferenceImplementations.cc` (reference). They are called through `CommutatorScalarTensor`.

**IMSRG(2) — on by default:**

| `comm_term_on` key | Function | Status |
|---|---|---|
| `"comm111st"` | `comm111st(X,Y,Z)` | production |
| `"comm121st"` | `comm121st(X,Y,Z)` | production |
| `"comm221st"` | `comm221st(X,Y,Z)` | production |
| `"comm122st"` | `comm122st(X,Y,Z)` | production |
| `"comm222_pp_hhst"` | `comm222_pp_hh_221st(X,Y,Z)` (combined with 221) | production |
| `"comm222_phst"` | `comm222_phst(X,Y,Z)` | production |

**IMSRG(3) n7 — off by default (`SetUseIMSRG3N7_Tensor`):**

`comm331st`, `comm231st`, `comm132st`, `comm232st`, `comm133st`, `comm223st` — all pass unit tests.

**IMSRG(3) full — off by default (`SetUseIMSRG3_Tensor`):**

`comm332_ppph_hhhpst`, `comm332_pphhst`, `comm233_pp_hhst`, `comm233_phst`, `comm333_ppp_hhhst`, `comm333_pph_hhpst` — all pass unit tests.

### 25.3 Why Tensor-Tensor is Harder

For $[X^{J_X}, Y^{J_Y}]^{J_Z}$ with both $J_X, J_Y > 0$:

- The output rank $J_Z$ ranges over $|J_X - J_Y| \leq J_Z \leq J_X + J_Y$ — a single commutator produces a **family** of output operators, not one.
- Every diagram picks up an extra recoupling coefficient (9j symbol) relative to the `st` case where $J_X=0$ reduces the 9j to a 6j.
- The `CommutatorTensorTensor` dispatcher must accept $J_Z$ as a parameter and construct `Z` with the appropriate rank.
- All loop indices over two-body channels must sum over $J_Z$ consistently.

### 25.4 What Would Be Needed to Implement `comm222tt`

Given a LaTeX formula for `comm222tt` (or any `tt` term), the implementation requires:

1. **New dispatcher** in `Commutator.hh`/`Commutator.cc`:
```cpp
// In Commutator.hh:
Operator CommutatorTensorTensor(const Operator &X, const Operator &Y, int J_Z);

// In Commutator.cc — route from Commutator():
else  // Both are tensor
{
    // Loop over valid J_Z values and accumulate, OR
    // take J_Z as a parameter if the caller knows which rank is wanted.
    return CommutatorTensorTensor(X, Y, /*J_Z*/);
}
```

2. **New `comm_term_on` keys**: `"comm222_pp_hhtt"`, `"comm222_phtt"`, etc.

3. **New toggle functions**: `SetUseIMSRG3_TensorTensor(bool)` analogous to `SetUseIMSRG3_Tensor`.

4. **The term functions** themselves (e.g. `comm222_pp_hhtt(X, Y, Z)`) in a new `TensorTensorCommutators.cc`, following the pattern of `TensorCommutators.cc` but with 9j symbols in place of 6j where the `ss`/`st` codes have 6j/nothing.

5. **Unit tests** in `UnitTest.cc` following the pattern of the `st` test block at line ~1761.

### 25.5 Locating the `st` Source as a Template

When implementing a `tt` term, the corresponding `st` function is the closest template:

| Target `tt` function | Template `st` function | Key change |
|---|---|---|
| `comm222_pp_hhtt` | `comm222_pp_hh_221st` | replace 6j with 9j; add $J_Z$ loop |
| `comm222_phtt` | `comm222_phst` | Pandya step gains an extra 9j factor |
| `comm111tt` | `comm111st` | 1-body recoupling gains a factor $\hat{J}_X \hat{J}_Y / \hat{J}_Z$ |

The `SQRT2` constant from `PhysicalConstants.hh` and `AngMom::SixJ`, `AngMom::NineJ`, `AngMom::phase`, `AngMom::Hat` are the angular-momentum tools used throughout the `st` implementations.
