# imsrg++ — AI-Friendly Codebase Reference (§1–19: Data Structures & API)

> **Split document.** This file covers data structures, conventions, and APIs (§1–19, ~18K tokens).  
> For implementing new terms, flow equations, and operator transforms, see [imsrg_howto.md](imsrg_howto.md) (§20–24, ~8K tokens).

**Language**: C++17  
**Author**: Ragnar Stroberg (original), extended by collaborators  
**Purpose**: J-coupled In-Medium Similarity Renormalization Group (IMSRG) for *ab initio* nuclear structure calculations  
**Linear algebra backend**: [Armadillo](http://arma.sourceforge.net/) (`arma::mat`, `arma::vec`, `arma::uvec`)  
**Python bindings**: via `pybind11` (`src/pyIMSRG.cc`)  
**Build system**: CMake (`CMakeLists.txt`, `build/`)  

---

## 1. Physics Background

### 1.1 What is IMSRG?

The IMSRG is a many-body method that unitarily transforms a nuclear Hamiltonian $H$ to a decoupled form using a continuous unitary transformation parameterized by a flow parameter $s$:

$$\frac{d H(s)}{ds} = [\eta(s), H(s)]$$

where $\eta(s)$ is the **anti-Hermitian generator** chosen to drive off-diagonal matrix elements to zero. The transformation is implemented via the **Magnus expansion**:

$$H(s) = e^{\Omega(s)} H(0) e^{-\Omega(s)}$$

where $\Omega(s)$ is an anti-Hermitian operator accumulated during the flow. The `Omega` equation is:

$$\frac{d\Omega}{ds} = \text{(series in nested commutators of } \eta \text{ and } \Omega)$$

### 1.2 Truncation Levels

- **IMSRG(2)**: Operators truncated at two-body level. The dominant approximation used in most calculations.
- **IMSRG(3)**: Operators truncated at three-body level (leading corrections, expensive).
- **MR-IMSRG**: Multi-reference variant for open-shell nuclei using a correlated reference state.
- **VS-IMSRG**: Valence-space IMSRG — decouples a valence shell from the core and Q-space.
- **EOM-IMSRG**: Equation-of-motion extension to compute excited states.

### 1.3 Normal Ordering

All operators are normal-ordered with respect to a reference state (typically a Hartree-Fock Slater determinant). A $k$-body operator normal-ordered to 2-body gives:

**2-body operator:**

$$E_0 = \sum_h (2j_h+1)\, n_h\, f_{hh} + \frac{1}{2}\sum_{hh'J}(2J+1)\,n_h n_{h'}\,\tilde{\Gamma}^J_{hh'hh'}$$

$$f_{ab} \mathrel{+}= \sum_h n_h \,\langle ah|\Gamma|bh \rangle$$

**3-body operator normal-ordered to 2-body:**

$$\Gamma^J_{ijkl} = V^J_{ijkl} + \sum_a n_a \sum_K \frac{2K+1}{2J+1} V^{(3),JJK}_{ijakla}$$

with prefactors $1/3$ for zero-body and $1/2$ for one-body when calling `DoNormalOrdering3`.

### 1.4 J-Coupling Convention

All two-body matrix elements are stored in the **J-coupled scheme**. The code uses doubled integers throughout for half-integer quantum numbers: `j2` = $2j$, `tz2` = $2t_z$.

The **tilde normalization** is used in flow equations:

$$\tilde{\Gamma}^J_{ijkl} = \sqrt{(1+\delta_{ij})(1+\delta_{kl})}\;\Gamma^J_{ijkl}$$

`GetTBME()` → returns $\tilde{\Gamma}$ (used in flow equations)  
`GetTBME_norm()` → returns $\Gamma$ (stored value)

### 1.5 Wigner-Eckart Convention

The code uses the **Edmonds/Brink-Satchler** convention (standard in nuclear physics):

$$\langle J\,M \,|\, \hat{O}^{(j)}_m \,|\, J'\,M' \rangle = (-1)^{2j}\, \frac{\langle J'M'\,jm\,|\,JM \rangle}{\sqrt{2J+1}}\, \langle J \| \hat{O}^{(j)} \| J' \rangle$$

- The $\sqrt{2J+1}$ is on the **bra side** (normalization of the bra state).
- **Scalar operators** (`rank_J=0`): stored as **non-reduced** matrix elements, i.e. the full $\langle JM|O|JM\rangle$.
- **Tensor operators** (`rank_J>0`): stored as **reduced** matrix elements $\langle J\|O^{(\lambda)}\|J'\rangle$.
- `MakeReduced()` → multiply by $\sqrt{2J+1}$ (one-body: $\sqrt{2j+1}$, two-body: $\sqrt{2J_{ch}+1}$)
- `MakeNotReduced()` → divide by $\sqrt{2J+1}$

---

## 2. Source File Map (`src/`)

| File | Role |
|---|---|
| `ModelSpace.hh/.cc` | Single-particle basis, channels, Kets, angular momentum caches |
| `Operator.hh/.cc` | Generic operator: 0+1+2+3-body, scalar or tensor |
| `TwoBodyME.hh/.cc` | Two-body matrix element storage and access |
| `ThreeBodyME.hh/.cc` | Three-body ME facade (delegates to storage backend) |
| `ThreeBodyStorage.hh` | Abstract base for 3BME storage backends |
| `ThreeBodyStorage_iso.hh/.cc` | Isospin-coupled 3BME storage |
| `ThreeBodyStorage_pn.hh/.cc` | Proton-neutron 3BME storage |
| `ThreeBodyStorage_no2b.hh/.cc` | NO2B projected 3BME storage |
| `ThreeBodyStorage_mono.hh/.cc` | Monopole-only 3BME storage |
| `ThreeLegME.hh/.cc` | 3-leg (dagger) operator matrix elements |
| `TwoBodyChannel.cc` | TwoBodyChannel initialization (declared in ModelSpace.hh) |
| `ThreeBodyChannel.cc` | ThreeBodyChannel initialization (declared in ModelSpace.hh) |
| `Commutator.hh/.cc` | IMSRG(2) commutator terms, Pandya transformation |
| `TensorCommutators.hh/.cc` | Commutator terms for tensor operators |
| `DaggerCommutators.hh/.cc` | Commutator terms for dagger (number-non-conserving) operators |
| `IMSRG3Commutators.hh/.cc` | IMSRG(3) commutator terms |
| `FactorizedDoubleCommutator.hh/.cc` | Factorized double-commutator corrections |
| `EOMFactorizedDoubleCommutator.hh/.cc` | EOM variant of factorized corrections |
| `BCH.hh/.cc` | Baker-Campbell-Hausdorff transformation |
| `Generator.hh/.cc` | Generator $\eta$ construction (Wegner, White, arctan, etc.) |
| `GeneratorPV.hh/.cc` | Parity-violating generator extension |
| `IMSRGSolver.hh/.cc` | ODE driver: Euler, Magnus, adaptive step |
| `IMSRGSolverPV.hh/.cc` | Parity-violating IMSRG solver |
| `HartreeFock.hh/.cc` | Hartree-Fock single-particle basis optimization |
| `HFMBPT.hh/.cc` | HF + many-body perturbation theory |
| `EOM.hh/.cc` | Equation-of-motion for excited states |
| `RPA.hh/.cc` | Random Phase Approximation |
| `AngMom.hh/.cc` | Angular momentum algebra (CG, 3j, 6j, 9j, 12j, Moshinsky) |
| `AngMomCache.hh/.cc` | Fast cached 6j lookup for 2-body commutators |
| `ReadWrite.hh/.cc` | I/O for operators, matrix elements (ME2J, OSLO, HDF5, etc.) |
| `imsrg_util.hh/.cc` | High-level utility functions for calculations |
| `IMSRGProfiler.hh/.cc` | Timing and counter profiler |
| `Jacobi3BME.hh/.cc` | Jacobi-to-lab-frame 3-body transformation |
| `M0nu.hh/.cc` | Neutrinoless double-beta decay matrix elements |
| `DarkMatterNREFT.hh/.cc` | Dark matter non-relativistic EFT operators |
| `Helicity.hh/.cc` | Helicity-basis operators |
| `ReferenceImplementations.hh/.cc` | Slow reference implementations for testing |
| `UnitTest.hh/.cc` | Unit tests |
| `Parameters.hh` | Input parameter parsing |
| `PhysicalConstants.hh` | Physical constants (MeV units) |
| `imsrg++.cc` | Main executable entry point |
| `pyIMSRG.cc` | Python bindings (pybind11) |

---

## 3. `ModelSpace` — The Foundational Data Structure

`ModelSpace` defines the entire single-particle and many-body basis. All other objects hold a pointer to a `ModelSpace`.

### 3.0 Key Typedef

```cpp
typedef unsigned long long int index_t;   // orbit index, ket index, channel index
```

`index_t` is used throughout the codebase for all integral indices into orbit, ket, and channel tables. Never use bare `int` for orbit/ket indices — use `index_t` or `size_t`.

### 3.1 Single-Particle Orbits — `struct Orbit`

```cpp
struct Orbit {
    int n;       // principal quantum number (0, 1, 2, ...)
    int l;       // orbital angular momentum
    int j2;      // 2*j  (doubled to keep integers)
    int tz2;     // 2*tz: proton = -1, neutron = +1
    double occ;      // occupation: hole=1, particle=0 (fractional for ensemble refs)
    double occ_nat;  // natural orbital occupation
    int cvq;     // classification: core=0, valence=1, qspace=2
    int index;   // unique global index
};
```

**Orbit string label format** (from `ModelSpace::Index2String`):
`{p|n}{n}{l-letter}{2j}` where l-letters are: s,p,d,f,g,h,i,j,k,l,m,n,o

| String | Meaning |
|---|---|
| `p0s1` | proton $0s_{1/2}$ |
| `p0d5` | proton $0d_{5/2}$ |
| `n1s1` | neutron $1s_{1/2}$ |
| `n0f7` | neutron $0f_{7/2}$ |
| `p0h11` | proton $0h_{11/2}$ |

**`cvq` classification** (important for VS-IMSRG):
- `cvq=0` → **core** orbit (below valence space; decoupled from valence)
- `cvq=1` → **valence** orbit (active space for shell-model diagonalization)
- `cvq=2` → **qspace** orbit (above valence; decoupled from valence)

The generator type `"shell-model"` drives off-diagonal elements between different `cvq` classes to zero.

### 3.2 Two-Body Kets — `struct Ket`

```cpp
struct Ket {   // |pq>  with p <= q (canonical ordering enforced)
    Orbit* op; Orbit* oq;
    size_t p, q;
    int dpq;           // 1 if p==q, else 0
    int phase_prefactor; // (-1)^{(jp+jq)/2 + 1}
    // Phase(J) = phase_prefactor * (-1)^J
    // Meaning: |pqJ> = Phase(J) * |qpJ>  (includes fermionic sign)
};
```

The exchange phase is: $|pqJ\rangle = (-1)^{j_p+j_q-J+1} |qpJ\rangle$

### 3.3 Three-Body Kets — `struct Ket3`

```cpp
struct Ket3 {  // |pqr; Jpq>  with canonical ordering handled externally
    size_t p, q, r;
    Orbit *op, *oq, *oR;
    int Jpq;   // intermediate coupling of orbits p and q
};
```

### 3.4 Two-Body Channels — `TwoBodyChannel`

Groups all two-body kets $|pqJ\rangle$ sharing the same conserved quantum numbers $(J, \pi, T_z)$:

```
TwoBodyChannel {
    int J;        // total angular momentum
    int parity;   // 0=even, 1=odd
    int Tz;       // total isospin projection (Tz = (tz_p + tz_q)/2)
    vector<int> KetList;   // local_idx -> global ket index
    vector<int> KetMap;    // global ket index -> local_idx  (-1 if not in channel)
    arma::uvec KetIndex_pp, KetIndex_hh, KetIndex_ph, ...  // pre-sorted subsets
    arma::vec  Ket_occ_hh, Ket_unocc_hh, Ket_occ_ph, ...  // precomputed occupations
}
```

Subsets: `pp`=particle-particle, `hh`=hole-hole, `ph`=particle-hole, `cc`=core-core, `vc`=valence-core, `qc`=qspace-core, `vv`=valence-valence, `qv`=qspace-valence, `qq`=qspace-qspace.

`TwoBodyChannel_CC`: the **cross-coupled** (Pandya-transformed) channels used for particle-hole diagrams. Index structure differs — $(J, \pi, T_z)$ for ph pairs.

### 3.5 Three-Body Channels — `ThreeBodyChannel`

Groups three-body kets $|pqr; J_{pq}\rangle$ by $(2J, \pi, 2T_z)$. Uses doubled integers for half-integer $J$.

### 3.6 Orbit Classification Sets

```cpp
set<index_t> holes;      // reference Slater determinant (n_i = 1)
set<index_t> particles;  // above reference (n_i = 0)
set<index_t> core;       // core for VS-IMSRG decoupling
set<index_t> valence;    // valence space orbits
set<index_t> qspace;     // above valence (excluded from VS)
set<index_t> proton_orbits, neutron_orbits, all_orbits;
```

### 3.7 Truncations

| Variable | Meaning |
|---|---|
| `Emax` | Single-particle: $e = 2n+l \leq E_\text{max}$ |
| `E2max` | Two-body: $e_i + e_j \leq E_{2\text{max}}$ |
| `E3max` | Three-body: $e_i + e_j + e_k \leq E_{3\text{max}}$ |
| `emax_3body_` | Can differ from `Emax` for 3-body space |
| `dE3max` | Relative energy cut for IMSRG(3) configurations |
| `occnat3cut` | Natural orbital occupation cut for 3-body states |
| `Lmax`, `Lmax2`, `Lmax3` | Angular momentum cuts for 1,2,3-body |

### 3.8 Precomputed Angular Momentum Quantities (Static, Global)

All cached as hash maps (`unordered_map<uint64_t, double>`):

| Cache | Content |
|---|---|
| `SixJList` | Wigner 6j symbols $\{j_1\,j_2\,j_3; J_1\,J_2\,J_3\}$ |
| `NineJList` | Wigner 9j symbols |
| `MoshList` | Moshinsky brackets for CM/relative transformation |
| `T3bList` | T-coefficients for Jacobi-to-lab 3-body transformation |
| `SixJCache_112112 six_j_cache_2b_` | Fast dedicated cache for 6j pattern dominating 2-body commutators |

`PandyaLookup`: precomputed channel index map for efficient Pandya transformations.

### 3.9 Predefined Valence Spaces (`ValenceSpaces`)

Static map of named shell-model spaces, e.g.:

| Name | Core | Valence orbits |
|---|---|---|
| `"s-shell"` | vacuum | $0s_{1/2}$ (p,n) |
| `"p-shell"` | He4 | $0p_{3/2}, 0p_{1/2}$ (p,n) |
| `"sd-shell"` | O16 | $0d_{5/2}, 0d_{3/2}, 1s_{1/2}$ (p,n) |
| `"fp-shell"` | Ca40 | $0f_{7/2}, 0f_{5/2}, 1p_{3/2}, 1p_{1/2}$ (p,n) |
| `"jj55"` | Sn100 | $0g_{7/2}, 1d_{5/2}, 1d_{3/2}, 2s_{1/2}, 0h_{11/2}$ (p,n) |
| ... | ... | ... |

### 3.10 `ModelSpace` Constructors

```cpp
// Minimal: HO basis up to emax, reference = Slater det of nucleus "He4","O16","Ca40",...
ModelSpace(int emax, std::string reference);

// With valence space: reference is closed-core nucleus, valence is named space
ModelSpace(int emax, std::string reference, std::string valence);
// Example: ModelSpace(4, "O16", "sd-shell")

// With separate 3-body emax (avoids memory blowup when emax_3body < emax):
ModelSpace(int emax, int emax_3body, std::string reference, std::string valence);

// With explicit hole list and valence list (as orbit strings):
ModelSpace(int emax, vector<string> hole_list, vector<string> valence_list);
// Example hole_list: {"p0s1", "p0p3", "p0p1", "n0s1", "n0p3", "n0p1"}

// With explicit core, hole, valence lists (MR-IMSRG style):
ModelSpace(int emax, vector<string> hole_list, vector<string> core_list, vector<string> valence_list);
```

Key setters called after construction:
```cpp
ms.SetHbarOmega(20.0);   // HO frequency in MeV (default 20)
ms.SetTargetMass(16);    // nuclear mass (for kinetic energy corrections)
ms.SetTargetZ(8);        // proton number
ms.SetE3max(14);         // if 3-body forces used
```

Look up orbit by quantum numbers:
```cpp
index_t i = ms.GetOrbitIndex(n, l, j2, tz2);     // returns orbit index
std::string s = ms.Index2String(i);               // returns "p0d5" etc.
Orbit& o = ms.GetOrbit(i);
int ch = ms.GetTwoBodyChannelIndex(J, parity, Tz);  // channel index from quantum numbers
TwoBodyChannel& tbc = ms.GetTwoBodyChannel(ch);
```

### 3.11 `SortedTwoBodyChannels` — OpenMP Load Balancing

```cpp
vector<unsigned int> SortedTwoBodyChannels;    // sorted by decreasing matrix dimension
vector<unsigned int> SortedTwoBodyChannels_CC; // same, for cross-coupled (Pandya) channels
```

All commutator OpenMP loops use these instead of `0..nch-1` to avoid load imbalance:

```cpp
int nch = modelspace->GetNumberTwoBodyChannels();
#pragma omp parallel for schedule(dynamic,1)
for ( int ich = 0; ich < nch; ich++ )
{
    int ch = modelspace->SortedTwoBodyChannels[ich];  // largest matrices first
    TwoBodyChannel& tbc = modelspace->GetTwoBodyChannel(ch);
    // ... heavy matrix multiply ...
}
```

Channels with many kets (large $J$, many orbits) get picked up first by threads, minimizing idle time.

`KetIndex_pp` etc. in `TwoBodyChannel` are `arma::uvec` so that Armadillo sub-column extraction `mat.cols(tbc.KetIndex_pp)` compiles to a single optimized BLAS call, enabling the pp/hh ladder diagram to be a pure matrix multiply:

```cpp
// pp-ladder: Z_pp += X.cols(pp) * diag(n_a*n_b) * Y.rows(pp)
// implemented as:
arma::mat X_pp = X_mat.cols(tbc.KetIndex_pp);
arma::mat Y_pp = Y_mat.rows(tbc.KetIndex_pp);   // .each_row scaled by occ
Z_mat += X_pp * Y_pp;
```

---

## 4. `Operator` — The Many-Body Operator Container

Every operator in the code — Hamiltonian, transition operator, EOM dagger, perturbation — is an instance of the same `Operator` class.

### 4.1 Body Structure

```cpp
class Operator {
    ModelSpace* modelspace;
    double      ZeroBody;       // E_0: scalar energy shift
    arma::mat   OneBody;        // f_{ij}: N_orb x N_orb matrix
    TwoBodyME   TwoBody;        // Gamma_{ijkl}^J: channel-blocked matrices
    ThreeBodyME ThreeBody;      // V_{ijklmn}^{J,T}: 3-body MEs
    ThreeLegME  ThreeLeg;       // (dagger operators only) Z_{ij,k}^J

    int rank_J;       // spherical tensor rank lambda
    int rank_T;       // isospin tensor rank
    int parity;       // 0=even, 1=odd
    int particle_rank;// 2 or 3
    int legs;         // = 2*particle_rank for normal ops, odd for daggers
    bool hermitian, antihermitian;
    bool is_reduced;  // true if storing reduced MEs
    index_t Q_space_orbit; // for dagger ops: the Q-space orbit index
}
```

### 4.2 Operator Classification

| Type | `rank_J` | `rank_T` | `parity` | `is_reduced` | `legs` | Example |
|---|---|---|---|---|---|---|
| Scalar Hamiltonian | 0 | 0 | 0 | false | 4 | $H$, $T$, $V_{NN}$ |
| Number operator | 0 | 0 | 0 | false | 4 | $\hat{N}$ |
| Gamow-Teller | 1 | 1 | 0 | true | 4 | $\hat{O}_{GT}$ |
| E2 transition | 2 | 0 | 0 | true | 4 | $\hat{Q}_{2m}$ |
| M1 transition | 1 | 0 | 1 | true | 4 | $\hat{\mu}$ |
| Dagger (EOM) | 0 | 0 | 0 | false | 3 | $A^\dagger_Q$ |
| 3-body Hamiltonian | 0 | 0 | 0 | false | 6 | $V_{3N}$ |

### 4.3 Arithmetic Overloads

All arithmetic acts body-by-body on `ZeroBody`, `OneBody`, `TwoBody`, `ThreeBody`:

```cpp
Op1 + Op2   // ZeroBody adds, OneBody adds, TwoBody.MatEl adds, etc.
Op * 3.0    // scales all bodies
3.0 * Op    // same (non-member)
Op1 - Op2
-Op
Op1 += Op2
Op *= 2.0
```

If one operator has higher `particle_rank`, the result is promoted automatically.

### 4.4 Normal Ordering Methods

```cpp
Operator DoNormalOrdering()              // auto-selects 2 or 3 body
Operator DoNormalOrdering2(sign, occupied_set)   // 2-body NO
Operator DoNormalOrdering3(sign, occupied_set)   // 3-body NO -> 2-body
Operator DoNormalOrderingDagger(sign, occupied)  // for dagger ops
Operator UndoNormalOrdering()            // sign=-1 version
Operator DoNormalOrderingCore()          // wrt core orbits
Operator UndoNormalOrderingCore()
Operator DoNormalOrderingFilledValence() // wrt filled valence
```

`sign = +1` → normal order, `sign = -1` → undo normal ordering.

### 4.5 One-Body Channel Structure

`OneBodyChannels[{l, 2j, 2tz}]` = set of orbit indices with those quantum numbers. Used to restrict one-body loops to symmetry-allowed pairs for tensor operators. Stored redundantly as `OneBodyChannels_vec` (flat vector) for O(1) access.

### 4.6 Perturbation Theory Utilities (built into Operator)

```cpp
double GetMP2_Energy()          // second-order MP energy correction
array<double,3> GetMP3_Energy() // third-order MP corrections
double GetMP2_3BEnergy()        // MP2 with 3-body forces
array<double,2> GetPPHH_Ladders() // pp/hh ladder diagrams
```

---

## 5. `TwoBodyME` — Two-Body Matrix Element Storage

### 5.1 Storage Layout

```cpp
map< array<size_t,2>, arma::mat > MatEl;
// Key: {ch_bra, ch_ket}
// Value: mat of size [N_kets(ch_bra) x N_kets(ch_ket)]
```

Only `ch_bra <= ch_ket` pairs are stored. The `ch_bra > ch_ket` half is reconstructed via hermiticity:

$$\Gamma^{J_1}_{ch_1, ch_2}(i,j) = (-1)^{J_1 - J_2} \Gamma^{J_2}_{ch_2, ch_1}(j,i) \quad \text{(hermitian)}$$

**Channel allocation rules** (set in `Allocate()`):

For a tensor operator with `rank_J=λ`, `rank_T=τ`, `parity=π`:
- Triangle$(J_{bra}, J_{ket}, \lambda)$ must be satisfied
- $|T_{z,bra} - T_{z,ket}| = \tau$
- $(p_{bra} + p_{ket} + \pi) \% 2 = 0$

For scalar ($\lambda=0, \tau=0, \pi=0$): only diagonal `{ch, ch}` blocks exist.

### 5.2 Access API — Three Levels

**Level 1: by orbit indices + channel index** (physics level)
```cpp
// Returns TILDE (unnormalized): sqrt(1+d_ab)(1+d_cd) * stored_value
double GetTBME     (ch_bra, ch_ket, a, b, c, d)
// Returns stored NORMALIZED value
double GetTBME_norm(ch_bra, ch_ket, a, b, c, d)
void   SetTBME     (ch_bra, ch_ket, a, b, c, d, value)
void   AddToTBME   (ch_bra, ch_ket, a, b, c, d, value)
```

**Level 2: by J quantum numbers** (convenience, infers channel)
```cpp
double GetTBME_J    (J_bra, J_ket, a, b, c, d)   // infers parity/Tz from orbits
double GetTBME_J_norm(J_bra, J_ket, a, b, c, d)
void   SetTBME_J   (J_bra, J_ket, a, b, c, d, v)
void   AddToTBME_J (J_bra, J_ket, a, b, c, d, v)
double GetTBME(J, parity, Tz, a, b, c, d)          // explicit quantum numbers
```

**Level 3: by local ket indices** (fastest, inner loops)
```cpp
double GetTBME_norm(ch_bra, ch_ket, ibra, iket)    // = mat(ibra, iket) directly
void   SetTBME     (ch_bra, ch_ket, ibra, iket, v)
void   AddToTBME   (ch_bra, ch_ket, ibra, iket, v)
void   AddToTBMENonHerm(...) // skips symmetry update — caller's responsibility
// GetME_pn_TwoOps: fetch same ME from two operators simultaneously (saves recoupling)
void GetTBME_J_norm_twoOps(OtherTBME, j_bra, j_ket, a,b,c,d, out1, out2)
```

**What Set/Add do automatically** (caller never needs to handle):
1. **Canonical ordering**: if $a > b$, applies exchange phase $(-1)^{j_a+j_b-J+1}$ and maps to $(\min, \max)$
2. **Channel swap**: if `ch_bra > ch_ket`, swaps and multiplies by $(-1)^{J_1-J_2}$
3. **Hermitian conjugate**: for `ch_bra == ch_ket`, simultaneously updates `mat(iket, ibra)` with $\pm$ sign

**Special accessors:**
```cpp
double GetTBMEmonopole(a, b, c, d)     // J-averaged: sum_J (2J+1)*V_J / sum_J (2J+1)
double Get_iso_TBME_from_pn(J, T, tz, a,b,c,d)  // pn -> isospin
void   Set_pn_TBME_from_iso(J, T, tz, a,b,c,d, v) // isospin -> pn
```

---

## 6. `ThreeBodyME` — Three-Body Matrix Element Storage

### 6.1 Coupling Scheme

Matrix elements stored in **unnormalized JT-coupled** form:

$$\langle (abJ_{ab}t_{ab})c;\, JT \,|\, V \,|\, (deJ_{de}t_{de})f;\, JT \rangle$$

with canonical ordering $a \geq b \geq c$ (bra) and $d \geq e \geq f$ (ket). Off-canonical permutations computed on the fly.

### 6.2 Storage Backend (Polymorphic)

`ThreeBodyME` is a facade holding a `unique_ptr<ThreeBodyStorage>`. Backend selected via `SetMode(string)`:

| Backend class | Mode string | Use case |
|---|---|---|
| `ThreeBodyStorage_iso` | `"isospin"` | Default on file read; isospin-coupled storage |
| `ThreeBodyStorage_pn` | `"pn"` | During IMSRG flow; proton-neutron basis |
| `ThreeBodyStorage_no2b` | `"no2b"` | NO2B projection; only $J_{ab}$-slice needed for NO |
| `ThreeBodyStorage_mono` | `"mono"` | Monopole-only; minimal memory |

`TransformToPN()`: converts isospin storage to pn (one-way; used before IMSRG(3) flow).

### 6.3 Access API

```cpp
// pn basis (most common in IMSRG)
ThreeBME_type GetME_pn(Jab, Jde, twoJ, a,b,c,d,e,f)
void SetME_pn        (Jab, Jde, twoJ, a,b,c,d,e,f, V)
void AddToME_pn      (Jab, Jde, twoJ, a,b,c,d,e,f, V)

// isospin basis
ThreeBME_type GetME_iso(Jab, Jde, twoJ, tab, tde, twoT, a,b,c,d,e,f)
void SetME_iso         (...)
void AddToME_iso       (...)

// fast channel+index access
ThreeBME_type GetME_pn_ch(ch_bra, ch_ket, ibra, iket)
void AddToME_pn_ch       (ch_bra, ch_ket, ibra, iket, V)
void SetME_pn_ch         (ch_bra, ch_ket, ibra, iket, V)

// fetch same ME from two operators simultaneously (saves recoupling in commutators)
vector<double> GetME_pn_TwoOps(Jab,Jde,twoJ,a,b,c,d,e,f, X, Y)

// special projections
ThreeBME_type GetME_pn_no2b(a,b,c,d,e,f, J2b)  // sum over J with (2J+1) weight
ThreeBME_type GetME_pn_mono(a,b,c,d,e,f)        // monopole
```

### 6.4 Permutation Handling (`ThreeBodyStorage` base)

When orbit order is non-canonical, `SortOrbits()` finds the permutation applied, then `RecouplingCoefficient()` computes:

$$\text{coeff} = (-1)^{\text{phase}} \times \{6j\}$$

The six permutations are enumerated as: `ABC` (identity), `BCA`, `CAB`, `ACB`, `BAC`, `CBA`.

`GetKetIndex_withRecoupling(Jab, twoJ, a,b,c, ibra_vec, coeff_vec)`: for a given $(a,b,c,J_{ab})$, returns all canonical kets that contribute and their recoupling weights.

---

## 7. `ThreeLegME` — Dagger Operator Storage

Used for **number-non-conserving** operators (odd `legs`), e.g. a dagger operator $A^\dagger_Q$ in EOM-IMSRG.

### 7.1 Storage Layout

```cpp
map<size_t, arma::mat> MatEl;
// Key: ch  (single two-body channel index)
// Value: mat of size [N_kets(ch) x N_all_orbits]
```

- **Bra side**: J-coupled two-body ket $|ijJ\rangle$ — same channel structure as `TwoBodyME`
- **Ket side**: flat orbit index $k$ over all orbits (no second channel)
- **Q orbit**: implicit, stored as `Q_space_orbit` on the parent `Operator`

### 7.2 Access

```cpp
double GetME     (ch, a, b, c)   // returns tilde: sqrt(1+d_ab) * stored
double GetME_norm(ch, a, b, c)   // stored value
double GetME_J   (J,  a, b, c)   // infers channel from J + orbital quantum numbers
void   SetME     (ch, a, b, c, v)
void   AddToME   (ch, a, b, c, v)
void   AddToME_J (J,  a, b, c, v)
```

`OneBody` for dagger: stored as `OneBody(i, 0)` — a single column. Q orbit encoded externally via `Q_space_orbit`, **not** as a column index.

---

## 8. `Commutator` Namespace — IMSRG(2) Flow Equations

The commutator $[X, Y] = Z$ is computed term by term. Naming convention: `commNMPss` where:
- `N` = particle rank of $X$ (in legs: 1-body=2, 2-body=4)
- `M` = particle rank of $Y$
- `P` = particle rank of $Z$
- `ss` = scalar-scalar

### 8.1 Scalar-Scalar Terms

| Function | Diagram | Description |
|---|---|---|
| `comm110ss` | $[f,f]\to E_0$ | 1b×1b → 0b |
| `comm220ss` | $[\Gamma,\Gamma]\to E_0$ | 2b×2b → 0b |
| `comm111ss` | $[f,f]\to f$ | 1b×1b → 1b |
| `comm121ss` | $[f,\Gamma]\to f$ | 1b×2b → 1b |
| `comm221ss` | $[\Gamma,\Gamma]\to f$ | 2b×2b → 1b |
| `comm122ss` | $[f,\Gamma]\to \Gamma$ | 1b×2b → 2b |
| `comm222_pp_hhss` | $[\Gamma,\Gamma]_{pp/hh}\to \Gamma$ | pp and hh ladders |
| `comm222_phss` | $[\Gamma,\Gamma]_{ph}\to \Gamma$ | ph (Pandya) term |
| `comm222_pp_hh_221ss` | combined pp/hh+221 | efficiency optimization |

### 8.2 Pandya Transformation

The ph term requires recoupling between pp and ph channels via the **Pandya transformation**:

$$\bar{X}^J_{\bar{i}j\bar{k}l} = -\sum_{J'} (2J'+1) \begin{Bmatrix} j_i & j_j & J \\ j_k & j_l & J' \end{Bmatrix} X^{J'}_{ijkl}$$

```cpp
void DoPandyaTransformation(Z, deque<mat>&, orientation)
void AddInversePandyaTransformation(deque<mat>&, Z)
// Single channel variant (parallelized):
void DoPandyaTransformation_SingleChannel(Z, mat&, ch_cc, orientation)
void DoPandyaTransformation_SingleChannel_XandY(X, Y, X_CC, Y_CC, ch_cc)
```

### 8.3 IMSRG(3) Terms

In `IMSRG3Commutators.cc`. Additional terms involving 3-body operators. Controlled by global flags:
```cpp
Commutator::SetUseIMSRG3(bool)      // include 3-body commutators
Commutator::SetUseIMSRG3N7(bool)    // include N^7 scaling terms
Commutator::SetUseIMSRG3_MP4(bool)  // MP4-level 3-body corrections
```
Individual terms can be toggled: `TurnOffTerm("comm330")`, `TurnOnTerm("comm330")`.

### 8.4 Top-Level `Commutator()` Dispatcher

The public entry point `Commutator(X, Y)` routes to the correct specialization based on `GetJRank()` and `IsNumberConserving()`:

```
X.IsNumberConserving() && Y.IsNumberConserving():
  X.JRank==0 && Y.JRank==0  →  CommutatorScalarScalar(X, Y)
  X.JRank==0 && Y.JRank!=0  →  CommutatorScalarTensor(X, Y)
  X.JRank!=0 && Y.JRank==0  →  -CommutatorScalarTensor(Y, X)   // [T,S]=-[S,T]
  X.JRank!=0 && Y.JRank!=0  →  std::exit(EXIT_FAILURE)         // NOT IMPLEMENTED
Y not number-conserving:  CommutatorScalarDagger(X, Y)  (or negative)
```

Before routing, if either operator `IsReduced()`, a non-reduced copy is made. The path for `both JRank!=0` prints a diagnostic and hard-exits — `CommutatorTensorTensor` does **not exist**.

### 8.5 Scalar-Tensor Commutators (`CommutatorScalarTensor`)

X is scalar (`JRank==0`), Y is tensor (`JRank=J_Y`). Z has the same rank as Y.

**IMSRG(2) `st` functions** (always on by default):

| Key | Function |
|---|---|
| `"comm111st"` | `comm111st(X,Y,Z)` |
| `"comm121st"` | `comm121st(X,Y,Z)` |
| `"comm221st"` | `comm221st(X,Y,Z)` |
| `"comm122st"` | `comm122st(X,Y,Z)` |
| `"comm222_pp_hhst"` | combined with 221 via `comm222_pp_hh_221st(X,Y,Z)` |
| `"comm222_phst"` | `comm222_phst(X,Y,Z)` |

**IMSRG(3) `st` n7 functions** (off by default, toggled by `SetUseIMSRG3N7_Tensor`):
`comm331st`, `comm231st`, `comm132st`, `comm232st`, `comm133st`, `comm223st`

**IMSRG(3) `st` full functions** (off by default, toggled by `SetUseIMSRG3_Tensor`):
`comm332_ppph_hhhpst`, `comm332_pphhst`, `comm233_pp_hhst`, `comm233_phst`, `comm333_ppp_hhhst`, `comm333_pph_hhpst`

All pass unit tests (marked `// PASS the unit test (J and T)` in `ReferenceImplementations.hh`).

**`tensor_transform_first_pass` mechanism:** On the first call for a given `(JRank, parities)` combination, `SetSingleThread(true)` is forced so that 6j/9j symbols are computed and cached without race conditions. Subsequent calls restore parallelism.

```cpp
if (Z.modelspace->tensor_transform_first_pass[Z.GetJRank()*4 + X.GetParity() + 2*Y.GetParity()])
    SetSingleThread(true);
// ... run all st terms ...
Z.modelspace->tensor_transform_first_pass.at(...) = false;
```

### 8.6 Tensor Commutators — File Layout

The `st` functions are **not in `Commutator.cc`** — they live in:
- `TensorCommutators.cc` — optimized production implementations
- `ReferenceImplementations.cc` — reference (loop-level) implementations

`Commutator.cc` also `#include`s `Commutator232.hh` (a separate header for the `comm232*` family).

Uses 6j symbols for recoupling. Works with **reduced matrix elements** (`is_reduced=true`).

### 8.7 Dagger Commutators

`DaggerCommutators.cc`: commutators involving dagger operators. Named by leg count (`sd` = scalar-dagger):

| Function | Formula | Notes |
|---|---|---|
| `comm211sd` | $[X^{(2)}, Y^{(1)}]^{(1)}$ | matrix × vector |
| `comm231sd` | $[X^{(2)}, Y^{(3)}]^{(1)}$ | monopole contraction |
| `comm413_233sd` | $[X^{(4)}, Y^{(1)}]^{(3)} + [X^{(2)}, Y^{(3)}]^{(3)}$ | combined |
| `comm431sd` | $[X^{(4)}, Y^{(3)}]^{(1)}$ | pp/hh ladder |
| `comm433sd_pphh` | $[X^{(4)}, Y^{(3)}]^{(3)}_{pp/hh}$ | matrix multiplication |
| `comm433sd_ph` | $[X^{(4)}, Y^{(3)}]^{(3)}_{ph}$ | Pandya transformation |

`comm433_pp_hh_431sd`: combines `comm431sd` and `comm433sd_pphh` by computing intermediate matrices once:
$$\mathcal{M}_{pp} = X^{(4)}\big|_{\text{cols}=pp} \cdot Y^{(3)}\big|_{\text{rows}=pp}, \quad \mathcal{M}_{hh} = X^{(4)}\big|_{hh} \cdot \text{diag}(n_a n_b) \cdot Y^{(3)}\big|_{hh}$$

### 8.8 Commutator Global Flags and `comm_term_on` Map

```cpp
namespace Commutator {
    bool single_thread;      // if true, disable OpenMP parallelism
    bool verbose;            // if true, print per-term timing info
    map<string, bool> comm_term_on;  // toggle individual terms on/off
}
```

Key control functions:
```cpp
Commutator::SetUseIMSRG3(true);           // include all IMSRG(3) ss terms
Commutator::SetUseIMSRG3N7(true);         // include N^7 scaling IMSRG(3) ss terms
Commutator::SetUseIMSRG3_MP4(true);       // include MP4-level 3b ss corrections
Commutator::SetUseIMSRG3_Tensor(true);    // include all IMSRG(3) st terms
Commutator::SetUseIMSRG3N7_Tensor(true);  // include N^7 IMSRG(3) st terms
Commutator::TurnOffTerm("comm220ss");     // disable a specific term
Commutator::TurnOnTerm("comm220ss");
Commutator::PrintSettings();              // dump all flags to stdout
```

Complete `comm_term_on` key list (default value shown):

**IMSRG(2) ss** (all `true`):
`comm110ss`, `comm220ss`, `comm111ss`, `comm121ss`, `comm221ss`, `comm122ss`, `comm222_pp_hhss`, `comm222_phss`

**IMSRG(2) st** (all `true`):
`comm111st`, `comm121st`, `comm221st`, `comm122st`, `comm222_pp_hhst`, `comm222_phst`

**IMSRG(3) ss n7** (all `false`):
`comm330ss`, `comm331ss`, `comm231ss`, `comm132ss`, `comm232ss`, `comm133ss`, `comm223ss`

**IMSRG(3) ss full** (all `false`):
`comm332_ppph_hhhpss`, `comm332_pphhss`, `comm233_pp_hhss`, `comm233_phss`, `comm333_ppp_hhhss`, `comm333_pph_hhpss`

**IMSRG(3) st n7** (all `false`):
`comm331st`, `comm231st`, `comm132st`, `comm232st`, `comm133st`, `comm223st`

**IMSRG(3) st full** (all `false`):
`comm332_ppph_hhhpst`, `comm332_pphhst`, `comm233_pp_hhst`, `comm233_phst`, `comm333_ppp_hhhst`, `comm333_pph_hhpst`

**Tensor-tensor (`tt`):** no keys exist — `CommutatorTensorTensor` is not implemented.

### 8.9 `ConstructScalarMpp_Mhh` — Intermediate Matrix for pp/hh Diagrams

```cpp
void ConstructScalarMpp_Mhh(const Operator &X, const Operator &Y, const Operator &Z,
                             TwoBodyME &Mpp, TwoBodyME &Mhh);
```

Builds intermediate matrices $M_{pp}$ and $M_{hh}$ used in both `comm222_pp_hhss` and (in some implementations) `comm222_phss`, avoiding redundant computation.

### 8.10 BCH Product Functions (in Commutator namespace)

```cpp
// These are used in BCH product e^X e^Y ~ e^Z:
void prod110ss(const Operator &X, const Operator &Y, Operator &Z); // 1b*1b -> 0b
void prod111ss(const Operator &X, const Operator &Y, Operator &Z); // 1b*1b -> 1b
void prod112ss(const Operator &X, const Operator &Y, Operator &Z); // 1b*2b (or 2b*1b) -> 2b
```

These compute the operator product $Z = X \cdot Y$ (not the commutator) and are used for the BCH expansion corrections beyond leading-order commutators.

---

## 9. `Generator` — Computing $\eta(s)$

The generator $\eta$ is anti-Hermitian and drives the flow toward decoupling. Different choices of $\eta$ correspond to different generator types.

### 9.1 Generator Types

| Type string | Formula | Properties |
|---|---|---|
| `"wegner"` | $\eta = [H_d, H_{od}]$ | Original Wegner, may hit poles |
| `"white"` | $\eta_{ia} = H_{ia} / \Delta_{ia}$ | White's approximation, fast |
| `"atan"` | $\eta_{ia} = \arctan(H_{ia}/\Delta_{ia})$ | Regularized, default choice |
| `"imaginarytime"` | $\eta_{ia} = -\text{sgn}(\Delta_{ia}) H_{ia}$ | Imaginary time evolution |
| `"qtransferatan1"` | momentum-transfer regulated | For open-shell systems |
| `"shell-model"` | drives off-diagonal core↔valence & qspace | VS-IMSRG |
| `"shell-model-3body"` | includes 3-body generator | VS-IMSRG(3) |
| `"hartree-fock"` | special HF generator | Initial basis optimization |
| `"1PA"` | one-particle approximation | EOM variant |

### 9.2 Denominator Partitioning

For White/arctan generators: $\eta_{ia} = f(\Delta_{ia})$ where

$$\Delta_{ia} = f_{ii} - f_{aa} \quad (\text{Möller-Plesset})$$
$$\Delta_{ia} = H_{ii} - H_{aa} \quad (\text{Epstein-Nesbet})$$

Set via `SetDenominatorPartitioning("Moller_Plesset")` or `"Epstein_Nesbet"`.

`denominator_delta`: level shift added to all denominators to prevent poles.  
`denominator_delta_index`: if set, applies the shift only near a specific orbit.

---

## 10. `BCH` Namespace — Magnus / BCH Transformation

```cpp
Operator BCH_Transform(Op, Omega)           // e^Omega H e^{-Omega}
Operator BCH_Product(X, Y)                  // e^X e^Y ~ e^Z via BCH series
Operator Standard_BCH_Transform(Op, Omega)  // standard nested commutator series
Operator Brueckner_BCH_Transform(Op, Omega) // Brueckner variant
```

The BCH transform is the series:

$$e^\Omega H e^{-\Omega} = H + [\Omega, H] + \frac{1}{2}[\Omega,[\Omega,H]] + \ldots$$

Controlled by thresholds:
```cpp
BCH::bch_transform_threshold  // stop when ||[Omega, ...term...]|| < threshold
BCH::bch_product_threshold
```

Optional corrections:
- `use_goose_tank_correction`: Goose-tank correction for large $\Omega$
- `use_brueckner_bch`: Brueckner resummation variant
- `use_factorized_correction`: factorized double-commutator approximation
- `only_2b_omega`: restrict $\Omega$ to 2-body even in IMSRG(3)

### 10.1 Factorized Double Commutator (`FactorizedDoubleCommutator` namespace)

An approximation to the double commutator $[\eta, [\eta, \Gamma]]$ that avoids full IMSRG(3) scaling:

```cpp
namespace Commutator::FactorizedDoubleCommutator {
    // Main routines (add to Z):
    void comm223_231(Eta, Gamma, Z);   // 1b contribution from [eta2, [eta2, Gamma2]]
    void comm223_232(Eta, Gamma, Z);   // 2b contribution from [eta2, [eta2, Gamma2]]

    // Intermediate variants (GooseTank corrections, TypeII/III):
    void comm223_231_chi1b(Eta, Gamma, Z);
    void comm223_231_chi2b(Eta, Gamma, Z);
    void comm223_232_chi1b(Eta, Gamma, Z);
    void comm223_232_chi2b(Eta, Gamma, Z);
    void comm223_132(Eta, Gamma, Z);   // cross term

    // Control flags (all default false, enable selectively):
    bool use_goose_tank_1b;     // Goose-tank 1-body correction
    bool use_goose_tank_2b;     // Goose-tank 2-body correction
    bool use_1b_intermediates;  // chi1b intermediates
    bool use_2b_intermediates;  // chi2b intermediates
    bool use_TypeII_1b;         // Type-II 1b correction
    bool use_TypeIII_1b;
    bool use_TypeII_2b;
    bool use_TypeIII_2b;
    bool use__GT_TypeI_2b;
    bool use__GT_TypeIV_2b;
}
```

This is activated by `BCH::use_factorized_correction = true` in `BCHSolver.hh`. It provides a computationally cheaper route to include some IMSRG(3) physics during the IMSRG(2) flow without carrying the full 3-body operator.

---

## 11. `IMSRGSolver` — ODE Driver

Controls the flow $ds$ integration and builds up $\Omega(s)$.

### 11.1 Key Fields

```cpp
Operator* H_0;             // initial Hamiltonian
deque<Operator> FlowingOps; // operators being simultaneously flowed
Operator Eta;              // current generator
deque<Operator> Omega;     // Magnus operator (may be split into multiple pieces)
Generator generator;       // generator prescription
double s, ds, smax;        // flow parameter, step size, maximum
double omega_norm_max;     // threshold: start new Omega when ||Omega|| > this
double eta_criterion;      // convergence criterion on ||eta||
string method;             // "magnus_euler", "magnus_modified_euler", "flow", ...
```

### 11.2 Integration Methods

| `method` string | Description |
|---|---|
| `"magnus_euler"` | Euler step on $\Omega$ |
| `"magnus_modified_euler"` | Midpoint/RK2 on $\Omega$ |
| `"flow"` | Direct flow (no Magnus accumulation) |
| `"flow_RK4"` | RK4 integration of flow equations |

### 11.3 Adaptive Step Size

When `magnus_adaptive=true`:
- Step size grows by at most `ds_max_growth_factor_ = 1.2` per step
- Backs off by `ds_backoff_factor_ = 0.5` if $|E_{MP2}(s+ds)| > |E_{MP2}(s)|$
- "Soft landing" phase: no growth once near convergence

### 11.4 Hunter-Gatherer Mode

When `hunter_gatherer=true`: instead of accumulating all $\Omega$ into one operator, multiple $\Omega$ pieces are kept and gathered (multiplied via BCH) when the norm exceeds `omega_norm_max`. Reduces approximation error from BCH truncation.

### 11.5 `FlowingOps` — Flowing Multiple Operators Simultaneously

The key feature: any number of additional operators can be flowed alongside the Hamiltonian at no extra BCH cost:

```cpp
IMSRGSolver solver(H_no);   // H_no goes into FlowingOps[0]
solver.AddOperator(R2);      // FlowingOps[1]
solver.AddOperator(GT_op);   // FlowingOps[2]
solver.Solve();
Operator H_flowed  = solver.GetOperator(0); // = GetH_s()
Operator R2_flowed = solver.GetOperator(1);
Operator GT_flowed = solver.GetOperator(2);
```

Each `AddOperator` call pushes onto `FlowingOps`. During the BCH transform at each step, **all** operators in `FlowingOps` are transformed with the same $\Omega$. The generator is always built from `FlowingOps[0]` (the Hamiltonian).

### 11.6 `IMSRGProfiler` and Scratch File Mechanism

`IMSRGProfiler` (in `IMSRGProfiler.hh`) tracks timing and FLOP counts per commutator term. Accessed via `solver.profiler`.

For very long flows where memory for multiple $\Omega$ pieces would be excessive, `FlushOmegaToScratch()` writes old $\Omega$ pieces to disk:

```cpp
solver.SetScratchDir("/scratch/tmp/");  // directory for omega dump files
// Or via ReadWrite:
solver.SetReadWrite(rw);  // rw.GetScratchDir() is used
```

During `Transform(Op)`, all on-disk and in-memory $\Omega$ pieces are applied sequentially.

### 11.7 Perturbative Triples Correction

```cpp
solver.SetPerturbativeTriples(true);  // adds (T) triples correction at end of flow
double E_T = solver.CalculatePerturbativeTriples();  // direct call
```

Computes the leading perturbative correction from 3-body diagrams using the converged IMSRG(2) $\Omega$.

---

## 12. `HartreeFock` — Basis Optimization

Performs a Hartree-Fock self-consistent field calculation to find the optimal single-particle basis before the IMSRG flow.

Key flow:
1. Build the HF Hamiltonian matrix from input $V_{NN}$ (and optionally $V_{3N}$ in NO2B approximation)
2. Diagonalize the Fock matrix iteratively
3. Rotate all operators into the HF basis via the unitary transformation $C$ (HF coefficients)

The HF transformation is:
$$f_{ij}^{HF} = \sum_{ab} C_{ia} f_{ab} C_{jb}^* + \sum_{abc, J} (2J+1) n_c\, C_{ia} C_{jb} \Gamma^J_{acbc}$$

After convergence, `UpdateEta()` sets $\eta$ to zero and restarts the IMSRG flow.

---

## 13. `EOM` — Equation-of-Motion Excited States

Computes excited states of the IMSRG-transformed Hamiltonian as:

$$|\Psi_\nu\rangle = \hat{O}_\nu |\Psi_0\rangle, \quad \hat{O}_\nu = \sum_i X_i^{(\nu)} a^\dagger_i + \sum_{ijk,J} X^{(\nu)}_{ij,k} (a^\dagger_i a^\dagger_j \tilde{a}_k)^J + \ldots$$

The EOM eigenvalue problem is solved in the space of 1-particle-1-hole ($1p1h$) and $2p1h$ excitations. Uses dagger operators (`ThreeLegME`) and the dagger commutators for computing the EOM matrix elements.

---

## 14. Angular Momentum Algebra (`AngMom` namespace)

```cpp
double CG(ja, ma, jb, mb, J, M)      // Clebsch-Gordan coefficient
double ThreeJ(j1,j2,j3, m1,m2,m3)    // Wigner 3j symbol
double SixJ(j1,j2,j3, J1,J2,J3)      // Wigner 6j symbol
double NineJ(j1,...,j9)               // Wigner 9j symbol
double NormNineJ(...)                 // normalized 9j: sqrt((2J1+1)(2J2+1)) * NineJ
double TwelveJ_1(...)                 // 12j symbol of the first kind
double Moshinsky(N,L,n,l,n1,l1,n2,l2,lam,B)  // Moshinsky bracket
bool   Triangle(j1, j2, j3)          // triangle inequality check
int    phase(x)                       // (-1)^x
```

`Moshinsky` bracket transforms from lab-frame $(n_1 l_1, n_2 l_2)$ to relative+CM $(nL, Nlam)$ with mass ratio parameter `B` (`mosh_beta_1 = pi/4` for equal masses).

`AngMomCache` / `SixJCache_112112`: specialized cache for the 6j pattern $\{j_1\,j_2\,J; j_3\,j_4\,J'\}$ with $(l_1, l_2)$ fixed at 1/2 steps apart — the dominant pattern in two-body commutators.

---

## 15. Physical Constants (`PhysConst` namespace)

All in MeV / fm units:

```cpp
HBARC         = 197.3269718    // MeV·fm
M_PROTON      = 938.2720813    // MeV/c²
M_NEUTRON     = 939.5654133    // MeV/c²
M_NUCLEON     = (Mp+Mn)/2
M_PION_CHARGED= 139.57018      // MeV/c²
NUCLEON_AXIAL_G = 1.27         // gA
PROTON_SPIN_G = 5.585690569
NEUTRON_SPIN_G= -3.826085
ALPHA_FS      = 1/137.035999   // fine structure constant
F_PI          = 92.2           // MeV, pion decay constant
PROTON_RCH2   = 0.707          // fm², from PDG 2020
NEUTRON_RCH2  = -0.1161        // fm², from PDG 2020
SQRT2, INVSQRT2, PI, SQRTPI
```

---

## 16. Key Numerical Conventions Summary

| Convention | Value/Rule |
|---|---|
| Half-integer quantum numbers | Stored as `j2 = 2j` (integer) |
| Proton | `tz2 = -1` |
| Neutron | `tz2 = +1` |
| Hole orbit | `occ = 1` |
| Particle orbit | `occ = 0` |
| Canonical 2-body ordering | `p <= q` (smaller index first) |
| Canonical 3-body ordering | `a >= b >= c` (larger index first) |
| Exchange phase | $\vert pqJ\rangle = (-1)^{j_p+j_q-J+1}\vert qpJ\rangle$ |
| Tilde normalization | $\tilde{\Gamma}_{ijkl} = \sqrt{(1+\delta_{ij})(1+\delta_{kl})}\,\Gamma_{ijkl}$ |
| Wigner-Eckart | Bra normalization: divide by $\sqrt{2J+1}$ |
| Scalar ME | Non-reduced (full $\langle JM\vert O\vert JM\rangle$) |
| Tensor ME | Reduced ($\langle J\Vert O^\lambda\Vert J'\rangle$) |
| OCC_CUT | $10^{-6}$: skip if $\vert n_a\bar{n}_b\vert$ below this |
| Energy unit | MeV throughout |
| Frequency | `hbar_omega` in MeV |

---

## 17. Quick Code Patterns for AI

### Reading a 2-body matrix element

```cpp
// By orbit indices (most readable):
double v = Op.TwoBody.GetTBME_J(J, J, a, b, c, d);   // tilde, sum over Tz/parity

// By channel + orbit indices:
int ch = modelspace.GetTwoBodyChannelIndex(J, parity, Tz);
double v = Op.TwoBody.GetTBME(ch, ch, a, b, c, d);    // tilde

// Fast inner loop (by local ket index):
TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
int ibra = tbc.GetLocalIndex(min(a,b), max(a,b));
int iket = tbc.GetLocalIndex(min(c,d), max(c,d));
double v = Op.TwoBody.GetTBME_norm(ch, ch, ibra, iket); // normalized
```

### Looping over a two-body channel

```cpp
int nch = modelspace.GetNumberTwoBodyChannels();
for (int ch = 0; ch < nch; ch++) {
    TwoBodyChannel& tbc = modelspace.GetTwoBodyChannel(ch);
    int J = tbc.J;
    int npq = tbc.GetNumberKets();
    for (int ibra = 0; ibra < npq; ibra++) {
        Ket& bra = tbc.GetKet(ibra);
        int i = bra.p, j = bra.q;
        for (int iket = ibra; iket < npq; iket++) {
            Ket& ket = tbc.GetKet(iket);
            int k = ket.p, l = ket.q;
            double me = Op.TwoBody.GetTBME_norm(ch, ch, ibra, iket);
        }
    }
}
```

### Normal ordering a Hamiltonian

```cpp
ModelSpace ms(emax, "O16", "sd-shell");
ms.SetHbarOmega(20.0);
// ... read in H ...
Operator H_bare(ms);   // scalar 2-body
// ... fill H_bare ...
Operator H_no = H_bare.DoNormalOrdering();  // wrt ms.holes
double E0 = H_no.ZeroBody;
```

### Setting up IMSRG flow

```cpp
IMSRGSolver solver(H_no);
solver.SetMethod("magnus_euler");
solver.generator.SetType("atan");
solver.SetDs(0.05);
solver.SetSmax(500.0);
solver.SetOmegaNormMax(0.25);
solver.Solve();
double E_gs = solver.GetH_s().ZeroBody;
```

### Flowing a transition operator

```cpp
Operator GT(ms, 1, 1, 0, 2); // rank_J=1, rank_T=1, parity=0, 2-body
// fill GT ...
solver.AddOperator(GT);
solver.Solve();
Operator GT_flowed = solver.GetOperator(0);
// GT_flowed contains reduced MEs <J||GT||J'> at s=smax
```

---

## 18. File I/O (`ReadWrite`)

`ReadWrite rw;` is instantiated first; all methods take operator and filename arguments.

### 18.1 Reading Two-Body Interactions

```cpp
// Darmstadt ME2J format (most common for SRG-evolved NN+3N interactions):
rw.ReadBareTBME_Darmstadt(filename, H_bare, E1max, E2max, lmax);

// ME2J gzipped:
// (automatically detected from .gz extension)

// Navratil format:
rw.ReadBareTBME_Navratil(filename, H_bare);

// Oslo group format:
rw.ReadTBME_Oslo(filename, H_bare);

// Jason format / Nathan format:
rw.ReadBareTBME_Jason(filename, H_bare);
rw.ReadTensorOperator_Nathan(file1b, file2b, op);
rw.ReadOperator_Nathan(file1b, file2b, op);
```

### 18.2 Reading Three-Body Interactions

```cpp
// Darmstadt 3N format:
rw.Read_Darmstadt_3body(filename, H_bare, E1max, E2max, E3max);

// HDF5 format (requires HDF5 support at compile time):
rw.Read3bodyHDF5(filename, op);
rw.Read3bodyHDF5_new(filename, op);
```

### 18.3 Writing Matrix Elements

```cpp
// ME2J format:
rw.Write_me2j(filename, op, emax, e2max, lmax);

// Three-body ME3J format:
rw.Write_me3j(filename, op, E1max, E2max, E3max);

// Shell-model codes:
rw.WriteNuShellX_int(op, filename);   // NuShellX interaction file
rw.WriteNuShellX_op(op, filename);    // NuShellX operator file
rw.WriteNuShellX_sps(op, filename);   // NuShellX single-particle file
rw.WriteAntoine_int(op, filename);    // ANTOINE format
rw.WriteAntoine_input(op, filename, A, Z);

// Valence 3-body for VS-IMSRG:
rw.WriteValence3body(threeBME, filename);
```

### 18.4 Generic Operator I/O (internal binary format)

```cpp
// Fast binary dump/restore of full Operator (for checkpointing):
rw.WriteOperator(op, filename);       // binary format
rw.ReadOperator(op, filename);
rw.WriteOperatorHuman(op, filename);  // text format (slow, for debugging)
rw.ReadOperatorHuman(op, filename);

rw.WriteTensorOneBody(filename, op, opname);
rw.WriteTensorTwoBody(filename, op, opname);
rw.WriteDaggerOperator(op, filename, opname);
```

### 18.5 Takayuki / Osaka Formats

```cpp
rw.ReadOneBody_Takayuki(filename, op);
rw.ReadTwoBody_Takayuki(filename, op);
rw.WriteOneBody_Takayuki(filename, op);
rw.WriteTwoBody_Takayuki(filename, op);
```

### 18.6 Reading Shell-Model Files

```cpp
rw.ReadNuShellX_int(op, filename);    // read NuShellX interaction
rw.ReadNuShellX_int_iso(op, filename);
rw.ReadNuShellX_sp(ms, filename);     // read single-particle space
```

---

## 19. `imsrg_util` — High-Level Operator Builders

The `imsrg_util` namespace (in `namespace imsrg_util`) provides factory functions for constructing all commonly needed operators. These are the main user-facing API.

### 19.1 Kinetic and CM Operators

```cpp
Operator Trel_Op(ModelSpace& ms);                 // relative kinetic energy T_rel
Operator TCM_Op(ModelSpace& ms);                  // CM kinetic energy T_CM
Operator HCM_Op(ModelSpace& ms);                  // HO CM Hamiltonian (for Lawson term)
Operator KineticEnergy_Op(ModelSpace& ms);         // p^2/2m
Operator KineticEnergy_RelativisticCorr(ms);      // relativistic T correction
Operator Trel_Masscorrection_Op(ms);
```

### 19.2 Radius and Density Operators

```cpp
Operator R2CM_Op(ModelSpace& ms);                 // <R_CM^2>
Operator Rp2_corrected_Op(ms, A, Z);             // proton radius^2 with finite-size correction
Operator Rn2_corrected_Op(ms, A, Z);             // neutron radius^2
Operator Rm2_corrected_Op(ms, A, Z);             // matter radius^2
Operator R2_p1_Op(ms);                           // one-body proton r^2
Operator R2_1body_Op(ms, option);                // "proton", "neutron", or "matter"
Operator R2_p2_Op(ms);                           // two-body correction to r^2
Operator R2_2body_Op(ms, option);
Operator DensityAtR(ms, R, pn);                  // density at radius R ("proton"/"neutron")
Operator FormfactorAtQ(ms, q, pn);               // form factor at momentum transfer q
Operator FourierBesselCoeff(ms, nu, R, index_list); // Fourier-Bessel coefficient
Operator RpSpinOrbitCorrection(ms);              // Darwin-Foldy spin-orbit correction
Operator E0Op(ms);                               // E0 operator
```

### 19.3 Electromagnetic Operators

```cpp
Operator ElectricMultipoleOp(ms, L);             // E_L operator (proton only)
Operator NeutronElectricMultipoleOp(ms, L);      // neutron E_L
Operator IntrinsicElectricMultipoleOp(ms, L);    // intrinsic E_L (CM removed)
Operator MagneticMultipoleOp(ms, L);             // M_L  (uses g_s,g_l)
Operator MagneticMultipoleOp_pn(ms, L, pn);     // proton or neutron separately
Operator IVDipoleOp(ms, rL, YL);                // isovector dipole
Operator ISDipoleOp(ms, rL, YL, Rms);           // isoscalar dipole
Operator SchiffOp(ms, rL, YL, Rms);             // Schiff moment operator
```

### 19.4 Gamow-Teller and Weak Operators

```cpp
Operator AllowedFermi_Op(ms);                   // Fermi transition operator
Operator AllowedGamowTeller_Op(ms);             // Gamow-Teller operator sigma*tau
Operator Sigma_Op(ms);                          // spin operator sigma
Operator Sigma_Op_pn(ms, pn);
Operator AxialCharge_Op(ms);                    // axial charge
```

### 19.5 Angular Momentum Operators

```cpp
Operator Isospin2_Op(ms);       // T^2
Operator TzSquared_Op(ms);      // T_z^2
Operator PSquaredOp(ms);        // total momentum P^2
Operator LdotS_Op(ms);          // L·S spin-orbit operator
Operator L2rel_Op(ms);          // L_rel^2
Operator LCM_Op(ms);            // CM angular momentum L_CM
Operator QdotQ_Op(ms);          // Q·Q quadrupole-quadrupole
Operator VQQ_Op(ms);            // Q·Q potential (pi/2 phase)
```

### 19.6 Occupation Number Operators

```cpp
Operator NumberOp(ms, n, l, j2, tz2);       // number in a specific orbit
Operator NumberOpAlln(ms, l, j2, tz2);      // sum over n
Operator NumberOpRef(ms);                    // total particle number
Operator  OneBodyDensity(ms, i, j);          // one-body density a†_i a_j
```

### 19.7 Potentials

```cpp
Operator VCoulomb_Op(ms, lmax);              // Coulomb potential
Operator VCentralCoulomb_Op(ms, lmax);       // central part of Coulomb
Operator TViolatingPotential_Op(ms, LECs);   // T-violating potential
Operator WoodsSaxon1b_Op(ms, V0, R0, a0);   // 1-body Woods-Saxon
Operator HOtrap_Op(ms, hw_trap);             // HO trap potential (for cold atoms)
```

### 19.8 String-Based Operator Constructor

```cpp
Operator OperatorFromString(ms, str);
// str can be: "Rp2", "Rm2", "R2CM", "GT", "Fermi", "E1", "E2", "M1", "Dipole",
//             "TCM", "Trel", "Sigma", "LdotS", "Isospin2", "R2p1", etc.
// Use this in scripts/Python to build operators by name.
```

### 19.9 Perturbation Theory (built into `imsrg_util`)

```cpp
double MBPT2_SpectroscopicFactor(Operator H, index_t p);
// Returns the leading-order spectroscopic factor for orbit p
```

### 19.10 Atomic Physics Operators (`imsrg_util::atomic_fs`, `::atomic_hfs`)

```cpp
namespace atomic_fs {
    Operator Darwin(ms, Z);          // Darwin term (relativistic)
    Operator RelativisticT(ms);      // relativistic kinetic energy
    Operator SpinOrbit(ms, Z);       // atomic spin-orbit
}
namespace atomic_hfs {
    Operator hQ(ms);                 // electric quadrupole hyperfine
    Operator hD(ms);                 // magnetic dipole hyperfine
    Operator NormalMassShift(ms, A);
    Operator SpecificMassShift(ms, A);
    Operator CombinedMassShift(ms, A);
}
```

---

