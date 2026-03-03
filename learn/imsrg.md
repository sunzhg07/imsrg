# imsrg++ — AI-Friendly Codebase Reference

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

### 8.4 Tensor Commutators

`TensorCommutators.cc`: commutator terms for tensor operators (rank $\Lambda \neq 0$). Uses 6j symbols for recoupling:

$$[X^{(\Lambda)}, Y^{(0)}]^{(\Lambda)} \ni \sum_{J_1 J_2} \sqrt{(2J_1+1)(2J_2+1)} \begin{Bmatrix} J_1 & J_2 & \Lambda \\ j_j & j_i & j_a \end{Bmatrix} X^{J_1}_{ab} Y^{J_2}_{cd}$$

Works entirely with **reduced matrix elements** (`is_reduced=true`).

### 8.5 Dagger Commutators

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

### 8.6 Commutator Global Flags

```cpp
namespace Commutator {
    bool single_thread;      // if true, disable OpenMP parallelism (useful for debugging)
    bool verbose;            // if true, print per-term timing info
    map<string, bool> comm_term_on;  // toggle individual terms on/off
}
```

Key control functions:
```cpp
Commutator::SetUseIMSRG3(true);        // include all IMSRG(3) terms
Commutator::SetUseIMSRG3N7(true);      // include N^7 scaling IMSRG(3) terms
Commutator::SetUseIMSRG3_MP4(true);    // include MP4-level 3b corrections
Commutator::TurnOffTerm("comm220ss");  // disable a specific term for debugging/testing
Commutator::TurnOnTerm("comm220ss");
Commutator::PrintSettings();
```

Terms that can be toggled via `comm_term_on`: `"comm110ss"`, `"comm220ss"`, `"comm111ss"`, `"comm121ss"`, `"comm221ss"`, `"comm122ss"`, `"comm222_pp_hhss"`, `"comm222_phss"`, `"comm330"`, etc.

### 8.7 `ConstructScalarMpp_Mhh` — Intermediate Matrix for pp/hh Diagrams

```cpp
void ConstructScalarMpp_Mhh(const Operator &X, const Operator &Y, const Operator &Z,
                             TwoBodyME &Mpp, TwoBodyME &Mhh);
```

Builds intermediate matrices $M_{pp}$ and $M_{hh}$ used in both `comm222_pp_hhss` and (in some implementations) `comm222_phss`, avoiding redundant computation.

### 8.8 BCH Product Functions (in Commutator namespace)

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
