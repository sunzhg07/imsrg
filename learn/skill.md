# IMSRG EOM Codebase — AI Knowledge Base

## 1. Project Layout
```
src/          C++ source (EOM, Commutator, ModelSpace, Operator, HartreeFock, IMSRGSolver, AngMom, …)
run/          Python scripts (mr_eom.py, lanczos.py, sr_eom.py, …)
build/        CMake build output; pyIMSRG.so = pybind11 module loaded by Python scripts
extern/       armadillo, boost, pybind11, half
```
Entry point: `src/imsrg++.cc` → HF → IMSRGSolver → EOM.

Build: `cd build && make -j8`  (force recompile: `touch ../src/FILE.cc && make -j8`)

Python env: `/Users/wolf/work/imsrg/.venv`; `PYTHONPATH=/Users/wolf/work/imsrg/build`

---

## 2. EOM Class (`src/EOM.cc`, `src/EOM.hh`)

### Constructors
| Signature | Mode |
|---|---|
| `EOM(Hs, rdm, J2, parity, itz)` | multi-reference (MR) |
| `EOM(Hs, tdm_file, J2, parity, itz)` | MR, reads TDM from file |
| `EOM(Hs, J2, parity, itz)` | single-reference (SR) |

### Configuration spaces (`ConstructConfigs`)
Five sectors stored flat in `eom_confs` (index `eom_dims`):
- `qv` — valence qp (cvq=2 → cvq=1)
- `ph` — particle-hole
- `ppvv` — 2p2h with two valence holes
- `pphv` — 2p2h mixed
- `pphh` — 2p2h with two core holes

`eom_confs` element type: `{index_t orb_i, index_t orb_j, int rank, int flat_index}`

`cvq` orbit values: `0` = core, `1` = valence, `2` = excited/virtual (above valence).

All sectors must be built before `ConstructNormMatrix` / `ConstructProjectMatrix`.

### Norm kernel (`ConstructNormMatrix`, `ComputeNorm`)
`Nkernel` (arma::mat, `eom_dims × eom_dims`) — overlap matrix restricted to `eom_confs` indices.
`ComputeNorm(Op1, Op2)` = `v1.t() * Nkernel * v2` using only `eom_confs` elements.

**Key fact**: `NormMultiref(Op1, Op2)` uses the full operator algebra path
(`GetVSEOM_ladder_multiref` + commutator + `GetVSEOM_Overlap_multiref`).
`ComputeNorm` and `NormMultiref` are **equivalent** — confirmed by debug run.
Both go negative at the same step for He8, proving the metric is **genuinely indefinite**
(not a construction bug) because He8 has fractional occupations (open-shell MR).

### Projection (`ConstructProjectMatrix`, `ProjectOprator`, `SqrtMat`)
`Prj_kernel` (arma::sp_mat) — projector onto positive-norm subspace of `Nkernel` blocks.
`SqrtMat` uses `arma::eig_sym` (NOT `arma::svd`):
- SVD returns `|λ|` → cannot distinguish negative eigenvalues → ghost directions survive in projector
- `eig_sym` keeps only `λ ≥ rel_tol * λ_max` (positive eigenvalues only)
- Old SVD version is **commented out** (not deleted) in `EOM.cc` with bug explanation, for reference.

`SqrtMat` zero-threshold: `1e-10` (was `1e-14` in original, changed by user to `1e-10`).

### `GetVSEOM_ladder` herm parameter
`herm=0` = anti-hermitian (sets operator to anti-hermitian before projection)
`herm=1` = hermitian (sets operator to hermitian before projection)
The correct value for the EOM amplitude operator is **`herm=1`**.
Using `herm=0` was the original bug; it was fixed in `HtcSingle`, `HtcMultiref`, `htc_single` (lanczos.py), and the `RunSR` initial random vector.

### Hamiltonian action
- SR: `HtcSingle(Hs, chi)` = `GetVSEOM_ladder_single(Commutator(Hs,chi), herm=1)`
- MR: `HtcMultiref(Hs, chi)` = `GetVSEOM_ladder_multiref(Commutator(Hs,chi), herm=1)`

### Diagonal commutator
`DcomMultiref(Hs, v)` — returns `(double, Operator)` pair; `.first` = diagonal H matrix element.

---

## 3. Solvers

### `LanczosSolve` (SR-EOM)
Standard real-symmetric Lanczos. Key fixes applied:
1. **Convergence flag**: `bool converged` — when set, skip the final full-subspace re-solve.
   Without this, the extra basis vector pushed at the last step causes `submat` OOB or introduces
   spurious low eigenvalues (e.g. `-16.706` when true eigenvalue is `3.316`).
2. **Bounds guard**: `dim = std::min((int)lanczos_vector.size(), max_iter)` in both the
   convergence-check submat and the final-solve submat.
   Root cause: loop runs `j = 0..max_iter-1`; at last step a new vector is pushed →
   `lanczos_vector.size() = max_iter+1` → `submat(0,0,max_iter,max_iter)` = OOB crash.

### `ArnoldiSolve` (MR-EOM)
Non-symmetric Krylov with indefinite metric. Current behaviour:
- Builds `hall` (Hamiltonian matrix) via **polarization identity**:
  `H_{ij} = 0.5*(NormMultiref(v_i, H1*v_j) + NormMultiref(v_j, H1*v_i)) + 0.5*(DcomMultiref(Hs, v_i+v_j) - diag_i - diag_j)`
  where `H1` part uses cached `h1v_cache` (saves one `HtcMultiref` call per pair) and `Dcom` part handles the double-commutator contribution.
- Double-pass CGS with `ProjectOprator` applied after each pass.
- **Norm guards** after GS:
  - `|ComputeNorm(w,w)| < null_tol * cn0` → null vector / breakdown, stop.
  - `|NormMultiref(w,w)| < bj_tol (1e-10)` → exact breakdown, stop.
  - `NormMultiref(w,w) < 0` → indefinite metric hit, solve subspace, stop.
- `null_tol = 1e-6` (relative, referenced to `cn0 = ComputeNorm(v0,v0)`).
- Restart was tried and removed — did not help because the metric itself is indefinite.

### Known root issue (He8 MR-EOM)
`bj = NormMultiref(w,w) → −119248` at step 8 (or step 28 with guards active).
Both `NormMultiref` and `ComputeNorm` give the same negative value → `Nkernel` is
correct; the MR norm form `<[T†,T]>_ρ` is genuinely indefinite for He8 open-shell reference.
True solution requires one of:
- **Generalized EVP** `H x = e N x` directly (dense `eom_dims × eom_dims`, feasible for He8 with ~738 configs). `Nkernel` and `Hkernel` are both available after `ConstructNormMatrix`. Use `arma::eig_gen` or `arma::eigs_gen`.
- **Biorthogonal (left/right) Lanczos** for indefinite metric — more complex to implement
- **Restrict to positive-definite sector** (drop `pphh`/`pphv` channels that contribute negative blocks)

Note: `eom_dims` for He8 p-shell = 738 (4 ph + 395 ppvv + 276 pphv + 66 pphh + 1 qv), so dense EVP is tractable.

---

## 4. pybind11 Bindings (`src/pyIMSRG.cc`)

EOM section (lines ~1185–1215). Currently exposed methods:
```python
EOM(Hs, rdm, J2, parity, itz)      # MR constructor
EOM(Hs, tdm_file, J2, parity, itz) # MR from file
EOM(Hs, J2, parity, itz)           # SR constructor
eom.Run(max_iter=200, state_want=6)
eom.ConstructConfigs()
eom.PrintConfigs()
eom.ConstructNormMatrix()
eom.ConstructProjectMatrix()
eom.force_decouple(H)
eom.ProjectOprator(Qin)
eom.GetVSEOM_Overlap_single(H1, H2)
eom.GetVSEOM_Overlap_multiref(H)
eom.GetVSEOM_ladder_single(H, herm)
eom.GetVSEOM_ladder_multiref(H, herm)
eom.ComputeNorm(Op1, Op2)    # added: Nkernel-based inner product
eom.NormMultiref(Op1, Op2)   # added: full operator algebra inner product
```

---

## 5. Live Code Changes vs HEAD (`git diff HEAD`)

### `src/EOM.cc` (13 lines changed)
1. `NormSingle`: old commented-out alternative implementation left as reference.
2. `HtcSingle`: `herm=0` → `herm=1`
3. `RunSR` initial vector: `herm=0` → `herm=1`
4. `LanczosSolve`: `bool converged` flag; skip final re-solve if converged; `dim` bounds cap.
5. `ArnoldiSolve`: `ComputeNorm` null-vector guard; `bj < 0` indefinite-metric stop.
6. `SqrtMat`: `s_max < 1e-10` threshold (was `1e-14`).

### `run/lanczos.py` (1 line changed)
- `htc_single`: `GetVSEOM_ladder_single(ht_plus, 0)` → `GetVSEOM_ladder_single(ht_plus, 1)`

### `src/pyIMSRG.cc`
- `ComputeNorm` and `NormMultiref` bindings added to EOM section.
- `PrintConfigs` binding added.

---

## 6. Debugging Workflow

```bash
# build
cd /Users/wolf/work/imsrg/build
touch ../src/EOM.cc && make -j8 2>&1 | grep -E "error:|EOM\.cc|pyIMSRG" | grep -v "ignoring"

# run MR-EOM
cd /Users/wolf/work/imsrg/run
source ../.venv/bin/activate
python3 mr_eom.py 2>&1 | tee /tmp/eom_debug.txt

# filter relevant output
grep -E "arnoldi|lanczos|DEBUG|dimension EOM|E_ref|E\(" /tmp/eom_debug.txt
```

Key diagnostic prints already in code:
- `dimension EOM ph/ppvv/pphv/pphh/all: ...` — from `ConstructConfigs`
- `arnoldi eigenvalues @ step N: ...` — every 5 steps
- `arnoldi: bj=X < 0 at step N (indefinite metric), stopping.` — metric failure
- `arnoldi: null vector (breakdown) at step N ...` — physical space exhausted
- `lanczos: energy converged` / `lanczos: exact breakdown at step N`

---

## 7. Architecture Summary

```
Python (mr_eom.py)
  └─ pyIMSRG.so (pybind11)
       └─ eom.Run(max_iter, state_want)
            ├─ dispatches to RunMR (is_multiref=true) or RunSR (is_multiref=false)
            │
            ├─ RunMR:
            │    ├─ ConstructConfigs       → eom_confs[], eom_dims
            │    ├─ ConstructNormMatrix    → Nkernel (arma::mat, eom_dims×eom_dims)
            │    ├─ ConstructProjectMatrix → Prj_kernel (arma::sp_mat) via SqrtMat(eig_sym)
            │    ├─ random projected initial vector (herm=1)
            │    └─ ArnoldiSolve           ← HtcMultiref, NormMultiref, DcomMultiref, ProjectOprator
            │
            └─ RunSR:
                 ├─ random initial vector GetVSEOM_ladder_single(random_op, herm=1)
                 └─ LanczosSolve          ← HtcSingle, NormSingle(=GetVSEOM_Overlap_single)
```

`NormSingle` is a thin wrapper: `return GetVSEOM_Overlap_single(T1, T2)` — no commutator path, purely overlap.
