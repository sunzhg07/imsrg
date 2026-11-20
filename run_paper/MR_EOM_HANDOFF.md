# MR-EOM benchmarks: status and open tasks (handoff, 2026-11-20)

Paper: github.com:sunzhg07/EOM-srg-doc, branch `mr-eom-revision` (`main.tex`).
Code: this repo, branch `devel`. All scripts below run from `run_paper/` with
`PYTHONPATH=../build` (pyIMSRG must be built for the Python you use).

Setup for every run: EM1.8/2.0, hw = 16 MeV, emax = 3, VS-IMSRG(2) with the
nucleus itself as reference, p-shell valence above a He4 core. References are
KSHELL p-shell eigenstates; their 1-, 2-, 3-body densities come from the
`kshell_threebody` transit build.

## Pipeline

| step | command |
|---|---|
| IMSRG + densities + MR-EOM | `NUCLEUS=He8\|Be8\|He6 PHASE=imsrg\|kshell\|eom python3 mr_rho3_full.py` (env `RDM3=1 MAX_ITER STATE_WANT JOBS KSHELL3_BIN`) |
| level-scheme figures | `NUCLEUS=He8\|Be8 python3 plot_mr_rho3.py` |
| He6 neutral bound check | `python3 compare_he6_fci.py` (needs `he6_fci/` logs) |

Outputs: `eom_result/{Nuc}_EM1.8_2.0_hw16_emax3_rho3_full.json`,
`*_rho3_levels.{pdf,png}`; densities in `run2/{nuc}_rho3/`; KSHELL work in
`kshell_{nuc}_rho3/e3/`; Hamiltonians in `hs_cache/Hs_{Nuc}_pshell_*_emax3.{dat,snt}`.

## Verified (do not redo)

1. `<Hbar>_rho` (rho1, rho2) = KSHELL eigenvalue + E_core to < 5e-6 MeV for all
   12 references (He8, Be8, He6; 0+_1, 0+_2, 2+_1, 2+_2). Hbar has no 3-body
   part, so rho3 changes it by exactly 0. Script: `rho3_test/check_eref.py`.
2. Norm kernel with rho3 is positive semidefinite (min eigenvalue ~1e-10) in all
   120 channels per nucleus (He8, Be8; 4 refs x J=0-4 x parity x dTz=0,1,2).
   Without rho3, 80/120 (He8) and 97/120 (Be8) channels have negative modes.
   Script: `rho3_test/norm_psd.py`.
3. Norm kernel = commutator `1/2 <[S_a+, S_b-]>` (rho3 off and on), relative
   diff <= 7e-13, for He8 gs ref: scalar and all 24 tensor channels
   (J=1-4, P=0,1, dTz=0,1,2). Log: `rho3_test/he8_gs_tensor_norm_summary.log`,
   script `rho3_test/test_norm_quick.py`.
   Use E3MAX=3 (not 6): only (p-shell)^3 kets meet rho3; identical numbers,
   ~65x faster.
4. sigma = ladder_{+1}([H, S]) is identical for S+ and S- (difference exactly
   0) in scalar, tensor, parity- and charge-changing channels: the decoupled H
   has no excitation-class components. Script: `rho3_test/test_sigma_spm.py`.
5. Arnoldi: every Ritz value equals the directly evaluated `(H2+H3)/N` to
   1e-13 with unit norm (all He8/Be8/He6 runs).
6. He6 -> Li6 (charge-changing, dTz=1, J=0-3): MR-EOM reproduces all six
   Li6 p-shell shell-model energies exactly (rel. to He6 gs: 1+ -3.95962,
   3+ 0.40792, 0+ 0.42353, 2+ 3.12422, 4.07629, 1+ 5.52042). With two valence
   nucleons rho3 = rho4 = 0, so the kernel is exact and the pool spans these
   states: this validates the whole chain (tensor norm, sigma, factorized H3,
   polarization, solver).

## Open tasks

### A. He6 neutral benchmark (main open item)
MR-EOM side is done: `eom_result/He6_EM1.8_2.0_hw16_emax3_rho3_full.json`
(8 neutral jobs; lowest roots, absolute MeV: 0+ from 0+_1 ref -20.79,
from 2+_1 ref -17.91; 2+ from 0+_1 ref -20.86, from 2+_1 ref -17.96; not
converged at 100 steps, residuals 0.01-1 MeV). Only the FCI side is missing.
Full-space diagonalization of the He6 evolved H, then the Ritz bound check.
- `cd he6_fci && sbatch he6_fci.slurm` (see `he6_fci/README.md`; KSHELL must
  be a standard build such as kshell-20230714-ver4: the `kshell_threebody`
  build stops after orbit 15 of the 20-orbit space). Dimension 3.86e6 (M=0).
- Then `python3 compare_he6_fci.py`. Expect: every MR-EOM root >= the k-th
  non-valence FCI level of the same J. A violation is a bug. The gap measures
  the incompleteness of the two-body operator pool.
- KSHELL does not add the `.snt` zero-body term; the compare script does.

### B. Decide the cause of the He8/Be8 discrepancy
MR-EOM roots lie 3-5 MeV (He8) and 1-4 MeV (Be8) above the lowest
non-valence evolved-H FCI states; He8 -> Li8 misses the Li8 p-shell states by
4-6 MeV. Two candidates remain (no bug, see Verified 6):
(i) rho3/rho4 contractions of [S_A, K3] (the induced three-body
intermediate), which are dropped; (ii) the two-body pool (Li8 states need
n->p inside the valence space, reachable only via pair transfer) plus the
H2-only Krylov space. Task A separates them: if He6 neutral roots are close to
FCI, (i) dominates in He8/Be8.

### C. Convergence of neutral runs
All He8/Be8 neutral runs stopped at MAX_ITER=100 with residual estimates
0.01-0.8 MeV. Rerun with `MAX_ITER=300` on the larger machine and update the
table/figures (`plot_mr_rho3.py`).

### D. Optional extra norm checks
Tensor kernel=commutator was run only on the He8 gs reference. Repeat for the
Be8 gs and the 2+ references with `rho3_test/test_norm_quick.py EMAX REF He4
p-shell 3 J P dTz` (E3MAX=3!).

### E. Paper (EOM-srg-doc, branch mr-eom-revision)
- Add the He6 results: the exact He6 -> Li6 test (Verified 6) and Task A.
- Results/Summary already state the two candidates of Task B; update once A
  is done.
- Abstract still describes the goal, not the outcome.
- He8 figure has an SR-EOM column from a different (closed-shell He8)
  Hamiltonian; keep or drop.
- `Iterative solution` subsection comes after the MR section but starts with
  the SR Lanczos case; consider order.

## Gotchas
- `test_rho3_norm.py` args: `EMAX REF CORE VAL E3MAX J2 PARITY ITZ`; pass them
  as separate arguments (a quoted "J P T" string silently falls back to J=0).
- FCI logs in `run_paper/log_*_fci_heb_{hf,hs}_*`: add the `.snt` zero-body
  term (`fci_heb_*.snt` header) to get absolute energies.
- Old `Hs_He8_FCI_...snt` / `hs` FCI for Be8 belong to the same Hs as the
  current runs (checked: valence states reappear to 1e-5 MeV).
