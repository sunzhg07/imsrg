# He6 full-space benchmark for MR-EOM

Full-space (emax=3, no many-body truncation) diagonalization of the He6
VS-IMSRG Hamiltonian, the same operator MR-EOM uses
(`hs_cache/Hs_He6_pshell_EM1.8_2.0_hw16_emax3.dat`, EM1.8/2.0, hw=16).

Why He6: with two valence neutrons rho3 = rho4 = 0, so the MR-EOM Hamiltonian
kernel has no truncation beyond that of H. The p-shell space is exactly
decoupled, so every MR-EOM root of a given J must lie at or above the
corresponding non-valence FCI level (Rayleigh-Ritz in a subspace).
`../compare_he6_fci.py` checks this and prints the gaps.

## Files

| file | content |
|---|---|
| `fci_hs_He6.snt` | H after UndoNormalOrdering to the vacuum, FCI model space (zero-body 3.59244929 MeV, not added by KSHELL) |
| `he6_fci_p.ptn` | 2p 4n, positive parity, no truncation (gen_partition option 0) |
| `he6_j0p.input`, `he6_j4p.input` | KSHELL inputs, M=0 and M=2, 30 eigenvalues, J-projected |
| `he6_fci.slurm` | riker batch script (OpenMP, 32 threads) |
| `write_snt.py` | regenerates `fci_hs_He6.snt` from the Hs cache file |

M-scheme dimension: 3.86e6 (M=0), 3.34e6 (M=2); about 6 GB of Lanczos vectors.
On a laptop with 4 threads this ran at about 20 s per Lanczos step.

## Run

The directory is self-contained; only KSHELL is needed.

```bash
rsync -a run_paper/he6_fci/ riker:~/codes/imsrg_sunzhg07/run_paper/he6_fci/
ssh riker 'cd ~/codes/imsrg_sunzhg07/run_paper/he6_fci && sbatch he6_fci.slurm'
# KSHELL_EXE=/path/to/kshell.exe sbatch he6_fci.slurm   to use another build
```

Note: the `kshell_threebody` build stops after reading orbit 15 of the
20-orbit model space; use a standard build (kshell-20230714-ver4 or the one in
`../kshell_he8/`).

## Compare

Copy `log_He6_fci_hs_j0p.txt` and `log_He6_fci_hs_j4p.txt` back here, then

```bash
cd run_paper
PYTHONPATH=../build python3 compare_he6_fci.py
```
