#!/bin/bash
#SBATCH -A nph123
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --gpus=1
#SBATCH --cpus-per-gpu=32
#SBATCH --mem-per-gpu=380G
#SBATCH -t 06:00:00
#SBATCH -J he6_fci
#SBATCH -o he6_fci.%j.out
#
# Full-space diagonalization of the He6 VS-IMSRG Hamiltonian (EM1.8/2.0,
# hw=16, emax=3, He6 reference, p-shell valence), J=0+ and J=2+.
# Benchmark for the neutral He6 MR-EOM (see README.md).
#
# M-scheme dimension 3.86e6 (M=0), 3.34e6 (M=2); ~6 GB of Lanczos vectors.
#
# Override the KSHELL binary with KSHELL_EXE=/path/to/kshell.exe sbatch he6_fci.slurm
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
module load openblas/0.3.33-omp
KSHELL_EXE=${KSHELL_EXE:-$HOME/codes/kshell-20230714-ver4/bin/kshell.exe}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_GPU:-32}
export GFORTRAN_UNBUFFERED_ALL=y   # progress visible in the logs while running
ulimit -s unlimited
echo "start $(date)  kshell=$KSHELL_EXE  threads=$OMP_NUM_THREADS"
for name in j0p j4p; do
  srun -n1 -c"$OMP_NUM_THREADS" "$KSHELL_EXE" he6_${name}.input > log_He6_fci_hs_${name}.txt 2>&1
  echo "$name exit $?  $(date)"
  grep -E "<H>:" log_He6_fci_hs_${name}.txt | tail -n 30 | head -n 12
done
echo "done $(date)"
