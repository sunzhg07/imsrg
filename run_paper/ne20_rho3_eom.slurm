#!/bin/bash
#SBATCH -A nph123
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --gpus=1
#SBATCH --cpus-per-gpu=32
#SBATCH --mem-per-gpu=380G
#SBATCH -t 06:00:00
#SBATCH -J ne20_rho3
#SBATCH -o ne20_rho3_eom.%j.out
# Ne20 emax4: VS-IMSRG -> KSHELL (1b/2b/3b densities of 0+_1) -> 0+ MR-EOM, rho3 off/on.
#   cd run_paper && sbatch ne20_rho3_eom.slurm
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
PY=/autofs/nccs-svm1_sw/riker/miniforge3/25.3.1-0/bin/python3
module load openblas/0.3.33-omp
export OMP_NUM_THREADS=32
echo "host $(hostname) start $(date)"
export PHASE=${PHASE:-all}; EMAX=${EMAX:-4} srun -n1 -c32 $PY -u ne20_rho3_eom.py
echo "done $(date)"
