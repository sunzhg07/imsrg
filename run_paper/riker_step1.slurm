#!/bin/bash
#SBATCH -A nph123
#SBATCH -p gpu
#SBATCH --gpus=1
#SBATCH -N 1
#SBATCH --cpus-per-gpu=32
#SBATCH --mem-per-gpu=380G
#SBATCH -t 06:00:00
#SBATCH -J step1_inputs
#SBATCH -o step1_inputs.%j.out
#
# Riker step 1 (compute_report.md §6): regenerate Mg24 (EMAX, default 4) and
# O14/N14 (emax 3) inputs at the current build.
#   cd run_paper && sbatch riker_step1.slurm          # EMAX=4
#   cd run_paper && EMAX=6 sbatch riker_step1.slurm   # Mg24 only matters for 6
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"

module load openblas/0.3.33-omp   # kshell.exe/transit.exe need libopenblas.so.0
PY=/autofs/nccs-svm1_sw/riker/miniforge3/25.3.1-0/bin/python3
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-$(nproc)}}
export EMAX=${EMAX:-4}
KBIN=$HOME/codes/kshell-20230714-ver4/bin

echo "host $(hostname)  EMAX=$EMAX  OMP=$OMP_NUM_THREADS  commit $(git rev-parse --short HEAD)"

# Mg24: VS-IMSRG (Hs, .snt, E1c/E1p), then KSHELL densities.
$PY mg24_vs_imsrg.py
W=kshell_mg24/e$EMAX
mkdir -p "$W"
if [[ ! -f $W/mg24_sd.ptn ]]; then
  snt=$(ls ../hs_cache/Hs_Mg24_sdshell_*_hw16_emax${EMAX}.snt)
  cp "$snt" "$W/mg24.snt"
  (cd "$W" && echo | $PY "$KBIN/gen_partition.py" mg24.snt mg24_sd.ptn 4 4 + > gen_mg24.log)
fi
$PY mg24_kshell_rdm.py

# O14 -> N14 GT (emax 3; independent of EMAX above).
EMAX=3 $PY o14_n14_gt.py

echo "done $(date)"
