#!/bin/bash
#SBATCH -A nph123
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --gpus=1
#SBATCH --cpus-per-gpu=32
#SBATCH --mem-per-gpu=380G
#SBATCH -t 06:00:00
#SBATCH -J he8_rho3
#SBATCH -o he8_rho3_full.%j.out
# (1) rho3 norm tests (kernel vs commutator) for He8 itz=0/1 at emax3
# (2) He8 full MR-EOM (8 jobs + He8->Li8 cx1/cx2/cx3), rho3 off/on  [he8_rho3_full.py]
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
PY=/autofs/nccs-svm1_sw/riker/miniforge3/25.3.1-0/bin/python3
module load openblas/0.3.33-omp
export OMP_NUM_THREADS=32
T=rho3_test
REF=../run2/he8_rho3/he8_gs_e3.ref
echo "start $(date)"
for a in "0 0 0" "2 0 0" "1 0 1" "2 0 1" "3 0 1"; do
  tag=$(echo $a | tr ' ' _)
  (cd $T && srun -n1 -c32 $PY test_rho3_norm.py 3 ../$REF He4 p-shell 6 $a > test_he8_e3_$tag.log 2>&1)
  echo "He8 e3 J2 P Tz = $a : $(grep 'max relative' $T/test_he8_e3_$tag.log)"
done
MAX_ITER=${MAX_ITER:-100} PHASE=eom EMAX=3 srun -n1 -c32 $PY -u he8_rho3_full.py
echo "He8 EOM exit $?  $(date)"
echo "done $(date)"
