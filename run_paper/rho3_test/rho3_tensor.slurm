#!/bin/bash
#SBATCH -A nph123
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --gpus=1
#SBATCH --cpus-per-gpu=32
#SBATCH --mem-per-gpu=380G
#SBATCH -t 02:00:00
#SBATCH -J rho3_tensor
#SBATCH -o rho3_tensor.%j.out
# emax=4 rho3 norm test (kernel matmul vs commutator), scalar and tensor,
# Ne20 USDB kshell_threebody rho3; covers ppvv qv<-vv and pphv vv<-cv.
W=/autofs/nccs-svm1_home1/wolf/codes/imsrg_sunzhg07/run_paper/rho3_test
PY=/autofs/nccs-svm1_sw/riker/miniforge3/25.3.1-0/bin/python3
module load openblas/0.3.33-omp
export OMP_NUM_THREADS=32
cd $W
for args in "0 0 0" "2 0 0" "1 1 0" "3 1 0" "4 0 0"; do
  tag=$(echo $args | tr ' ' _)
  srun -n1 -c32 $PY test_rho3_norm.py 4 ne20/ne20_kshell.ref O16 sd-shell 6 $args > test_emax4_$tag.log 2>&1
  echo "J2 P Tz = $args : $(grep -E '^max' test_emax4_$tag.log)"
done
