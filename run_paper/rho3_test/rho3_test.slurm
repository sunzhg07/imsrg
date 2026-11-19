#!/bin/bash
#SBATCH -A nph123
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --gpus=1
#SBATCH --cpus-per-gpu=32
#SBATCH --mem-per-gpu=380G
#SBATCH -t 01:00:00
#SBATCH -J rho3_test
#SBATCH -o rho3_test.%j.out
# (1) OSM Ne20/USDB 3-body density (cross-check of the kshell TDTD -> .ref conversion)
# (2) emax=4 rho3 norm test (kernel vs commutator), covers ppvv qv<-vv and pphv vv<-vc
set -x
W=/autofs/nccs-svm1_home1/wolf/codes/imsrg_sunzhg07/run_paper/rho3_test
PY=/autofs/nccs-svm1_sw/riker/miniforge3/25.3.1-0/bin/python3
module load openblas/0.3.33-omp
export OMP_NUM_THREADS=32

cd $W/osm
srun -n1 -c32 ./sm.exe sm.dat > osm.log 2>&1
echo "OSM exit $?"

cd $W
srun -n1 -c32 $PY test_rho3_norm.py 4 ne20/ne20_kshell.ref O16 sd-shell 6 > test_emax4_kshell.log 2>&1
echo "emax4 kshell exit $?"
if [ -f osm/ne20_osm.ref ]; then
  srun -n1 -c32 $PY test_rho3_norm.py 4 osm/ne20_osm.ref O16 sd-shell 6 > test_emax4_osm.log 2>&1
  echo "emax4 osm exit $?"
fi
