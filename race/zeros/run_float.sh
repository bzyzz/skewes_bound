#!/bin/bash
# Floating-point zero location for the four characters of the race 24;11 (chi(11) = -1).
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=${PYTHON:-python}
for qn in "3 2" "4 3" "8 5" "24 11"; do
  set -- $qn
  echo "=== $1.$2 start $(date)"
  $PY zeros/float_zeros.py $1 $2 100100 20 data/float_$1_$2.npy
  echo "=== $1.$2 done $(date)"
done
