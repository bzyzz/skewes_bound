#!/bin/bash
# Certified zeros of L(s, chi_-8) (8.3) and L(s, chi_12) (12.11) up to 5e4 for the team races mod 8 and mod 12.
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=${PYTHON:-python}
for qn in "8 3" "12 11"; do
  set -- $qn
  echo "=== $1.$2 start $(date)"
  $PY zeros/float_zeros.py $1 $2 50100 16 data/float_$1_$2.npy
  $PY zeros/arb_certify.py $1 $2 data/float_$1_$2.npy data/arb_$1_$2 16
  $PY zeros/turing.py data/arb_$1_$2/certified.npz 40 $(seq -s, 5000 5000 50000)
  echo "=== $1.$2 done $(date)"
done
