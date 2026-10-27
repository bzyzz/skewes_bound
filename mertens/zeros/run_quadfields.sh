#!/bin/bash
# Zeros of L(s, chi_5) and L(s, chi_13) (Q(sqrt 5), Q(sqrt 13)): float location to 20100, Arb certification, Turing.
cd "$(dirname "$0")/.."
Z=../race/zeros   # shared zero tools (float location, Arb certification, Turing)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=${PYTHON:-python}
for qn in "5 4" "13 12"; do
  set -- $qn
  echo "=== $1.$2 float start $(date)"
  $PY $Z/float_zeros.py $1 $2 20100 16 zeros/float_$1_$2.npy
  echo "=== $1.$2 certify start $(date)"
  $PY $Z/arb_certify.py $1 $2 zeros/float_$1_$2.npy zeros/arb_$1_$2 16
  $PY $Z/turing.py zeros/arb_$1_$2/certified.npz 40 $(seq -s, 2000 2000 20000)
  echo "=== $1.$2 done $(date)"
done
