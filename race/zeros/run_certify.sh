#!/bin/bash
# Arb certification + Turing for each character once its float zeros are ready.
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=${PYTHON:-python}
TS=$(seq -s, 5000 5000 100000)
for qn in "3 2" "4 3" "8 5" "24 11"; do
  set -- $qn
  until grep -q "=== $1.$2 done" data/run_float.log; do sleep 60; done
  if grep -q "=== 24.11 done" data/run_float.log; then P=24; else P=6; fi
  echo "=== certify $1.$2 with $P procs, start $(date)"
  $PY zeros/arb_certify.py $1 $2 data/float_$1_$2.npy data/arb_$1_$2 $P
  $PY zeros/turing.py data/arb_$1_$2/certified.npz 40 $TS
  echo "=== certify $1.$2 done $(date)"
done
