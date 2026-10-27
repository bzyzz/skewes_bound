#!/bin/bash
# Extend the certified zeros of L(s, chi_-8) (8.3) to T = 1e5 (team race mod 8, region near 10^28.13).
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=${PYTHON:-python}
$PY zeros/float_zeros.py 8 3 100100 16 data/float_8_3.npy
$PY zeros/arb_certify.py 8 3 data/float_8_3.npy data/arb_8_3 16
$PY zeros/turing.py data/arb_8_3/certified.npz 40 $(seq -s, 5000 5000 100000)
echo "=== 8.3 to 1e5 done $(date)"
