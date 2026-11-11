#!/bin/bash
# Consumer-side recheck of every stored zero interval (verify_endpoints.py), for all characters.
# usage: run_endpoints.sh [nproc=16]     log: data/run_endpoints.log
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=${PYTHON:-python}
NP=${1:-16}
for f in data/arb_*/certified.npz; do
  echo "=== $f start $(date)"
  $PY zeros/verify_endpoints.py "$f" "$NP" || echo "=== $f FAILED"
done
echo "=== all done $(date)"
