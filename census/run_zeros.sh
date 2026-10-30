#!/bin/bash
# Float zeros (heuristic) for every primitive character of conductor <= DMAX, up to height T.
cd "$(dirname "$0")"
PY=${PYTHON:-python}; T=${T:-5000}; NP=${NP:-20}; DMAX=${DMAX:-30}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
$PY list_chars.py $DMAX > chars_$DMAX.txt
while read d n par real; do
  [ "${REAL_ONLY:-1}" = 1 ] && [ "$real" = 0 ] && continue
  out=zeros/z_${d}_${n}_T${T}.npy
  [ -f "$out" ] && continue
  echo "=== $d.$n $(date +%T)"
  $PY float_zeros_any.py $d $n $T $NP $out | head -1
done < chars_$DMAX.txt
echo "ALL DONE $(date)"
