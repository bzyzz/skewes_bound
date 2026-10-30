#!/bin/bash
# Heuristic team-race scans (team_scan.py) for the distinct team races with q <= 30 that are open at 1e9.
cd "$(dirname "$0")"
PY=${PYTHON:-python}; export OMP_NUM_THREADS=${NT:-16}
for job in "5 3000" "15 3000" "20 3000" "21 3000" "28 3000" "12 1000" "24 5000"; do
  set -- $job
  echo "=== q=$1 UMAX=$2 $(date +%T)"
  $PY team_scan.py $1 $2 2000 10000 0.3 > scans/team_q$1.out 2>&1
  $PY - "$1" "$2" <<'P'
import json, sys
r = json.load(open(f"scans/team_q{sys.argv[1]}_U{sys.argv[2]}.json"))
print("  chars", r["chars"], "candidates", r["n_candidates"], "best_jp", r["best_jp"])
print("  first", r["first"])
P
done
echo "SCANS DONE $(date)"
