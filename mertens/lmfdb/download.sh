#!/bin/bash
# Download LMFDB / Platt zeta-zero files listed in $1 into data/, verify md5, 4 in parallel.
cd "$(dirname "$0")"
mkdir -p data
LIST=$1
cat "$LIST" | xargs -P 4 -I{} sh -c '
  f={}
  if [ -f data/$f ] && grep -q "$(md5sum data/$f | cut -d" " -f1) \*$f" md5.txt; then exit 0; fi
  for try in 1 2 3; do
    curl -sL -b "human=1" --max-time 900 https://beta.lmfdb.org/data/riemann-zeta-zeros/$f -o data/$f.part && mv data/$f.part data/$f
    if grep -q "$(md5sum data/$f | cut -d" " -f1) \*$f" md5.txt; then echo "ok $f"; exit 0; fi
    echo "md5 mismatch $f (try $try)"
  done
  echo "FAILED $f"'
echo "DOWNLOAD DONE"
