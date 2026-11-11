#!/bin/bash
# Download listed LMFDB zero files; any exhausted retry or dispatcher failure is fatal.
set -euo pipefail
cd "$(dirname "$0")"
if [ "$#" -ne 1 ] || [ ! -f "$1" ]; then
  echo "Usage: download.sh EXISTING_FILE_LIST" >&2
  exit 2
fi
mkdir -p data
xargs -P 4 -I{} sh -c '
  set -eu
  f=$1
  case "$f" in ""|*[!A-Za-z0-9._-]*) echo "Invalid filename: $f" >&2; exit 2;; esac
  checksum() {
    if command -v md5sum >/dev/null 2>&1; then
      md5sum "$1" | cut -d" " -f1
    else
      md5 -q "$1"
    fi
  }
  matches() {
    [ -f "$1" ] || return 1
    actual=$(checksum "$1") || return 1
    grep -Fqx "$actual *$f" md5.txt
  }
  if matches "data/$f"; then exit 0; fi
  for try in 1 2 3; do
    if curl -fsSL -b "human=1" --max-time 900 \
        "https://beta.lmfdb.org/data/riemann-zeta-zeros/$f" -o "data/$f.part" \
        && matches "data/$f.part"; then
      mv "data/$f.part" "data/$f"
      echo "ok $f"
      exit 0
    fi
    echo "Download/checksum failed: $f (try $try)" >&2
  done
  echo "FAILED $f" >&2
  exit 1
' sh '{}' < "$1"
echo "DOWNLOAD DONE"
