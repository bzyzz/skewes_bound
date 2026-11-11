#!/usr/bin/env python3
"""Hash manifest of the numerical evidence behind the certificates (zero files and per-batch Arb sums).

The evidence directories are too large for the git repository (about 5 GB). This script fixes their
content: for each directory it hashes every file (SHA-256), writes the full list "sha256  size  name"
to <out>/<label>.sha256.gz, and records in evidence_manifest.json the number of files, the total size
and the SHA-256 of that (uncompressed) list. A copy of the data can then be checked against the
manifest with `make_manifest.py --verify`.

Directories (label: path, what it holds):
  zero_cache_low:  $MERTENS_ZERO_CACHE     Arb zeros n <= 2,325,572 (text, one ball per line)
  zero_cache_ext:  ../zero_cache_ext       Arb zeros 2,325,573 <= n <= 13,955,181
  zero_cache_arb:  $ARB_ZERO_DIR           Arb zeros 13,955,182 <= n <= 195,565,267 (npz batches of 2000)
  arb_batches:     arb_batches             per-batch partial sums of S_K written by arb_recert.py
The LMFDB files used by the Mertens 10^9 certificate are fixed by mertens/lmfdb/md5.txt (LMFDB's own list).

usage: make_manifest.py [--verify] [out_dir=evidence_manifest]
"""
import gzip
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DIRS = {
    "zero_cache_low": os.environ.get("MERTENS_ZERO_CACHE", os.path.join(ROOT, "..", "mertens_bound", "certification", "zero_cache")),
    "zero_cache_ext": os.path.join(ROOT, "zero_cache_ext"),
    "zero_cache_arb": os.environ.get("ARB_ZERO_DIR", os.path.join(ROOT, "zero_cache_arb")),
    "arb_batches": os.path.join(HERE, "arb_batches"),
}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def listing(folder):
    """The text 'sha256  size  name\\n' for every regular file directly in `folder`, sorted by name."""
    names = sorted(n for n in os.listdir(folder) if os.path.isfile(os.path.join(folder, n)))
    lines, total = [], 0
    for n in names:
        p = os.path.join(folder, n)
        size = os.path.getsize(p)
        total += size
        lines.append(f"{sha256_file(p)}  {size}  {n}\n")
    return "".join(lines), len(names), total


def main():
    args = [a for a in sys.argv[1:] if a != "--verify"]
    verify = "--verify" in sys.argv[1:]
    out = args[0] if args else os.path.join(HERE, "evidence_manifest")
    summary_path = os.path.join(HERE, "evidence_manifest.json")
    if verify:
        want = json.load(open(summary_path))
        ok = True
        for label, folder in DIRS.items():
            text, n, total = listing(folder)
            got = dict(files=n, bytes=total, list_sha256=hashlib.sha256(text.encode()).hexdigest())
            same = all(want[label][k] == got[k] for k in got)
            ok &= same
            print(f"{label}: {'OK' if same else 'MISMATCH'} ({n} files, {total} bytes)", flush=True)
        sys.exit(0 if ok else 1)
    os.makedirs(out, exist_ok=True)
    summary = {}
    for label, folder in DIRS.items():
        text, n, total = listing(folder)
        with gzip.open(os.path.join(out, f"{label}.sha256.gz"), "wt") as f:
            f.write(text)
        summary[label] = dict(files=n, bytes=total, list_sha256=hashlib.sha256(text.encode()).hexdigest())
        print(f"{label}: {n} files, {total} bytes, list sha256 {summary[label]['list_sha256']}", flush=True)
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=1)
    print("->", summary_path)


if __name__ == "__main__":
    main()
