# Mertens' first theorem: an explicit first sign change

Let $\Delta_1(x) = \sum_{p\le x} (\log p)/p - \log x - E$, where $E = -\gamma - \sum_p \log p/(p(p-1)) = -1.33258\ldots$.

- Rosser and Schoenfeld showed that $\Delta_1(x) > 0$ for $x \le 10^8$.
- Lay (2015) proved that $\Delta_1$ changes sign infinitely often.
- No explicit bound for the first sign change was known. This is Problem 14 of the comparative prime number theory problem list (PIMS 2024).

This directory certifies such a bound with the two-scale kernel of the main repository. The statements, the derivation and the reviews are in [`THEOREM.md`](THEOREM.md).

## Results

| Result | Statement | Zeros used | Certificate |
|---|---|---|---|
| Main | $\Delta_1(x) < 0$ for some $x$ with $\log x \in [495.7028331177, 495.7028331467]$, so for some $x < 10^{215.281005099}$ | all $2{,}846{,}548{,}032$ zeros below $10^9$: ours (Arb) below $8\cdot10^7$, D. J. Platt's (LMFDB) above | `lmfdb/first_cert_lmfdb_T1000000000.json` |
| Self-contained | $\Delta_1(x) < 0$ for some $x$ with $\log x \in [1859.129042158, 1859.129060642]$, so for some $x < 10^{807.4095}$ | our $2{,}325{,}572$ Arb zeros below $1.3\cdot10^6$ | `criterion/first_cert_T1300000.json` |
| Quadratic fields (Problem 8) | the first integer $x$ with $\Delta_K(x) < 0$ is $28$ for $\mathbb{Q}(\sqrt5)$, $16$ for $\mathbb{Q}(\sqrt{13})$, and at most $148$ for every field with $\lvert D\rvert \le 300$ | none (direct evaluation in Arb) | `scan/fields_first_negative.json` |

**For comparison:**
- Büthe (Acta Arith. 2015) proved that the first sign changes in Mertens' second and third theorems occur before $e^{495.702833165}$.
- The main result lies in the same region and is slightly lower. The first theorem is harder there: the same zero sum must exceed $1$ rather than $1 - 2/\log x$, and with the zeros below $10^9$ it does so only on an interval of width $2.3\cdot10^{-8}$ in $\log x$.

## Reproducing

Run from the repository root (see the main [`README.md`](../README.md) for the Python environment):

```bash
make mertens-build    # libplatt.so (LMFDB decoder) and libscan.so (heuristic scans)
make mertens-fields   # quadratic fields, certified in Arb (< 1 s)
make mertens-certify  # the 10^807.41 certificate (~2 min on 16 cores; needs the zero caches, see ../README.md, Data)
```

The main certificate also needs Platt's zeros between $8\cdot10^7$ and $10^9$ (33 GB, 440 files) and our Arb zeros up to $8\cdot10^7$ (`zero_cache_arb/`):

```bash
make mertens-download            # LMFDB files listed in lmfdb/files_to_1e9.txt into lmfdb/data/, md5-checked
make mertens-lmfdb WORKERS=24    # the 10^215.281 certificate (~70 min on 24 cores)
```

Set `LMFDB_ZETA_DIR` to use LMFDB files stored elsewhere.

The checks made on Platt's data:
- the md5 sum of every file;
- index contiguity across all sources;
- $N(8\cdot10^7)$ and $N(10^9)$, recomputed independently with Arb's `arb_zeta_nzeros`;
- the seam at zero $\#195{,}565{,}267$.

An independent review also checked the decoder (`lmfdb/platt.c`) bit for bit against a separate pure-Python decoder on four files, and 826 zeros against Arb's `acb_dirichlet_zeta_zeros`.

**Provenance.** `lmfdb/first_cert_lmfdb_T1000000000.json` records the SHA-256 of `certify_first.py` and `certify_first_lmfdb.py` as they were before this directory was moved into the repository. The move changed only path constants and a docstring; the other recorded hashes (`certio.py`, `error_bounds.py`, `platt.c`, `platt.py`, `libplatt.so` built with `-O2`, `kernel_fast.py`) are current. `criterion/first_cert_T1300000.json` was regenerated after the move, with identical results.

## Layout

```text
THEOREM.md                     statements, inputs, derivation, certified numbers, reviews
criterion/certify_first.py     Arb criterion for Delta_1 with our zeros (uses ../certification/{certio,error_bounds}.py)
criterion/search_first.py      floating-point parameter search
lmfdb/certify_first_lmfdb.py   the 10^9 certificate: our Arb zeros + Platt's LMFDB zeros
lmfdb/platt.c, platt.py        exact decoder for Platt's zero files
lmfdb/kernel_fast.py           rigorous series for the mollifier transform
lmfdb/download.sh, md5.txt, files_to_1e9.txt   download and verification of the LMFDB files
lmfdb/float_pass.py, search_215.py             floating-point pass and parameter search near 10^215.28
scan/                          heuristic scans (wide.py, refine.py), validation against primes, quadratic fields
sim/sim215.py                  pre-computation simulation: how many zeros a certificate near 10^215.28 needs
zeros/                         zeros of L(s, chi_5), L(s, chi_13) to 2e4 (heuristics for the quadratic fields)
```
