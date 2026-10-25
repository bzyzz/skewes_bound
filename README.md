# A two-scale compactly supported kernel for Lehman's method

Certified computations, in Arb ball arithmetic, for:
- upper bounds on **Skewes's number** $\Xi := \min\{x \ge 2 : \pi(x) > \mathrm{li}(x)\}$;
- the first sign change of $\theta(x) - x$;
- the **prime number race** $\pi(x;24,11) - \pi(x;24,1)$.

All three use a new two-scale kernel in Lehman's explicit-formula method.

The kernel is $K = K_{\mathrm{core}} * \psi_\delta$:
- A low-order cosine-power core of half-width $\eta_1$ sets the averaging window.
- A much narrower $\cos^{2M}$ mollifier makes negligible the contribution of any zeros off the critical line above the verified height $A$ ($3\cdot 10^{12}$ for $\zeta$, by Platt–Trudgian 2021).

This decouples the window half-width $\eta_0$ from the zero cutoff $T_0$. We work at $\eta_0 T_0 \approx 12$; the truncated Gaussian kernels of all earlier Skewes computations needed $\eta_0 T_0 \approx 114$–$136$. Statements, inputs and proof outlines are in [`THEOREM.md`](THEOREM.md) and [`race/THEOREM.md`](race/THEOREM.md). A manuscript is in preparation.

## Results

| Result | Statement | Zeros used | Status |
|---|---|---|---|
| Skewes, main theorem | $\pi(x) > \mathrm{li}(x)$ for some $x \in [e^{727.951335292}, e^{727.951335596}]$, so $\Xi < e^{727.951335596}$ | 195,565,267 ($T_0 = 8\cdot 10^7$) | certified in Arb |
| $\theta(x) - x$ | $\theta(x) > x$ for some $x \in [e^{727.9513323698}, e^{727.9513328538}]$ | 118,488,123 ($T_0 = 5\cdot 10^7$) | certified in Arb |
| Prime race mod 24 | $\pi(x;24,11) < \pi(x;24,1)$ for some $x$ with $\log x \in [111.5456338, 111.5460262]$, so the first sign change is $< 10^{48.4439}$ | 675,495 zeros of four Dirichlet $L$-functions ($T_0 = 10^5$) | certified in Arb |
| Smith (2016), RH removed | $\Xi < e^{727.951335402030}$ | Platt's sum over $1.04\cdot 10^{11}$ zeros | **conditional** on one tabulated value in Smith's MSc thesis |

**For comparison:**
- **Skewes:** the best published bounds are $\Xi < e^{727.951336105}$ unconditionally and $e^{727.951335621}$ under RH (Saouter–Trudgian–Demichel 2015, $5.25\cdot 10^8$ zeros).
- **$\theta(x) - x$:** Platt–Trudgian (2016) give $\theta(x) > x$ for some $x < e^{727.951332668}$, using about $2\cdot 10^{10}$ zeros.
- **Prime race:** for $\pi(x;24,11) < \pi(x;24,1)$ the previous bound was $10^{353}$ (Ford–Hudson 2001).

**Main theorem.**
- Every zero below $T_0 = 8\cdot10^7$ is isolated with Arb:
  - $n \le 13{,}955{,}181$ from the Arb zero caches;
  - the rest re-isolated with `acb_dirichlet_zeta_zeros` by `certification/arb_recert.py`, about 15 hours on 64 cores.
- The Arb-only certificates are `certification/s_cert*_arb.json`. With them the envelope margin is $I(\omega,\eta_0) \ge 2.32\cdot 10^{-5}$, giving more than $2.75\cdot 10^{153}$ consecutive integers with $\pi(x) > \mathrm{li}(x)$.
- The certificates without the `_arb` suffix used floating-point Riemann–Siegel brackets (`fast_cert_zeros.hpp`) above $6.8\cdot 10^6$. They are kept for comparison and give the same bounds with slightly smaller margins.

**Prime race.** See [`race/README.md`](race/README.md).
- The four $L$-functions have conductors 3, 4, 8 and 24.
- Their zeros are certified in Arb and proved complete up to $T_0 = 10^5$ by Turing's method.
- The certified zero intervals are shipped, so the whole certificate can be re-checked from this repository.

## Reproducing

Requirements:
- Python 3.11, [python-flint](https://github.com/flintlib/python-flint) ≥ 0.9.0 (tested with 0.9.0 / FLINT 3.6.0), and numpy.
- A C compiler with OpenMP, only for the floating-point zero location and the heuristic scans.

```bash
conda env create -f environment.yml && conda activate skewes   # or: pip install -r requirements.txt
make test          # regression tests for outward rounding and failure modes
make check-only    # re-check all shipped certificates (Skewes float and Arb, Smith) against the Arb error bounds
make race-check    # race: Turing completeness of the shipped zeros, then the Arb certificate (~15 s on 8 cores)
```

Full recomputation needs the zero data below:

```bash
make arb-recert WORKERS=64   # isolate all zeros up to 8e7 in Arb and sum S_K per batch (~15 h)
make arb-combine             # combine into certification/*_arb.json (checks coverage and N(T0))
make race-zeros              # locate and certify the race zeros, then run Turing (~80 core-hours)
```

The θ certificate and the conditional Smith refinement are checked with

```bash
cd certification
python theta_certify.py 727.9513329060 1.618e-6 2e-9 60 14,72,14 6800000
python theta_certify.py 727.9513326118 2.4e-7 2e-9 60 4,52,44 50000000    # uses zero_cache_arb/
python certify_smith_unconditional.py
```

### Provenance

Each certificate records the SHA-256 of the code that produced it (`code_sha256`). The recorded hashes are not all current:
- For `certify_S.py` and `theta_certify.py`, they predate the edit that replaced absolute paths by repository-relative ones. That edit changed path constants only.
- Three early certificates (`s_cert_T3e6`, `s_cert_T4p3e6`, `s_cert_T6p8e6_st`) predate the last revision of `error_bounds.py`.

`check_theorem.py` always recomputes the error terms with the current `error_bounds.py`, and every shipped certificate passes with the current code.

### Data

The zeta zero data are too large for git and are regenerated by the listed scripts; no download is provided. The shipped certificates can be re-checked without them (`make check-only`).

| Data | Size | Produced by |
|---|---|---|
| Arb zero balls, $n \le 2{,}325{,}572$ | 0.6 GB | [mertens_bound](https://github.com/bzyzz/mertens_bound) (`certification/zero_cache/`); point `$MERTENS_ZERO_CACHE` at it (default `../mertens_bound/certification/zero_cache`) |
| Arb zero balls, $2{,}325{,}573 \le n \le 13{,}955{,}181$ | 0.4 GB | `build_zero_cache.py` → `zero_cache_ext/` |
| Arb zero balls, $13{,}955{,}182 \le n \le 195{,}565{,}267$ | 4.2 GB | `certification/arb_recert.py` → `zero_cache_arb/` (per-batch sums in `certification/arb_batches/`) |
| floating-point zero brackets above $6.8\cdot 10^6$ | 2.2 GB | `build_and_scan_50m.py`, `extend_to_80m.py` → `chunks_bin/` (cross-check only) |
| Race: certified zero intervals, $0 < \gamma \le 10^5$ | 11 MB | `race/zeros/` pipeline → `race/data/arb_*/certified.npz` (**shipped**) |

## Layout

```text
THEOREM.md                             Skewes / theta / Smith: statements, inputs, certified numbers
certification/
  certio.py                            outward-safe Arb serialisation, provenance, checks
  error_bounds.py                      hat K and the error terms R_2..R_5 in Arb
  certify_S.py                         S_K(omega, T0) over all zeros below T0
  arb_recert.py, arb_combine.py        Arb isolation of all zeros up to 8e7 and the Arb-only certificates
  check_theorem.py                     combines a certificate with the error bounds
  theta_certify.py                     theta(x) > x criterion, all in Arb
  certify_smith_unconditional.py       conditional refinement of Smith (2016)
  test_certio.py                       regression tests
  s_cert*_arb.json                     Arb-only certificates for pi(x) - li(x)
  s_cert*.json, s_certificate.json     certificates with floating-point zeros above 6.8e6 (comparison)
  theta_cert_*.json                    certificates for theta(x) - x
  smith_unconditional_certificate.json certificate for the conditional Smith refinement
race/                                  pi(x;24,11) - pi(x;24,1): zeros, Turing, criterion, certificate (see race/README.md)
build_zero_cache.py                    Arb zero cache builder
fast_cert_zeros.hpp, fast_cert_lib.cpp floating-point Riemann-Siegel zero bracketing (cross-check)
build_and_scan_50m.py, extend_to_80m.py driver for the floating-point zeros above 6.8e6
nzeros_edges_*.json                    Arb zero counts N(T) at the 2000-unit chunk edges
global_scan_nufft.py, fast_skewes.cpp  heuristic NUFFT scan for earlier candidate regions
optimize_params.py, extract_zeros.py   parameter search and zero extraction helpers
```

## License

MIT; see [LICENSE](LICENSE).
