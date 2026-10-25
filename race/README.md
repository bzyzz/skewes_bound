# First sign change of π(x;24,11) − π(x;24,1)

**Theorem.** π(x;24,11) < π(x;24,1) for some x with log x ∈ [111.5456338, 111.5460262]. Hence the first sign change of this prime number race occurs before 10^48.4439.

The previous bound was 10³⁵³ (Ford–Hudson, Acta Arith. 100 (2001), Corollary 4). They had located this region (their Table 3, ω ≈ 111.546) but could not certify it. Together with their other results, every race π(x;24,b) − π(x;24,1), with (b,24) = 1 and b ≠ 1, now changes sign below 10^70.021. The limiting race is b = 5.

The proof uses the two-scale kernel of this repository (see [`../THEOREM.md`](../THEOREM.md) §1) inside a Lehman/Ford–Hudson explicit-formula criterion. It needs the zeros up to T₀ = 10⁵ of the four L-functions whose characters satisfy χ(11) = −1. The statement, the inputs and the proof outline are in [`THEOREM.md`](THEOREM.md).

| Character (Conrey) | Conductor | Parity | Zeros 0 < γ ≤ 10⁵ | GRH height used (Platt 2016) |
|---|---|---|---|---|
| 3.2 | 3 | odd | 155,553 | 3.3·10⁷ |
| 4.3 | 4 | odd | 160,131 | 2.5·10⁷ |
| 8.5 | 8 | even | 171,163 | 1.25·10⁷ |
| 24.11 | 24 | even | 188,648 | 4.1·10⁶ |

Certified result (`criterion/cert_T100000_w111.5458300.json`):
- parameters: ω = 111.5458300, η₁ = 1.2·10⁻⁴, δ = 7.62·10⁻⁵, M = 40, κ₆ core;
- U(ω) = −0.0050 ± 5.1·10⁻⁴ < 0.

## Pipeline

1. **`zeros/float_zeros.py`** locates the zeros of Z(t,χ) in double precision (Euler–Maclaurin), with a near-miss rescan for close pairs. It is a locator only.
2. **`zeros/arb_certify.py`** certifies each zero in Arb (20 digits): Z(γ ± 10⁻⁸) must have strictly opposite signs, and the certified intervals must be disjoint.
3. **`zeros/turing.py`** proves the list complete:
   - It certifies L(σ,χ) > 0 on [½, 2], so that N(T) = θ(T)/π + S(T).
   - It applies Turing's method with max(Trudgian 2011, Rumely 1993) for ∫S at T = 5000, 10000, …, 10⁵.
   - The result is that every zero with 0 < |γ| ≤ 10⁵ is simple, lies on the critical line, and is in the list.
4. **`criterion/race_criterion.py`** evaluates the error terms R₃, R₄, R₅, ε₀ and the main term in Arb. R₄ uses Bennett–Martin–O'Bryant–Rechnitzer (2021); R₅ is bounded by integration by parts above Platt's GRH heights.
5. **`criterion/certify.py`** computes the zero sums S_χ in Arb over the certified intervals, checks every hypothesis, and writes the certificate.

`criterion/search.py` is the floating-point parameter search used to choose ω, η₁, δ and M. `heuristics/` holds the exploratory scans that located the region and compared it with Ford–Hudson. They are not rigorous; regenerate their float zeros with `heuristics/zeros_L.py`.

## Reproducing

From the repository root:

```bash
make race-check     # Turing completeness of the shipped zeros, then the Arb certificate (~15 s on 8 cores)
make race-zeros     # recompute the zeros from scratch (~80 core-hours), then Turing
```

Or directly:

```bash
cd race
for c in 3_2 4_3 8_5 24_11; do python zeros/turing.py data/arb_$c/certified.npz 40 $(seq -s, 5000 5000 100000); done
cd criterion && python certify.py 100000 111.5458300 1.2e-4 7.62e-5 40 0,0,1 8
```

`data/arb_<q>_<n>/certified.npz` holds the certified zero intervals [lo, hi] (exact float64 endpoints at which Z was evaluated), plus the parameters `q, n, eps, dps`. `turing_<q>_<n>.json` holds the Turing certificates.

## Status

- An independent line-by-line review of the proof outline, the code and the cited inputs (2026-10-25) found no error affecting the result. The corrections it prompted are recorded in `THEOREM.md`:
  - two misprinted constants in Ford–Hudson's explicit formula (effect below 10⁻²¹);
  - a missing justification for the interchange of the zero sum with the integral, now added.
- Expert review is still pending.
