# Gaussian baseline for the θ(x) comparison (non-rigorous)

These scripts produced the Gaussian column of the paper's θ(x) comparison table. They run a floating-point search with the method of Platt–Trudgian (2016, Theorem 3 and §3.1.6), on the same zeros as the two-scale certificates. The output is the best endpoint found in the search domain. It is neither a certificate nor a proven optimum.

- `explore.py` loads the float zeros and evaluates the zero sum on an ω-grid. It uses binned Taylor expansions in long double.
- `explore_gauss.py` runs the search. For each `T0` it scans:
  - `T0/sqrt(alpha)` in {3, 3.25, …, 4.5};
  - `A` ≤ 3·10¹², taking the smallest value (in 1% steps) for which the error term R4 of Theorem 3 is below 1e-10, with `eta = 2A/alpha`;
  - `omega` on a grid of step 2e-10 over 727.951334 + [−2.2, 2.0]·10⁻⁶.

  It then applies the §3.1.6 trimming and reports the smallest `omega + eta0` with a positive margin.
- `gauss3_*.log` are the runs reported in the paper: `python explore_gauss.py <T0>`.

**Data.** `explore.py` reads float zeros (long-double ordinates to 8·10⁷) through the absolute path `S` near its top. These files (3.2 GB) are not part of the repository. To rerun, point `S` to float zeros in the same format: `indep/zeros.npz` with integer and fractional parts `gi`, `gf`, and `indep2/g_hi_ld.npy`.
