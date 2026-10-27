# Sign changes in Mertens' first theorem below 10^215.2811 and 10^807.41 — statements, inputs, proof outline

Status: **draft; independent review done (2026-10-25), no error affecting the result.** Certificates: `criterion/first_cert_T1300000.json` (Theorem 1, §1–§5) and `lmfdb/first_cert_lmfdb_T1000000000.json` (Theorem 2, §8).

Notation:
- Δ₁(x) := Σ_{p≤x} (log p)/p − log x − E, with E = −γ − Σ_p log p/(p(p−1)) = −1.33258227573322…
- θ(x) = Σ_{p≤x} log p, ψ(x) = Σ_{p^m≤x} log p, and e(t) := θ(t) − t.
- ρ = β + iγ runs over the nontrivial zeros of ζ(s), counted with multiplicity.
- A = 3·10¹².

Background: Rosser and Schoenfeld found Δ₁(x) > 0 for 0 < x ≤ 10⁸. Lay (2015) proved that Δ₁ changes sign infinitely often. No explicit bound for the first sign change was known; this is Problem #14 of the CPNT problem list (PIMS 2024).

## 1. Statement

**Theorem.** Δ₁(x) < 0 for some x with log x ∈ [1859.129042158, 1859.129060642]. Hence the first sign change of Δ₁ occurs before e^{1859.129060642} < 10^{807.4095}.

## 2. Inputs

- **[cited] Explicit formula.** ψ₀(x) = x − Σ_ρ x^ρ/ρ − log 2π − ½ log(1 − x⁻²) for x > 1, the sum taken symmetrically in |γ|.
- **[derived] Integrated explicit formula.** ∫_x^∞ (ψ(t) − t)t⁻² dt = Σ_ρ x^{ρ−1}/(ρ(ρ−1)) − ∫_x^∞ c(t)t⁻² dt, where c(t) = log 2π + ½ log(1 − t⁻²). The sum converges absolutely.
  - Integrate the truncated explicit formula over [x, X] first. The truncation error O(t log²(tT)/T + log t·min(1, t/(T⟨t⟩))) integrated against t⁻² over the finite range tends to 0 as T → ∞.
  - Then let X → ∞: each term X^{ρ−1}/(ρ(ρ−1)) tends to 0 (β < 1), dominated by 1/|ρ(ρ−1)|, which is summable.
  - A uniform truncation bound cannot be integrated over [x, ∞) directly, because the error integrated against t⁻² diverges.
- **[cited] Platt–Trudgian, Math. Comp. 85 (2016), (9) (arXiv:1407.1914, checked in the original).** |θ(y) − y| ≤ 1.5423·10⁻⁹·y for y ≥ e²⁰⁰. This is the same input as in the θ(x) − x criterion (`../certification/theta_certify.py`).
  - (9) rests on Faber–Kadiri's tables, and Faber–Kadiri later published a corrigendum.
  - The argument only needs some ε ≲ 3·10⁻³, so any published bound of this form suffices.
- **[cited] Rosser–Schoenfeld (1962), Theorem 9.** θ(y) < 1.01624·y for y > 0.
- **[cited] Platt–Trudgian, Bull. LMS 53 (2021).** All zeros with 0 < |γ| ≤ A are simple and lie on the critical line.
- **[cited] Lehman's zero-density lemmas**, through `error_bounds.bound_R4_*`: the same R₄ as for Skewes (`../THEOREM.md` §2).
- **[computed] Zeros.** The 2,325,572 zeros with 0 < γ ≤ T₀ = 1.3·10⁶, as Arb balls from the `mertens_bound` zero cache (`$MERTENS_ZERO_CACHE`; see `../README.md`, Data).
  - `arb_zeta_nzeros` certifies N(T₀) = 2,325,572.
  - γ_N < T₀ < γ_{N+1} is checked with `acb_dirichlet_zeta_zeros`.

## 3. Exact representation [derived]

From Σ_{p≤x} (log p)/p = ∫_{2⁻}^x dθ(t)/t = θ(x)/x + ∫_2^x θ(t)t⁻² dt and E = 1 − log 2 + ∫_2^∞ e(t)t⁻² dt (both sides tend to 0 by the PNT), for every x ≥ 2:

Δ₁(x) = e(x)/x − ∫_x^∞ e(t) t⁻² dt.

Now write θ = ψ − θ(√·) − R_m with R_m(t) := Σ_{m≥3} θ(t^{1/m}), so that 0 ≤ R_m(t) ≤ 1.01624(t^{1/3} + log₂(t)·t^{1/4}). Put x = eᵘ and F(u) := √x·Δ₁(x). Then

F(u) = −Σ_ρ x^{ρ−½}/(ρ−1) − θ(√x)/√x + √x∫_x^∞ θ(√t)t⁻² dt + (remainder terms),

where the zero terms of e(x)/√x and of the integral combine through 1/ρ + 1/(ρ(ρ−1)) = 1/(ρ−1). The remainder terms come from c(t) and R_m.

Upper bounds, valid for √x ≥ e²⁰⁰:
- **Bias:** with ε = 1.5423·10⁻⁹, −θ(√x)/√x + √x∫θ(√t)t⁻² ≤ −(1−ε) + 2(1+ε) = 1 + 3ε.
- **Remainder:**
  - −c(x)/√x ≤ x^{−5/2};
  - √x∫c(t)t⁻² ≤ log(2π)·x^{−1/2};
  - −R_m(x)/√x ≤ 0;
  - √x∫R_m(t)t⁻² ≤ 1.53x^{−1/6} + 1.36(log₂ x + 2)x^{−1/4}.
  
  Their sum is R_small.

**Validation.**
- At ω = 16.5, 17.3 and 18.1, the kernel-smoothed F(u) computed directly from the primes agrees with 1 − 2 Re Σ_{γ>0} K̂(γ) x^{iγ}/(ρ−1) up to an offset of +0.024 to +0.032 (`scan/validate_first.py`). The offset is mostly the cube-of-primes term ≈ 0.5·x^{−1/6}, plus the fluctuation of θ(√x) − √x.
- The independent review made an exact check. It computed the non-zero part B(u) (θ(√x), all R_m terms, the c(t) terms, and the tail integrals) from primes up to 3.2·10⁸. Then smoothed F = B + Z, with Z the zero sum over γ ≤ 1.3·10⁶, holds to within 2·10⁻⁸ at ω = 14.3, 16.5, 17.3, 18.1, 19.0 and 19.4.
- This pins down the sign, the bias, the weight 1/(ρ−1) and the phase convention e^{+iγω}.

## 4. The criterion [derived]

K = K_core ∗ ψ_δ is the two-scale kernel (`../THEOREM.md` §1): K ≥ 0, supp K = [−η₀, η₀], ∫K = 1. Put I(ω) := ∫K(u−ω)F(u)du. Then I(ω) ≤ U, where

U = 1 + 3ε + S_M(ω,T₀) + R₄ + R₅ + R_small,  S_M = −2 Σ_{0<γ≤T₀} K̂(γ) Re(e^{iγω}/(ρ−1)).

If U < 0, then F(u) < 0, i.e. Δ₁(eᵘ) < 0, for some u ∈ [ω−η₀, ω+η₀].

- **Interchange of Σ_ρ and ∫K.** The truncated explicit formula for ψ (Davenport, Ch. 17) has an error bounded uniformly on the window. So the symmetric partial sums are bounded and converge off the prime powers, and dominated convergence applies. The integral term already converges absolutely.
- **Zeros with 0 < |γ| ≤ A (on the line).** ∫K(u−ω)e^{iγu}du = K̂(γ)e^{iγω} exactly, so the zeros with |γ| ≤ T₀ give S_M. Since |ρ−1| ≥ |γ|, the zeros with T₀ < |γ| ≤ A contribute at most 2Σ_{T₀<γ≤A} |κ̂_core(η₁γ)|/γ ≤ R₄.
- **Zeros with |γ| > A (possibly off the line).** Integrate ∫K(u−ω)e^{(ρ−½)u}du by parts 2M times; there are no boundary terms. This gives modulus ≤ e^{u_hi/2}·2C_M(Mπ/δ)^{2M}|γ|^{−2M}. Dividing by |ρ−1| ≥ |γ| and summing with Σ_{γ>A} γ^{−2M−1} < A^{−2M} log A gives R₅ = 4C_M log A·e^{u_hi/2}(Mπ/(Aδ))^{2M}. This is the R₅^θ of the θ criterion.

Implemented in `criterion/certify_first.py`, adapted from `../certification/theta_certify.py`.

## 5. Certified numbers [computed]

Parameters:
- T₀ = 1.3·10⁶ (N = 2,325,572), ω = 1859.1290514;
- η₁ = 9.24·10⁻⁶ (η₁T₀ = 12.0), δ = 2·10⁻⁹, M = 494;
- λ = (0.04, 0.52, 0.44).

| Term | Value |
|---|---|
| S_M(ω, T₀) | −1.01044645655603 ± 8.2·10⁻¹⁵ |
| 1 + 3ε | 1.0000000046 |
| R₄ (min of envelope 9.84·10⁻⁵ and Stieltjes 7.50·10⁻⁵) | 7.495·10⁻⁵ |
| R₅ | 6.8·10⁻¹⁷⁴ |
| R_small | 4.1·10⁻¹³⁵ |
| **U** | **−0.0103715 < 0** |

Window: log x ∈ [1859.129042158, 1859.129060642], so log₁₀ x < 807.409492183.

## 6. Heuristics

- The true first sign change is probably near u ≈ 495.7028 (x ≈ 10^{215.28}). That is the region where Büthe (2015) certified the first sign change for Mertens' second and third theorems.
- The full-weight sum there is: +0.00267 with zeros to 1.3·10⁶, +0.00155 with zeros to 6.8·10⁶, and −0.00018 with all 195.6M Arb zeros to 8·10⁷. The negative part is only about 3·10⁻⁸ wide.
- The first theorem is harder than the second, because the second theorem's threshold is lower by 2/u (≈ 0.004 at u = 495.7).
- No other region with u < 1859 is negative in the wide scan (`scan/wide.py`, `scan/refine.py`, zeros to 1.3·10⁶).
- **What a certificate at 10^{215.28} would need** (`sim/sim215.py`). Keep the 195.6M zeros below 8·10⁷ exact, and model the higher zeros as a random-phase Gaussian process (scale 1.0 or 1.3, 300 realizations each). The probability that the criterion succeeds is:

  | T₀ | 8·10⁷ | 2·10⁸ | 5·10⁸ | 10⁹ | ≥ 3·10⁹ |
  |---|---|---|---|---|---|
  | Probability | 0 | 0–0.01 | 0.90–0.94 | 0.99–1.00 | 1.00 |

  - The certified endpoint would be log x ≈ 495.7028331, i.e. 10^{215.2810}.
  - The true sum is negative there in every realization.
  - Zeros up to 10⁹ (2.9·10⁹ zeros) are available from LMFDB, which distributes Platt's rigorously isolated zeros up to 3.06·10¹⁰: 35 GB up to 10⁹. Summing them in Arb takes roughly 16 core-hours.

## 7. Review and open items

**Independent review (2026-10-25)** by a fresh agent, which did not see these conclusions:
- re-derived §3–§4;
- reproduced the Arb certificate exactly (S_M and U identical to the JSON);
- recomputed S_M in long double with its own kernel code (agreement within 1.3·10⁻¹⁰);
- checked R₄ empirically: the actual tail over (1.3·10⁶, 6.8·10⁶] is 7.11·10⁻⁵ plus about 2·10⁻⁶ beyond, so the bound 7.495·10⁻⁵ holds with about 2% slack;
- spot-checked the zero cache against `acb_dirichlet_zeta_zeros`.

It found no error affecting U < 0. All analytic error terms together are about 7.5·10⁻⁵, under 1% of the margin. Corrections made here:
- the interval in §1 is now rounded outward;
- the integrated explicit formula now has its justification (§2);
- the validation wording is corrected (§3).

Open items:
- **Code hardening (applied 2026-10-25):** `certify_first.py` now requires every summed zero ball to lie strictly below T₀, and the number of zeros summed to equal N(T₀) exactly. The certificate was regenerated and is unchanged.
- Expert review.

## 8. Sharper result with Platt's zeros (2026-10-25)

**Theorem 2.** Δ₁(x) < 0 for some x with log x ∈ [495.7028331177, 495.7028331467]. Hence the first sign change of Δ₁ occurs before 10^{215.281005099}. This is below Büthe's 2·10²¹⁵ = 10^{215.301} for the second and third theorems.

**Inputs.**
- Our 195,565,267 Arb zeros with γ < 8·10⁷.
- D. J. Platt's zeros (LMFDB, `lmfdb/`) for 195,565,267 < n ≤ N(10⁹) = 2,846,548,032. Each zero is the exact ball t₀ + Z·2⁻¹⁰¹ ± 2⁻¹⁰².
  - Platt's isolation proves that these zeros lie on the critical line, and his Turing check proves the list complete (Platt, *Isolating some non-trivial zeros of zeta*, Math. Comp. 86 (2017); LMFDB).
  - This is an additional trust input that Theorem 1 does not need.

**Checks.**
- The md5 sum of every file was verified.
- On 1,603,659 overlapping zeros near 10⁶, LMFDB agrees with our Arb zeros to 2.3·10⁻¹⁰ (float comparison).
- On 20,000 random overlapping zeros just below 8·10⁷, the LMFDB balls overlap our Arb balls (Arb comparison).
- Indices are contiguous from 195,565,268 to 2,846,548,032. The first zero above T₀ is γ = 1000000000.1156508900…

**Parameters.** T₀ = 10⁹, ω = 495.7028331322, η₁ = 1.4·10⁻⁸ (η₁T₀ = 14), δ = 5·10⁻¹⁰, M = 135, λ = (0.04, 0.52, 0.44). The mollifier factor uses a rigorous series in place of the 135-term product (`lmfdb/kernel_fast.py`).

**Certified numbers** (`lmfdb/first_cert_lmfdb_T1000000000.json`, 3,973 s on 24 workers):

| Term | Value |
|---|---|
| S_M | −1.00024 ± 5.7·10⁻⁶ |
| R₄ (Stieltjes) | 4.787·10⁻⁵ |
| R₅ | 3.9·10⁻³⁸ |
| R_small | 2.0·10⁻³⁶ |
| **U** | **≤ −1.956·10⁻⁴** |

The margin is thin: about 4 times R₄, and 44 times the radius of U. It matches the floating-point prediction (min U = −1.97·10⁻⁴). Code: `lmfdb/certify_first_lmfdb.py`, `platt.c`, `platt.py`, `kernel_fast.py`. **Independent review (2026-10-25)**, by a fresh agent without access to these conclusions, found no error that affects U < 0.
- **Decoder:** bit-exact against an independent pure-Python decoder on 4 files, including the file at the 8·10⁷ seam and the file containing T₀. All 440 md5 sums match. The files have no trailing bytes, the t₀ values are exact integers, and blocks and files are contiguous from 78,146,000 to 1,000,046,000.
- **Zeros against Arb:** 826 zeros at all heights agree with `acb.zeta_zeros` to within 2⁻¹⁰² (maximum 0.9998·2⁻¹⁰²). Arb's `zeta_nzeros` confirms N(8·10⁷) = 195,565,267 and N(10⁹) = 2,846,548,032 (also OEIS A072080).
- **Mollifier series:** encloses the direct product on [0, 0.5] and beyond. The tail bound and its range are correct.
- **Independent float recomputation:** S_M = −1.00024458, giving U ≈ −1.967·10⁻⁴, inside the certified ball.
- **Radius:** most of the S_M radius (4.4·10⁻⁶ in U) comes from the 16-digit text caches for 2.3M ≤ n ≤ 13.96M. It costs about 2% of the margin.

**Hardening (applied 2026-10-25, certificate regenerated in 70 min, result unchanged).** The certificate itself now does the following:
- certifies N(8·10⁷) = 195,565,267 and N(10⁹) = 2,846,548,032 with Arb's `arb_zeta_nzeros`, independently of Platt's counts, and checks with `acb_dirichlet_zeta_zeros` that each cutoff lies strictly between consecutive zeros;
- checks that Platt's zero #195,565,267 overlaps our Arb zero #195,565,267, so that no zero is dropped or doubled at the seam;
- hashes `libplatt.so`;
- attaches radius 2⁻¹⁰¹ to Platt's zeros.

The regenerated certificate has U ≤ −1.956·10⁻⁴ and X < 10^{215.281005099}.
