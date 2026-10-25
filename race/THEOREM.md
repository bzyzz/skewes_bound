# A negative value of π(x;24,11) − π(x;24,1) below 10^48.5 — statement, inputs, proof outline

Status: **draft for expert checking.** The certificate is `criterion/cert_T100000_w111.5458300.json` (2026-10-25).

Each input is marked **[cited]**, **[derived]** (proved here) or **[computed]** (by a script in this repository).

Notation:
- Δ(x) := π(x;24,11) − π(x;24,1).
- X := {3.2, 4.3, 8.5, 24.11} (Conrey labels): the primitive real characters χ with χ(11) = −1. They have conductors q_χ = 3, 4, 8, 24 and parities δ_χ = 1, 1, 0, 0.
- For χ ∈ X, ρ = β + iγ runs over the nontrivial zeros of L(s, χ).

## 1. Statement

**Theorem.** π(x;24,11) < π(x;24,1) for some x with log x ∈ [111.5456338, 111.5460262]. Hence the first sign change of Δ satisfies x < e^{111.5460262} < 10^{48.4439}. Ford–Hudson (2001) gave 10³⁵³ for this race. Combined with their other mod-24 entries, every Δ_{24,b,1} now changes sign below 10^{70.02} (set by b = 5).

## 2. Inputs

- **(FH) [cited] Ford–Hudson (2001), (2.10)–(2.13).**
  - (2.10), **corrected**: for χ real primitive with no zeros in (0, 1),
    S(χ;x) := Σ_{p^m ≤ x} χ(p)^m/m = −Σ_ρ li(x^ρ) + ∫_x^∞ dy/(y^{1−δ}(y²−1) log y) − (1−δ) log log x + C_χ.
    - The integral lies in (0, 1/x) for x > 10.
    - The constant is C_χ = log L(0,χ) for odd χ and C_χ = log L′(0,χ) − γ_E for even χ. These are the residues of log L(s,χ)·x^s/s at s = 0 in Perron's formula, using (1/2πi)∫x^s log s ds/s = −log log x − γ_E.
    - FH print +(1−δ) log log x and, for even χ, K₀ = C − log(τπL(1,χ)/(2q)). Both disagree with this derivation and with a direct numerical test: the smoothed S(χ;x) computed from primes for x ≈ e^{16.5}–e^{18.1} gives C = −1.09, −0.70, −0.69, +0.2 ± 0.2 for 3.2, 4.3, 8.5, 24.11, against the predicted −1.099, −0.693, −0.704, +0.252. Their odd constant K₁ = log(τL(1)/(iπ)) = log L(0,χ) is correct.
    - The code (`C_const`) computes |C_χ| rigorously from Arb's Taylor series of L(s,χ) at s = 0: 1.0986, 0.6931, 0.7035, 0.2524. The values L′(0,χ₈) = log(1+√2) and L′(0,χ₂₄) = log(5+2√6) are the expected regulators.
    - The constants enter multiplied by u e^{−u/2} ≈ 7·10⁻²³.
  - (2.11): π_{q,a}(x) = φ(q)⁻¹ Σ_χ χ̄(a)S(χ;x) − Σ_{p^m ≤ x, m ≥ 2, p^m ≡ a} 1/m.
  - (2.12): Σ_{p^m ≤ x, m ≥ 3} 1/m ≤ 1.3x^{1/3}/log x for x ≥ e³⁰.
  - (2.13), **corrected**: |S(χ*;x) − S(χ;x)| ≤ Σ_{p | q, p ∤ q*} Σ_{m ≤ log x/log p} 1/m ≤ Σ_{p} (1 + log log x − log log p). This is ≤ log log x + 1.37 for p = 2 and ≤ log log x + 0.91 for p = 3.
    - FH's final step "≤ log log x + 1 − log 2" is false for p = 2, 3.
    - The code uses the corrected constant 1.37 (fixed 2026-10-25; the certificate was regenerated and is unchanged at printed precision).
- **(D) [cited] Dusart, Ramanujan J. 45 (2018), Thm 5.1, as quoted in Revers (J. Number Theory 2025, §4).** π(y) ≤ (y/log y)(1 + 1/log y + 2/log²y + 7.32/log³y) for y ≥ 4·10⁹.
  - The primary source has not been re-read here.
  - The weaker version with 7.59, valid for y > 1, would add 1.25·10⁻⁵ to the main term, which is harmless.
- **(P) [cited] Platt, Math. Comp. 85 (2016), arXiv:1305.3087.** GRH holds for all primitive characters of modulus q ≤ 400,000 up to height t₀ = max(10⁸/q, 7.5·10⁷/q + 200) for even q, and max(10⁸/q, 3.75·10⁷/q + 200) for odd q. That gives 3.33·10⁷, 2.5·10⁷, 1.25·10⁷ and 4.17·10⁶; we use the rounded-down values A_χ = 3.3·10⁷, 2.5·10⁷, 1.25·10⁷ and 4.1·10⁶. (Checked against Theorem 7.1 in the body of Platt's paper; since 10⁸/q dominates for q | 24, the even/odd distinction does not matter.)
- **(B) [cited] Bennett–Martin–O'Bryant–Rechnitzer, Math. Comp. 90 (2021), Thm 1.1.** Let N(T,χ) = #{ρ: |γ| ≤ T}, T ≥ 5/7 and ℓ = log(q(T+2)/2π) > 1.567. Then |N(T,χ) − (T/π) log(qT/2πe) + χ(−1)/4| ≤ 0.22737ℓ + 2 log(1+ℓ) − 0.5.
- **(T) [cited] Turing's method.**
  - Trudgian, Math. Comp. 80 (2011), Thm 3.3: |∫_{t₁}^{t₂} S(t,χ)dt| ≤ 1.975 + 0.084 log(qt₂/2π) for t₂ > t₁ > 50.
  - We also enforce Rumely's (1993) constants (1.8397, 0.1242) as recomputed by Trudgian, and use the larger of the two bounds.
  - BMOR flag an error in Trudgian's 2014 and 2015 papers, but not in the 2011 one. Our Turing margins exceed 10×.
- **(Z) [computed] Zeros.** `zeros/arb_certify.py` and `zeros/turing.py` give, for each χ ∈ X:
  - L(σ,χ) > 0 on [½, 2], so there are no real zeros in (0, 1) (by the functional equation) and N(T) = θ(T)/π + S(T);
  - disjoint intervals of width 2·10⁻⁸, each containing a sign change of Z(t,χ);
  - N(T₀,χ)/2 = F(T₀) at T₀ = 10⁵ (and every 5000 below).
  
  Consequently, the zeros with 0 < |γ| ≤ T₀ are exactly ½ ± iγ_j, simple, with each γ_j in its interval.

## 3. The explicit inequality [derived]

Put x = eᵘ with u ≥ 100. Units mod 24 satisfy p² ≡ 1, so no p² ≡ 11. Also χ(1) = 1, and the characters with χ(11) = χ(1) cancel. Combining (2.11), (2.13) and (2.10):

8Δ(x) = 2 Σ_{χ∈X} Σ_ρ li(x^ρ) + 4(π(√x) − 2) + ϑ(ε₁(x)).

Here ε₁(x) = 2Σ_X(1/x + (1−δ_χ) log log x + |C_χ|) + 2Σ_X n_χ(log log x + 1.37) + 8·1.3x^{1/3}/log x, with n_χ = 1, 1, 1, 0 (the primes dividing 24 but not q_χ).

Multiply by u e^{−u/2}, use (D) for π(e^{u/2}), and drop −8u e^{−u/2} < 0:

E(u) := 8u e^{−u/2} Δ(eᵘ) ≤ 8Φ_D(u) + 2 Σ_X Σ_ρ u e^{−u/2} li(e^{ρu}) + ε₀(u),

where Φ_D(u) = 1 + 2/u + 8/u² + 58.56/u³ and

ε₀(u) = u e^{−u/2}·[2Σ_X(e^{−u} + (1−δ_χ) log u + |C_χ|) + 2Σ_X n_χ(log u + 1.37)] + 10.4e^{−u/6}.

On the window we evaluate the decreasing factors at u_lo and the increasing ones at u_hi.

## 4. The criterion [derived]

The kernel K = K_core ∗ ψ_δ is the one in [`../THEOREM.md`](../THEOREM.md) §1: K ≥ 0, supp K = [−η₀, η₀], ∫K = 1, and K̂(γ) = κ̂_core(η₁γ)κ̂_{2M}(δγ). Define I(ω) := ∫K(u−ω)E(u)du. If I(ω) < 0, then E < 0 somewhere on [ω−η₀, ω+η₀].

**Interchange of Σ_ρ and ∫K.**
- The sum Σ_ρ li(x^ρ) converges only conditionally (symmetric partial sums over |γ| ≤ T).
- The truncated explicit formula for ψ(x,χ) (Davenport, Ch. 19, (13)–(14)) has an error O(x log²(qxT)/T + log x·min(1, x/(T⟨x⟩))), uniformly in x. Partial summation against dt/log t transfers this to S(χ;x).
- So the partial sums of Σ_ρ u e^{−u/2}li(e^{ρu}) are uniformly bounded on [u_lo, u_hi] and converge for u off the prime powers, a null set. Dominated convergence then gives ∫K·Σ_ρ = Σ_ρ ∫K.
- The same justification applies to the Skewes criterion in `../THEOREM.md`.

Integrating §3 against K gives I(ω) ≤ U(ω) := 8Φ_D(u_lo) + 2Σ_X S_χ(ω,T₀) + R₃ + R₄ + R₅ + ε₀, where S_χ = Σ_{0<|γ|≤T₀} K̂(γ)e^{iγω}(1/ρ + 1/(ρ²ω)). The error terms are as follows.

- **R₃ (zeros with |γ| ≤ A_χ, on the line by (Z) and (P)).**
  - Lehman's expansion li(e^w) = e^w(1/w + 1/w² + 2∫e^{−t}(w−t)^{−3}dt) gives u e^{−u/2}li(e^{ρu}) = e^{iγu}(1/ρ + 1/(ρ²u)) + ϑ(2/(|γ|³u²)).
  - In addition, |1/(ω+y) − 1/ω| ≤ |y|/u_lo² and ∫|y|K ≤ η₀/2.
  - Hence R₃ = (2/u_lo²)Σ_X[2Σ|γ|⁻³ + (η₀/2)Σγ⁻²], where the power sums are exact over the certified zeros plus a Stieltjes tail from (B).
- **R₄ (zeros with T₀ < |γ| ≤ A_χ).**
  - φ(t) := H(η₁t)(1/t + 1/(t²ω)) majorises |K̂(γ)(1/ρ + 1/(ρ²ω))|, where H = Σλ_mΦ_{2m} is nonincreasing for η₁t > 3π.
  - Stieltjes with (B) for the one-sided count N⁺ = N/2 gives Σ_{T₀<γ≤A}φ(γ) ≤ ∫φ·(1/2π)log(qt/2π)dt + E(T₀)φ(T₀) + ½∫φE′dt.
  - R₄ = 4Σ_X[that bound].
- **R₅ (zeros with |γ| > A_χ, possibly off the line).**
  - Put h_ρ(z) = z e^{−ρz}li(e^{ρz}), analytic for |z−u| ≤ R := u_lo/2. There |arg z| ≤ π/6, so |Im ρz| ≥ ½|ρz| and |h_ρ| ≤ 2/|γ|.
  - Integrating by parts 2M times (K ∈ C^{2M−1}, with ‖K^{(j)}‖₁ ≤ 2C_M(Mπ/δ)^j) and applying Cauchy and Leibniz gives |J_ρ| ≤ e^{u_hi/2}·4C_M|γ|^{−2M−1}(Mπ/δ + 2M/R)^{2M}.
  - Summing over |γ| > A_χ with N(t) ≤ (t/π)log(qt/2π) gives R₅.
- **Main term.** Φ_D is decreasing, so ∫K·8Φ_D ≤ 8Φ_D(u_lo).
- **Over-counting.** The power sums in R₃ run over all zeros, which is more than needed and therefore safe. For T₀ < |γ| ≤ A_χ, both R₃ and R₄ use β = ½ from (P).

Implemented in `criterion/race_criterion.py` (error terms) and `criterion/certify.py` (S_K in Arb over the certified zeros, every hypothesis checked).

## 5. Certified numbers [computed]

Parameters: T₀ = 10⁵, ω = 111.5458300, λ = (0, 0, 1) (core κ₆), η₁ = 1.2·10⁻⁴ (η₁T₀ = 12), δ = 7.62·10⁻⁵, M = 40, η₀ = 1.962·10⁻⁴.

| χ | zeros 0 < γ ≤ 10⁵ (Turing) | S_χ(ω, T₀) | Σ γ⁻² (all zeros) | Σ \|γ\|⁻³ |
|---|---|---|---|---|
| 3.2 | 155,553 | −0.6607 ± 7.4e-5 | 0.1134022 | 0.0072340 |
| 4.3 | 160,131 | −0.9986 ± 5.4e-5 | 0.1560297 | 0.0139330 |
| 8.5 | 171,163 | −0.5196 ± 5.9e-5 | 0.2365263 | 0.0269167 |
| 24.11 | 188,648 | −1.9038 ± 3.5e-5 | 0.5624742 | 0.1316496 |

| Term | Value (upper bound) |
|---|---|
| 2Σ_χ S_χ | −8.165 ± 8.0e-4 |
| main term 8Φ_D(u_lo) | 8.1489203 |
| R₃ | 5.78e-5 |
| R₄ (envelope) | 1.1305e-2 |
| R₅ | 9.13e-6 |
| ε₀ | 8.8e-8 |
| **U** | **−0.0050 ± 5.1e-4 < 0** |

Each zero is certified only to ±10⁻⁸. Tighter zero intervals, and the Stieltjes refinement of R₄, would widen the margin. They are not needed.

Validation:
- The float formula behind E(u) reproduces the known first sign changes of 8;5,1, 3;2,1, 24;13,1 and 8;7,1 (`heuristics/validate.py`).
- The Arb S_χ agrees with an independent float sum (to 1e-7 on a 5000-zero test).

Independent review (2026-10-25, fresh agent with no access to these conclusions): §3–§4, the code and the cited inputs were checked line by line, and the certified numbers were recomputed independently.
- Recomputed: U ≈ −0.00503, R₄ by Simpson quadrature 0.011255 (the code's 0.011305 is conservative), R₅ and the main term exactly.
- An explicit-formula test from actual primes confirms the sign and the coefficient +2 of the zero term.
- Findings: the two FH constants (corrected above; effect < 10⁻²¹); the interchange justification (added above); write-up details (added). One latent code issue: no check that T₀ < A_χ, which holds here.

Code fixes from the review (ε₀ with the true C_χ and the constant 1.37; `require(T0 < A_χ)`) were applied on 2026-10-25, and the certificate was regenerated: U = −0.0050 ± 5.05·10⁻⁴, the same as before at printed precision.

Open items before any claim:
- re-reading Dusart Thm 5.1 in the original;
- a reproduction from scratch;
- expert review.
