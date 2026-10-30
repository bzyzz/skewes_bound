# skewes_verify

A Lean 4 + Mathlib (v4.35.0) check of the paper's main theorems. For each one, Lean derives the
stated bound from the paper's cited analytic results and its Arb-certified numbers. Lean also
proves the paper's kernel lemmas: the normalisation and transform of the two-scale kernel
(Lemma 3.1(a),(c) and eq. (3.4)).

```
lake exe cache get && lake build
lake env lean AxiomCheck.lean     # prints the axioms each theorem depends on
```

`#print axioms` reports only `propext`, `Classical.choice` and `Quot.sound` for every theorem
below. No file contains `sorry`, `admit`, `native_decide` or a new `axiom`. These standard axioms
do not discharge the hypotheses listed below; those remain inputs.

## The theorems

| Paper | Lean | Conclusion |
|---|---|---|
| Thm 1.1 | `Race.race` | `∃ x, log x ∈ [111.5456338, 111.5460262] ∧ π(x;24,11) < π(x;24,1) ∧ x < 10^48.4439` |
| Cor 1.2 | `Race.race24` | every `b ∈ {5,7,11,13,17,19,23}` has `x < 10^70.021` with `π(x;24,b) < π(x;24,1)` |
| Thm 1.3(a) | `Team.team8` | `∃ x, log x ∈ [64.764452, 64.764772] ∧ π(x;8,1) > (π(x;8,3)+π(x;8,5)+π(x;8,7))/3 ∧ x < 10^28.1270` |
| Thm 1.3(b) | `Team.team12` | the same modulo 12, `log x ∈ [339.631918, 339.632698]`, `x < 10^147.5007` |
| Thm 1.4 | `Mertens.mertens215` | `∃ x, log x ∈ [495.7028331177, 495.7028331467] ∧ Σ_{p≤x} log p/p < log x + E ∧ x < 10^215.2811` |
| Thm 1.5 | `Mertens.mertens807` | `∃ x, log x ∈ [1859.129042158, 1859.129060642] ∧ Δ₁(x) < 0 ∧ x < 10^807.4095` |
| Thm 1.6 | `Skewes.skewes` | `∃ x, log x ∈ [727.951335292, 727.951335596] ∧ li(x) < π(x) ∧ 2 ≤ x < 1.397166434·10^316` |

All the objects in these statements are genuine (`Defs.lean`). `π(x)` is `Nat.primeCounting ⌊x⌋₊`,
and `π(x;q,b)` counts the primes `p ≤ x` with `p % q = b`. `li` is the principal-value integral.
`E = -γ - Σ_p log p/(p(p-1))` uses Mathlib's `eulerMascheroniConstant`. The kernel `K` is the
convolution of eq. (3.3), built from `Real.cos`.

## What is proved in Lean, and what is assumed

| Step | Status |
|---|---|
| Lemma 3.1(a): `κ_{2m} ≥ 0`, `∫ κ_{2m} = 1` with `C_m = 4^m (m!)²/(2(2m)!)` (`Kernel.lean`) | proved |
| `K ≥ 0`, `∫ K = 1`, `K = 0` outside `[-η₀, η₀]`, so `∫_window K(u-ω) du = 1` (`Kernel.lean`) | proved |
| Lemma 3.1(c): `κ̂_{2m}(ξ) = (sin ξ/ξ) ∏ (kπ)²/((kπ)²-ξ²)`, `\|κ̂_{2m}\| ≤ 1`, and the envelope `Φ_{2m}` for `ξ > mπ` (`Fourier.lean`) | proved |
| Eq. (3.4): `K̂(γ) = κ̂_core(η₁γ) κ̂_{2M}(δγ)` (`Fourier.lean`) | proved |
| The averaging step, `∫ K Φ ≤ Φ(u_lo)`, and the exact value of the bias `8Φ(u_lo)`, `3Φ(u_lo)`, `Φ(ω)` (`Averaging.lean`, theorem files) | proved |
| The remainders `r(u)` of Prop. 5.3 and the prime-power term of Lemma 6.1, bounded on the window. These bounds are cruder than the paper's `r(u_lo)` and `R₂`, but they suffice (`Remainders.lean`) | proved |
| The windows `[ω-η₀, ω+η₀]`, every sign check `U < 0` / `L > 0` / `I > 0`, and all the conversions to powers of ten, including `e^{727.951335596} < 1.397166434·10^316` (`Decimal.lean`) | proved |
| The side conditions of Theorems 4.7, 5.5, 6.2 at the paper's parameters (`τ < 1`, `η₁T₀ > 3π`, …) (`SideConditions.lean`) | proved |
| `hExpl`: the explicit inequality (Lemma 4.3, §4.5, Prop. 5.3, Lemma 6.1): Landau / Riemann–von Mangoldt explicit formulae, Dusart, Platt–Trudgian `θ`, Rosser–Schoenfeld | **hypothesis** (cited + derived on paper) |
| `hZint`, `hZ`: the zero sum integrated against `K` is at most `S + R₃ + R₄ + R₅` (proofs of Theorems 4.7, 5.5, 6.2). These rely on term-by-term integration, (G)RH verification heights (Platt, Platt–Trudgian), zero counting (BMOR, Lehman), Lemma 2.1 and Lemmas 3.2–3.3 | **hypothesis** (cited + derived on paper) |
| `hS`, `hR3`, `hR4`, `hR5`, `hε`: the certified zero sums and error bounds | **hypothesis** (Arb certificates) |
| Corollary 1.2 for `b ≠ 11`: Ford–Hudson's windows | **hypothesis** (cited) |

In `hExpl` and `hZ`, the zero sum `Z(u)` is an arbitrary real function. Lean does not
identify it with `Σ_ρ …` over zeros of `ζ` or of Dirichlet `L`-functions. So Lean certifies that
each conclusion follows from these inputs. It does not certify the inputs. Lemma 3.1(c) and
eq. (3.4), which `hZ` and the Arb sums rely on, are proved.

Formalizing the inputs would need, in Lean: explicit formulae for `ψ(x, χ)` and `π(x)` with
truncation errors, the zeros of `ζ` and `L(s, χ)` with multiplicity and counting bounds, Lemma 3.2
(`2M`-fold integration by parts plus Cauchy estimates), and a verified evaluation of `K̂` at
millions of zeros. That is far beyond current Mathlib.

## Certified numbers used

Each numerical hypothesis is the upper (or lower) end of a printed Arb ball. Where a certificate
prints a rounded upper bound without a radius, the hypothesis raises the last printed digit by one.

| Theorem | Source | Hypotheses |
|---|---|---|
| 1.1 | `race/criterion/cert_T100000_w111.5458300.json` | `two_sum_S`, `R3`, `R4`, `R5`, `eps0` |
| 1.3(a) | `race/team/team_cert_q8_T100000_w64.764612.json` | `sum_S`, `R3`, `R4`, `R5`, `eps0` |
| 1.3(b) | `race/team/team_cert_q12_T50000_w339.632308.json` | `sum_S`, `R3`, `R4`, `R5`, `eps0` |
| 1.4 | `mertens/lmfdb/first_cert_lmfdb_T1000000000.json` | `S_M`, `R4_stieltjes`, `R5` |
| 1.5 | `mertens/criterion/first_cert_T1300000.json` | `S_M`, `R4_stieltjes`, `R5` |
| 1.6 | `certification/s_certificate_arb.json`, `certification/error_bounds.py` | `S_K`; `R₃`, `R₄,env`, `R₅` from Table 6 |

Not formalized: the `2.75·10^153` consecutive integers in Theorem 1.6 (which use [SD10, Thm 6.2]),
Theorem 7.2 (`θ(x) > x`, which does not improve the literature), and the conditional Theorem 8.1.

## Build notes

`.lake/packages` is a symlink to `../../mertens_bound/lean/.lake/packages`, which holds the
same Mathlib build (v4.35.0). On a fresh clone, delete the link and run `lake exe cache get`.
