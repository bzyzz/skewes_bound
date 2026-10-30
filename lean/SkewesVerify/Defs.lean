import Mathlib

/-!
# The objects in the theorems

All of these are Mathlib objects or direct definitions from them:

* `piR x = π(x)`, via `Nat.primeCounting`;
* `piMod q b x = π(x; q, b)`, the number of primes `p ≤ x` with `p ≡ b (mod q)`;
* `li x`, the principal value `lim_{ε→0⁺} (∫_0^{1-ε} + ∫_{1+ε}^x) dt / log t`;
* `mertensE = -γ - Σ_p log p / (p (p - 1))` and `Delta1 x = Σ_{p ≤ x} log p / p - log x - E`;
* `Phi u = 1 + 2/u + 8/u² + 58.56/u³` (eq. (2.1)).
-/

open Real Finset Filter

namespace SkewesVerify

/-- `π(x)`. -/
noncomputable def piR (x : ℝ) : ℝ := (Nat.primeCounting ⌊x⌋₊ : ℝ)

/-- `π(x; q, b)`: the number of primes `p ≤ x` with `p % q = b`. -/
noncomputable def piMod (q b : ℕ) (x : ℝ) : ℝ :=
  (((range (⌊x⌋₊ + 1)).filter fun p => p.Prime ∧ p % q = b).card : ℝ)

/-- `li(x)`, as a principal value at `t = 1`. -/
noncomputable def li (x : ℝ) : ℝ :=
  limUnder (nhdsWithin 0 (Set.Ioi 0)) fun ε : ℝ =>
    (∫ t in (0 : ℝ)..(1 - ε), 1 / log t) + ∫ t in (1 + ε)..x, 1 / log t

/-- Mertens' constant `E = -γ - Σ_p log p / (p (p - 1))` in `Σ_{p ≤ x} log p / p = log x + E + o(1)`. -/
noncomputable def mertensE : ℝ :=
  -eulerMascheroniConstant - ∑' p : Nat.Primes, log (p : ℝ) / ((p : ℝ) * ((p : ℝ) - 1))

/-- `Δ₁(x) = Σ_{p ≤ x} log p / p - log x - E`. -/
noncomputable def Delta1 (x : ℝ) : ℝ :=
  (∑ p ∈ (range (⌊x⌋₊ + 1)).filter Nat.Prime, log (p : ℝ) / (p : ℝ)) - log x - mertensE

/-- `Φ(u) = 1 + 2/u + 8/u² + 58.56/u³`, eq. (2.1). -/
noncomputable def Phi (u : ℝ) : ℝ := 1 + 2 / u + 8 / u ^ 2 + 5856 / 100 / u ^ 3

theorem Phi_antitoneOn : AntitoneOn Phi (Set.Ioi 0) := by
  intro a ha b hb hab
  simp only [Set.mem_Ioi] at ha hb
  unfold Phi
  have h1 : 2 / b ≤ 2 / a := div_le_div_of_nonneg_left (by norm_num) ha hab
  have h2 : 8 / b ^ 2 ≤ 8 / a ^ 2 :=
    div_le_div_of_nonneg_left (by norm_num) (by positivity) (pow_le_pow_left₀ ha.le hab 2)
  have h3 : 5856 / 100 / b ^ 3 ≤ 5856 / 100 / a ^ 3 :=
    div_le_div_of_nonneg_left (by norm_num) (by positivity) (pow_le_pow_left₀ ha.le hab 3)
  linarith

theorem one_lt_Phi {u : ℝ} (hu : 0 < u) : 1 < Phi u := by
  unfold Phi; have : 0 < 2 / u := by positivity
  have : 0 ≤ 8 / u ^ 2 := by positivity
  have : 0 ≤ 5856 / 100 / u ^ 3 := by positivity
  linarith

theorem continuousOn_Phi : ContinuousOn Phi (Set.Ioi 0) := by
  unfold Phi
  refine ContinuousOn.add (ContinuousOn.add (ContinuousOn.add continuousOn_const ?_) ?_) ?_ <;>
    refine ContinuousOn.div continuousOn_const (by fun_prop) fun u hu => ?_ <;>
    simp only [Set.mem_Ioi] at hu <;> positivity

end SkewesVerify
