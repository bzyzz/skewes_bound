import SkewesVerify.Race
import SkewesVerify.Team
import SkewesVerify.Mertens
import SkewesVerify.Skewes

/-!
# The side conditions of the criteria, at the paper's parameters

Theorems 4.7, 5.5 and 6.2 come with conditions on their parameters. These conditions are
premises of the cited criteria, not hypotheses of the Lean theorems. Here we check that the
parameters of §4.4, §4.5, §5.3 and §6.3 satisfy them.

* Race (Theorem 4.7): `u_lo ≥ 100`, `η₁T₀ > 3π`, `T₀ < min A_χ`, and
  `τ_χ = (Mπ/δ + 4M/u_lo) / A_χ < 1` for the smallest height `A_χ = 4.1·10⁶`.
* Team races: the same `τ` condition for `A = 1.25·10⁷` (mod 8) and `A = 8.3·10⁶` (mod 12),
  and `u_lo ≥ 45`.
* Mertens (Theorem 5.5): `u_lo ≥ 400`, `η₁T₀ > 3π`, `2πe < T₀ < A`, `τ = Mπ/(Aδ) < 1`.
* Skewes (Theorem 6.2): `ω - η₀ > 400`, `η₀ ≤ ω/100`, `2πe < T₀ ≤ A`, `η₁T₀ > 3π`, and
  `τ = Mπ/(Aδ) + 2M/(300A) < 1`.
-/

open Real

namespace SkewesVerify

theorem two_pi_e_lt {T : ℝ} (hT : 18 ≤ T) : 2 * π * exp 1 < T := by
  have := exp_one_lt_d9
  have := pi_lt_d2
  nlinarith [pi_pos, exp_pos 1]

theorem Race.side_conditions :
    100 ≤ ω - P.η₀ ∧ 3 * π < P.η₁ * 100000 ∧ (100000 : ℝ) < 4100000 ∧
      ((P.M : ℝ) * π / P.δ + 4 * P.M / (ω - P.η₀)) / 4100000 < 1 := by
  simp only [KParams.η₀, P, ω]
  have := pi_lt_d2
  refine ⟨by norm_num, by linarith, by norm_num, ?_⟩
  rw [div_lt_one (by norm_num)]; norm_num; nlinarith

theorem Team.side_conditions :
    45 ≤ ω8 - P8.η₀ ∧ 3 * π < P8.η₁ * 100000 ∧
      ((P8.M : ℝ) * π / P8.δ + 4 * P8.M / (ω8 - P8.η₀)) / 12500000 < 1 ∧
    45 ≤ ω12 - P12.η₀ ∧ 3 * π < P12.η₁ * 50000 ∧
      ((P12.M : ℝ) * π / P12.δ + 4 * P12.M / (ω12 - P12.η₀)) / 8300000 < 1 := by
  simp only [KParams.η₀, P8, P12, ω8, ω12]
  have := pi_lt_d2
  refine ⟨by norm_num, by linarith, ?_, by norm_num, by linarith, ?_⟩ <;>
    (rw [div_lt_one (by norm_num)]; norm_num; nlinarith)

theorem Mertens.side_conditions :
    (400 ≤ ω215 - P215.η₀ ∧ 3 * π < P215.η₁ * 10 ^ 9 ∧ 2 * π * exp 1 < 10 ^ 9 ∧
      (10 : ℝ) ^ 9 < 3 * 10 ^ 12 ∧ (P215.M : ℝ) * π / (3 * 10 ^ 12 * P215.δ) < 1) ∧
    (400 ≤ ω807 - P807.η₀ ∧ 3 * π < P807.η₁ * 1300000 ∧ 2 * π * exp 1 < 1300000 ∧
      (1300000 : ℝ) < 3 * 10 ^ 12 ∧ (P807.M : ℝ) * π / (3 * 10 ^ 12 * P807.δ) < 1) := by
  simp only [KParams.η₀, P215, P807, ω215, ω807]
  have := pi_lt_d2
  refine ⟨⟨by norm_num, by linarith, two_pi_e_lt (by norm_num), by norm_num, ?_⟩,
    ⟨by norm_num, by linarith, two_pi_e_lt (by norm_num), by norm_num, ?_⟩⟩ <;>
    (rw [div_lt_one (by norm_num)]; norm_num; nlinarith)

theorem Skewes.side_conditions :
    400 < ω - P.η₀ ∧ P.η₀ ≤ ω / 100 ∧ 2 * π * exp 1 < 80000000 ∧
      (80000000 : ℝ) ≤ 3 * 10 ^ 12 ∧ 3 * π < P.η₁ * 80000000 ∧
      (P.M : ℝ) * π / (3 * 10 ^ 12 * P.δ) + 2 * P.M / (300 * (3 * 10 ^ 12)) < 1 := by
  simp only [KParams.η₀, P, ω]
  have := pi_lt_d2
  refine ⟨by norm_num, by norm_num, two_pi_e_lt (by norm_num), by norm_num, by linarith, ?_⟩
  norm_num
  nlinarith

end SkewesVerify
