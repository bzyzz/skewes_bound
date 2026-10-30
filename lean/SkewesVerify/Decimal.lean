import Mathlib

/-!
# Decimal forms of the bounds

The paper states each window as `log x ∈ [a, b]` and then converts `e^b` to a power of ten
(`e^{111.5460262} < 10^{48.4439}` and so on). This file proves those conversions.

The only input is Mathlib's Taylor bound `Real.exp_bound'` for `exp` on `[0, 1]`. From it we get
`log 10 > 2.30258509299404` and `log 1.397166434 > 0.3344462099`, both checked with `norm_num`.
-/

open Real Finset

namespace SkewesVerify

/-- A lower bound for `log 10` with fourteen correct decimals. -/
noncomputable def L10 : ℝ := 230258509299404 / 10 ^ 14

theorem exp_L10_div_eight_le :
    exp (L10 / 8) ≤ (∑ m ∈ range 22, (L10 / 8) ^ m / m.factorial) +
      (L10 / 8) ^ 22 * (22 + 1) / ((22).factorial * 22) :=
  Real.exp_bound' (by unfold L10; norm_num) (by unfold L10; norm_num) (by norm_num)

theorem exp_L10_lt : exp L10 < 10 := by
  have h8 : exp L10 = exp (L10 / 8) ^ 8 := by
    rw [← exp_nat_mul]; congr 1; ring
  rw [h8]
  refine lt_of_le_of_lt (pow_le_pow_left₀ (exp_pos _).le exp_L10_div_eight_le 8) ?_
  unfold L10
  norm_num [Finset.sum_range_succ, Nat.factorial]

theorem L10_lt_log_ten : L10 < log 10 :=
  (lt_log_iff_exp_lt (by norm_num)).mpr exp_L10_lt

/-- `e^a < 10^b` as soon as `a ≤ b · 2.30258509299404`. -/
theorem exp_lt_ten_rpow {a b : ℝ} (hb : 0 < b) (h : a ≤ b * L10) : exp a < (10 : ℝ) ^ b := by
  rw [rpow_def_of_pos (by norm_num)]
  exact exp_lt_exp.mpr (by nlinarith [L10_lt_log_ten])

/-- A lower bound for `log 1.397166434`. -/
noncomputable def L1397 : ℝ := 3344462099 / 10 ^ 10

theorem exp_L1397_lt : exp L1397 < 1397166434 / 10 ^ 9 := by
  refine lt_of_le_of_lt (Real.exp_bound' (n := 20) (by unfold L1397; norm_num)
    (by unfold L1397; norm_num) (by norm_num)) ?_
  unfold L1397
  norm_num [Finset.sum_range_succ, Nat.factorial]

/-- `e^a < 1.397166434 · 10^316` as soon as `a ≤ 0.3344462099 + 316 · 2.30258509299404`. -/
theorem exp_lt_skewes_decimal {a : ℝ} (h : a ≤ L1397 + 316 * L10) :
    exp a < 1397166434 / 10 ^ 9 * 10 ^ 316 := by
  have h316 : exp (316 * L10) < 10 ^ 316 := by
    rw [show (316 : ℝ) * L10 = ((316 : ℕ) : ℝ) * L10 by norm_num, exp_nat_mul]
    exact pow_lt_pow_left₀ exp_L10_lt (exp_pos _).le (by norm_num)
  calc exp a ≤ exp (L1397 + 316 * L10) := exp_le_exp.mpr h
    _ = exp L1397 * exp (316 * L10) := exp_add _ _
    _ < 1397166434 / 10 ^ 9 * 10 ^ 316 :=
        mul_lt_mul'' exp_L1397_lt h316 (exp_pos _).le (exp_pos _).le

end SkewesVerify
