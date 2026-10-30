import SkewesVerify.Averaging

/-!
# The explicit remainders

* `mertensR u`, the term `r(u)` of Proposition 5.3 with `x = e^u`:
  `e^{-5u/2} + log(2π) e^{-u/2} + 1.53 e^{-u/6} + 1.36 (u / log 2 + 2) e^{-u/4}`;
* `skewesRem u`, the last term of Lemma 6.1:
  `u (log 2 + (2 / log 2) e^{u/3}) / (e^{u/2} Φ(u))`.

On a window `[lo, hi]` with `lo > 0` we bound each of them by an explicit constant times
`(1 + lo/120)^{-20}`. The paper's `R₂` and `r(u_lo)` are sharper, but these bounds suffice here.
-/

open Real Set

namespace SkewesVerify

/-- `r(u)` of Proposition 5.3, with `x = e^u`. -/
noncomputable def mertensR (u : ℝ) : ℝ :=
  exp (-(5 * u / 2)) + log (2 * π) * exp (-(u / 2)) + 153 / 100 * exp (-(u / 6)) +
    136 / 100 * (u / log 2 + 2) * exp (-(u / 4))

/-- The last term of Lemma 6.1. -/
noncomputable def skewesRem (u : ℝ) : ℝ :=
  u * (log 2 + 2 / log 2 * exp (u / 3)) / (exp (u / 2) * Phi u)

theorem log_two_pi_le_two : log (2 * π) ≤ 2 := by
  rw [log_le_iff_le_exp (by positivity)]
  have he : (2.7182818283 : ℝ) ^ 2 ≤ exp 2 := by
    rw [show (2 : ℝ) = 1 + 1 by norm_num, exp_add, ← sq]
    exact pow_le_pow_left₀ (by norm_num) exp_one_gt_d9.le 2
  nlinarith [pi_lt_d2]

theorem exp_neg_div_six_le {u : ℝ} (hu : 0 ≤ u) : exp (-(u / 6)) ≤ 1 / (1 + u / 120) ^ 20 := by
  have := exp_neg_le_inv_pow (x := u / 6) (by positivity) 20 (by norm_num)
  convert this using 3; push_cast; ring

theorem exp_neg_le_of_le {a b : ℝ} (h : a ≤ b) : exp (-b) ≤ exp (-a) := exp_le_exp.mpr (by linarith)

theorem continuousOn_mertensR {s : Set ℝ} : ContinuousOn mertensR s := by
  unfold mertensR; fun_prop

theorem continuousOn_skewesRem {s : Set ℝ} (hs : ∀ u ∈ s, 0 < u) : ContinuousOn skewesRem s := by
  unfold skewesRem
  refine ContinuousOn.div (by fun_prop) ((continuous_exp.comp (by fun_prop)).continuousOn.mul
    (continuousOn_Phi.mono fun u hu => hs u hu)) fun u hu => ?_
  exact (mul_pos (exp_pos _) (lt_trans one_pos (one_lt_Phi (hs u hu)))).ne'

/-- `r(u) ≤ (5 + 1.36 (hi/0.69 + 2)) (1 + lo/120)^{-20}` for `0 ≤ lo ≤ u ≤ hi`. -/
theorem mertensR_le {lo hi u : ℝ} (hlo : 0 ≤ lo) (hu : u ∈ Icc lo hi) :
    mertensR u ≤ (5 + 136 / 100 * (hi / (69 / 100) + 2)) * (1 / (1 + lo / 120) ^ 20) := by
  have hu0 : 0 ≤ u := hlo.trans hu.1
  set e := exp (-(u / 6))
  have he0 : 0 ≤ e := (exp_pos _).le
  have h1 : exp (-(5 * u / 2)) ≤ e := exp_neg_le_of_le (by linarith)
  have h2 : exp (-(u / 2)) ≤ e := exp_neg_le_of_le (by linarith)
  have h4 : exp (-(u / 4)) ≤ e := exp_neg_le_of_le (by linarith)
  have hl2 : 69 / 100 < log 2 := lt_trans (by norm_num) log_two_gt_d9
  have hq : u / log 2 ≤ hi / (69 / 100) := by
    calc u / log 2 ≤ hi / log 2 := div_le_div_of_nonneg_right hu.2 (by linarith)
      _ ≤ hi / (69 / 100) := div_le_div_of_nonneg_left (hlo.trans (hu.1.trans hu.2))
          (by norm_num) hl2.le
  have hq0 : 0 ≤ u / log 2 + 2 := by have : 0 ≤ u / log 2 := div_nonneg hu0 (by linarith); linarith
  have hlog := log_two_pi_le_two
  have hlogpos : 0 ≤ log (2 * π) := log_nonneg (by nlinarith [pi_gt_three])
  have hmono : e ≤ 1 / (1 + lo / 120) ^ 20 := by
    calc e ≤ exp (-(lo / 6)) := exp_neg_le_of_le (by linarith [hu.1])
      _ ≤ _ := exp_neg_div_six_le hlo
  have hsum : mertensR u ≤ (5 + 136 / 100 * (hi / (69 / 100) + 2)) * e := by
    unfold mertensR
    have := mul_le_mul hlog h2 (exp_pos _).le (by norm_num)
    have := mul_le_mul (mul_le_mul_of_nonneg_left (by linarith : u / log 2 + 2 ≤ hi / (69 / 100) + 2)
      (by norm_num : (0 : ℝ) ≤ 136 / 100)) h4 (exp_pos _).le (mul_nonneg (by norm_num) (by linarith))
    nlinarith
  calc mertensR u ≤ _ := hsum
    _ ≤ _ := by
        have hhi : 0 ≤ hi := hlo.trans (hu.1.trans hu.2)
        exact mul_le_mul_of_nonneg_left hmono (by positivity)

/-- The remainder of Lemma 6.1 is at most `4 hi (1 + lo/120)^{-20}` for `0 < lo ≤ u ≤ hi`. -/
theorem skewesRem_le {lo hi u : ℝ} (hlo : 0 < lo) (hu : u ∈ Icc lo hi) :
    skewesRem u ≤ 4 * hi * (1 / (1 + lo / 120) ^ 20) := by
  have hu0 : 0 < u := lt_of_lt_of_le hlo hu.1
  have hΦ := one_lt_Phi hu0
  have hl2 : 69 / 100 < log 2 := lt_trans (by norm_num) log_two_gt_d9
  have hl2' : log 2 < 1 := lt_trans log_two_lt_d9 (by norm_num)
  have hc : 2 / log 2 ≤ 3 := by rw [div_le_iff₀ (by linarith)]; linarith
  have hnum : 0 ≤ log 2 + 2 / log 2 * exp (u / 3) := by
    have : 0 < log 2 := by linarith
    positivity
  -- divide by `Φ(u) > 1` and split `e^{u/3} / e^{u/2} = e^{-u/6}`
  have step1 : skewesRem u ≤ u * (log 2 + 2 / log 2 * exp (u / 3)) * exp (-(u / 2)) := by
    unfold skewesRem
    rw [div_le_iff₀ (mul_pos (exp_pos _) (by linarith)), mul_assoc _ (exp _), ← mul_assoc
      (exp _), ← exp_add, show -(u / 2) + u / 2 = 0 by ring, exp_zero, one_mul]
    exact le_mul_of_one_le_right (mul_nonneg hu0.le hnum) hΦ.le
  have step2 : u * (log 2 + 2 / log 2 * exp (u / 3)) * exp (-(u / 2)) =
      u * (log 2 * exp (-(u / 2)) + 2 / log 2 * exp (-(u / 6))) := by
    rw [show -(u / 6) = u / 3 + -(u / 2) by ring, exp_add]; ring
  have he : exp (-(u / 2)) ≤ exp (-(u / 6)) := exp_neg_le_of_le (by linarith)
  have hmono : exp (-(u / 6)) ≤ 1 / (1 + lo / 120) ^ 20 :=
    (exp_neg_le_of_le (by linarith [hu.1])).trans (exp_neg_div_six_le hlo.le)
  have step3 : log 2 * exp (-(u / 2)) + 2 / log 2 * exp (-(u / 6)) ≤ 4 * exp (-(u / 6)) := by
    have := exp_pos (-(u / 6))
    nlinarith [mul_le_mul hl2'.le he (exp_pos _).le zero_le_one]
  calc skewesRem u ≤ _ := step1
    _ = _ := step2
    _ ≤ hi * (4 * exp (-(u / 6))) :=
        mul_le_mul hu.2 step3 (by positivity) (hlo.le.trans (hu.1.trans hu.2))
    _ ≤ hi * (4 * (1 / (1 + lo / 120) ^ 20)) :=
        mul_le_mul_of_nonneg_left (by linarith) (hlo.le.trans (hu.1.trans hu.2))
    _ = 4 * hi * (1 / (1 + lo / 120) ^ 20) := by ring

end SkewesVerify
