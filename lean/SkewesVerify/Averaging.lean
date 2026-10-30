import SkewesVerify.Kernel
import SkewesVerify.Defs

/-!
# The averaging step shared by all criteria

If `f ≤ g` almost everywhere on the window and `∫ K(u - ω) g(u) du < 0`, then `f(u) < 0` for some
`u` in the window, because `K ≥ 0` (and symmetrically for `>`). No integrability of `f` is needed.

We also prove the main-term bound `∫ K(u - ω) Φ(u) du ≤ Φ(ω - η₀)` (`Φ` is decreasing and
`∫ K = 1`), and a crude bound `e^{-x} ≤ (1 + x/n)^{-n}` used for the small remainders.
-/

open Real MeasureTheory Set

namespace SkewesVerify

theorem exists_neg_of_integral_neg {w f g : ℝ → ℝ} {s : Set ℝ} (hs : MeasurableSet s)
    (hw : ∀ u ∈ s, 0 ≤ w u) (hfg : ∀ᵐ u, u ∈ s → f u ≤ g u)
    (h : ∫ u in s, w u * g u < 0) : ∃ u ∈ s, f u < 0 := by
  by_contra hne
  push Not at hne
  have : 0 ≤ ∫ u in s, w u * g u := setIntegral_nonneg_ae hs (hfg.mono fun u hu hus =>
    mul_nonneg (hw u hus) ((hne u hus).trans (hu hus)))
  linarith

theorem exists_pos_of_integral_pos {w f g : ℝ → ℝ} {s : Set ℝ} (hs : MeasurableSet s)
    (hw : ∀ u ∈ s, 0 ≤ w u) (hfg : ∀ᵐ u, u ∈ s → g u ≤ f u)
    (h : 0 < ∫ u in s, w u * g u) : ∃ u ∈ s, 0 < f u := by
  obtain ⟨u, hu, hf⟩ := exists_neg_of_integral_neg (f := fun u => -f u) (g := fun u => -g u)
    (w := w) hs hw (hfg.mono fun _ hu hus => neg_le_neg (hu hus)) (by
      simp only [mul_neg, integral_neg]; linarith)
  exact ⟨u, hu, by linarith⟩

namespace KParams

variable (P : KParams)

theorem measurableSet_window (ω : ℝ) : MeasurableSet (P.window ω) := measurableSet_Icc

theorem K_shift_nonneg (ω u : ℝ) : 0 ≤ P.K (u - ω) := P.K_nonneg _

theorem integrableOn_K_mul {ω : ℝ} {g : ℝ → ℝ} (hg : ContinuousOn g (P.window ω)) :
    IntegrableOn (fun u => P.K (u - ω) * g u) (P.window ω) :=
  (P.integrable_K_shift ω).integrableOn.mul_continuousOn hg isCompact_Icc

theorem integrableOn_K_const (ω c : ℝ) :
    IntegrableOn (fun u => P.K (u - ω) * c) (P.window ω) :=
  P.integrableOn_K_mul continuousOn_const

theorem integral_K_add {ω : ℝ} {f g : ℝ → ℝ}
    (hf : IntegrableOn (fun u => P.K (u - ω) * f u) (P.window ω))
    (hg : IntegrableOn (fun u => P.K (u - ω) * g u) (P.window ω)) :
    ∫ u in P.window ω, P.K (u - ω) * (f u + g u) =
      (∫ u in P.window ω, P.K (u - ω) * f u) + ∫ u in P.window ω, P.K (u - ω) * g u := by
  rw [← integral_add hf hg]; congr 1; ext u; ring

theorem integrableOn_K_add {ω : ℝ} {f g : ℝ → ℝ}
    (hf : IntegrableOn (fun u => P.K (u - ω) * f u) (P.window ω))
    (hg : IntegrableOn (fun u => P.K (u - ω) * g u) (P.window ω)) :
    IntegrableOn (fun u => P.K (u - ω) * (f u + g u)) (P.window ω) :=
  (hf.add hg).congr_fun (fun u _ => by simp only [Pi.add_apply]; ring) (P.measurableSet_window ω)

theorem integral_K_const_mul {ω : ℝ} (c : ℝ) (f : ℝ → ℝ) :
    ∫ u in P.window ω, P.K (u - ω) * (c * f u) = c * ∫ u in P.window ω, P.K (u - ω) * f u := by
  rw [← integral_const_mul]; congr 1; ext u; ring

theorem integrableOn_K_const_mul {ω : ℝ} (c : ℝ) {f : ℝ → ℝ}
    (hf : IntegrableOn (fun u => P.K (u - ω) * f u) (P.window ω)) :
    IntegrableOn (fun u => P.K (u - ω) * (c * f u)) (P.window ω) :=
  IntegrableOn.congr_fun (Integrable.const_mul hf c) (fun u _ => by ring)
    (P.measurableSet_window ω)

theorem integral_K_const (ω c : ℝ) : ∫ u in P.window ω, P.K (u - ω) * c = c := by
  rw [integral_mul_const, P.integral_window, one_mul]

theorem continuousOn_Phi_window {ω : ℝ} (h : 0 < ω - P.η₀) : ContinuousOn Phi (P.window ω) :=
  continuousOn_Phi.mono fun _ hu => lt_of_lt_of_le h hu.1

/-- `∫ K(u - ω) Φ(u) du ≤ Φ(ω - η₀)`. -/
theorem integral_K_Phi_le {ω : ℝ} (h : 0 < ω - P.η₀) :
    ∫ u in P.window ω, P.K (u - ω) * Phi u ≤ Phi (ω - P.η₀) := by
  calc ∫ u in P.window ω, P.K (u - ω) * Phi u
      ≤ ∫ u in P.window ω, P.K (u - ω) * Phi (ω - P.η₀) := by
        refine setIntegral_mono_on (P.integrableOn_K_mul (P.continuousOn_Phi_window h))
          (P.integrableOn_K_const ω _) (P.measurableSet_window ω) fun u hu => ?_
        refine mul_le_mul_of_nonneg_left ?_ (P.K_nonneg _)
        exact Phi_antitoneOn (mem_Ioi.mpr h) (mem_Ioi.mpr (lt_of_lt_of_le h hu.1)) hu.1
    _ = Phi (ω - P.η₀) := P.integral_K_const ω _

/-- If `g ≤ c` on the window, then `∫ K(u - ω) g(u) du ≤ c`. -/
theorem integral_K_le_of_le {ω c : ℝ} {g : ℝ → ℝ} (hg : ContinuousOn g (P.window ω))
    (hc : ∀ u ∈ P.window ω, g u ≤ c) : ∫ u in P.window ω, P.K (u - ω) * g u ≤ c := by
  calc ∫ u in P.window ω, P.K (u - ω) * g u ≤ ∫ u in P.window ω, P.K (u - ω) * c :=
        setIntegral_mono_on (P.integrableOn_K_mul hg) (P.integrableOn_K_const ω _)
          (P.measurableSet_window ω) fun u hu => mul_le_mul_of_nonneg_left (hc u hu) (P.K_nonneg _)
    _ = c := P.integral_K_const ω _

end KParams

/-- `e^{-x} ≤ (1 + x/n)^{-n}` for `x ≥ 0`. -/
theorem exp_neg_le_inv_pow {x : ℝ} (hx : 0 ≤ x) (n : ℕ) (hn : 0 < n) :
    exp (-x) ≤ 1 / (1 + x / n) ^ n := by
  have hpos : 0 < 1 + x / n := by positivity
  have h1 : (1 + x / n) ^ n ≤ exp x := by
    calc (1 + x / n) ^ n ≤ exp (x / n) ^ n :=
          pow_le_pow_left₀ hpos.le (by linarith [add_one_le_exp (x / n)]) n
      _ = exp x := by rw [← exp_nat_mul]; congr 1; field_simp
  rw [exp_neg, one_div]
  exact inv_anti₀ (by positivity) h1

end SkewesVerify
