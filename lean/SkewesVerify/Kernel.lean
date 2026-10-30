import Mathlib

/-!
# The two-scale kernel (paper §3)

We define `κ_{2m}` (eq. (3.1)), the core and mollifier of eq. (3.3) and `K = K_core ⋆ ψ_δ`, and
prove the facts used by every criterion in the paper:

* Lemma 3.1(a): `κ_{2m} ≥ 0` and `∫ κ_{2m} = 1`, with `C_m = 4^m (m!)² / (2 (2m)!)`;
* `K ≥ 0`, `∫ K = 1`, and `K(y) = 0` for `|y| > η₀ = η₁ + δ`;
* hence `∫_{ω-η₀}^{ω+η₀} K(u - ω) du = 1` (`integral_window`).

The normalisation follows from Mathlib's reduction formula `integral_cos_pow`.
-/

open Real MeasureTheory Set

namespace SkewesVerify

/-- `C_m = 4^m (m!)² / (2 (2m)!)`. -/
noncomputable def C (m : ℕ) : ℝ := 4 ^ m * (m.factorial : ℝ) ^ 2 / (2 * ((2 * m).factorial : ℝ))

/-- `κ_{2m}(x) = C_m cos^{2m}(πx/2) 1_{[-1,1]}(x)`. -/
noncomputable def kappa (m : ℕ) : ℝ → ℝ :=
  indicator (Icc (-1) 1) fun x => C m * cos (π * x / 2) ^ (2 * m)

/-- `W m = ∫_{-π/2}^{π/2} cos^{2m}`. -/
noncomputable def W (m : ℕ) : ℝ := ∫ x in -(π / 2)..(π / 2), cos x ^ (2 * m)

theorem C_pos (m : ℕ) : 0 < C m := by unfold C; positivity

theorem C_succ (m : ℕ) : C (m + 1) = C m * ((2 * m + 2) / (2 * m + 1)) := by
  unfold C
  rw [show 2 * (m + 1) = (2 * m + 1) + 1 by ring, Nat.factorial_succ, Nat.factorial_succ,
    Nat.factorial_succ]
  push_cast
  field_simp
  ring

theorem W_zero : W 0 = π := by simp [W]

theorem W_succ (m : ℕ) : W (m + 1) = (2 * m + 1) / (2 * m + 2) * W m := by
  unfold W
  rw [show 2 * (m + 1) = 2 * m + 2 by ring, integral_cos_pow]
  simp [cos_neg, cos_pi_div_two]

/-- `C_m W_m = π / 2`, i.e. Lemma 3.1(a) in the variable `t = πx/2`. -/
theorem C_mul_W (m : ℕ) : C m * W m = π / 2 := by
  induction m with
  | zero => simp [C, W_zero]; ring
  | succ m ih =>
    rw [C_succ, W_succ]
    have h1 : (2 * (m : ℝ) + 1) ≠ 0 := by positivity
    have h2 : (2 * (m : ℝ) + 2) ≠ 0 := by positivity
    calc C m * ((2 * m + 2) / (2 * m + 1)) * ((2 * m + 1) / (2 * m + 2) * W m)
        = C m * W m := by field_simp
      _ = π / 2 := ih

theorem kappa_nonneg (m : ℕ) (x : ℝ) : 0 ≤ kappa m x := by
  unfold kappa
  refine indicator_nonneg (fun y _ => mul_nonneg (C_pos m).le ?_) x
  rw [pow_mul]; positivity

theorem kappa_eq_zero {m : ℕ} {x : ℝ} (hx : 1 < |x|) : kappa m x = 0 := by
  unfold kappa
  refine indicator_of_notMem (fun h => ?_) _
  exact absurd (abs_le.mpr ⟨h.1, h.2⟩) (not_le.mpr hx)

theorem integrable_kappa (m : ℕ) : Integrable (kappa m) :=
  ((by fun_prop : Continuous fun x : ℝ => C m * cos (π * x / 2) ^ (2 * m)).integrableOn_Icc
    ).integrable_indicator measurableSet_Icc

/-- **Lemma 3.1(a)**: `∫ κ_{2m} = 1`. -/
theorem integral_kappa (m : ℕ) : ∫ x, kappa m x = 1 := by
  unfold kappa
  rw [integral_indicator measurableSet_Icc, integral_Icc_eq_integral_Ioc,
    ← intervalIntegral.integral_of_le (by norm_num), intervalIntegral.integral_const_mul]
  have hsub : (∫ x in (-1 : ℝ)..1, cos (π * x / 2) ^ (2 * m)) = (π / 2)⁻¹ * W m := by
    have := intervalIntegral.integral_comp_mul_left (a := -1) (b := 1)
      (fun x => cos x ^ (2 * m)) (c := π / 2) (by positivity)
    simp only [smul_eq_mul] at this
    rw [show (π / 2) * (-1) = -(π / 2) by ring, mul_one] at this
    rw [W, ← this]
    congr 1; ext x; congr 2; ring
  rw [hsub, ← mul_assoc, mul_comm (C m), mul_assoc, C_mul_W]
  field_simp

/-! ### Scaling, the core and the mollifier -/

/-- `f_η(y) = η⁻¹ f(y/η)`. -/
noncomputable def scale (η : ℝ) (f : ℝ → ℝ) (y : ℝ) : ℝ := η⁻¹ * f (y / η)

theorem integrable_scale {η : ℝ} (hη : 0 < η) {f : ℝ → ℝ} (hf : Integrable f) :
    Integrable (scale η f) :=
  (hf.comp_div hη.ne').const_mul _

theorem integral_scale {η : ℝ} (hη : 0 < η) (f : ℝ → ℝ) : ∫ y, scale η f y = ∫ y, f y := by
  unfold scale
  rw [integral_const_mul, Measure.integral_comp_div, abs_of_pos hη, smul_eq_mul, ← mul_assoc,
    inv_mul_cancel₀ hη.ne', one_mul]

theorem scale_eq_zero {η : ℝ} (hη : 0 < η) {f : ℝ → ℝ} (hf : ∀ x, 1 < |x| → f x = 0) {y : ℝ}
    (hy : η < |y|) : scale η f y = 0 := by
  unfold scale
  rw [hf, mul_zero]
  rw [abs_div, abs_of_pos hη, lt_div_iff₀ hη, one_mul]
  exact hy

/-- The core profile `λ₁κ₂ + λ₂κ₄ + λ₃κ₆` (before scaling). -/
noncomputable def coreProfile (l₁ l₂ l₃ : ℝ) (x : ℝ) : ℝ :=
  l₁ * kappa 1 x + l₂ * kappa 2 x + l₃ * kappa 3 x

/-- The kernel parameters of eq. (3.3). -/
structure KParams where
  l₁ : ℝ
  l₂ : ℝ
  l₃ : ℝ
  η₁ : ℝ
  δ : ℝ
  M : ℕ
  hl₁ : 0 ≤ l₁
  hl₂ : 0 ≤ l₂
  hl₃ : 0 ≤ l₃
  hsum : l₁ + l₂ + l₃ = 1
  hη₁ : 0 < η₁
  hδ : 0 < δ

namespace KParams

variable (P : KParams)

/-- `η₀ = η₁ + δ`. -/
def η₀ : ℝ := P.η₁ + P.δ

/-- `K_core(y) = η₁⁻¹ Σ λ_m κ_{2m}(y/η₁)`. -/
noncomputable def core : ℝ → ℝ := scale P.η₁ (coreProfile P.l₁ P.l₂ P.l₃)

/-- `ψ_δ(y) = δ⁻¹ κ_{2M}(y/δ)`. -/
noncomputable def moll : ℝ → ℝ := scale P.δ (kappa P.M)

/-- `K = K_core ⋆ ψ_δ`, eq. (3.3). -/
noncomputable def K : ℝ → ℝ := convolution P.core P.moll (ContinuousLinearMap.lsmul ℝ ℝ) volume

theorem core_nonneg (y : ℝ) : 0 ≤ P.core y := by
  unfold core scale coreProfile
  have := P.hη₁
  have := P.hl₁; have := P.hl₂; have := P.hl₃
  have := kappa_nonneg 1 (y / P.η₁); have := kappa_nonneg 2 (y / P.η₁)
  have := kappa_nonneg 3 (y / P.η₁)
  positivity

theorem moll_nonneg (y : ℝ) : 0 ≤ P.moll y := by
  unfold moll scale
  have := P.hδ; have := kappa_nonneg P.M (y / P.δ)
  positivity

theorem integrable_core : Integrable P.core :=
  integrable_scale P.hη₁ ((((integrable_kappa 1).const_mul _).add
    ((integrable_kappa 2).const_mul _)).add ((integrable_kappa 3).const_mul _))

theorem integrable_moll : Integrable P.moll := integrable_scale P.hδ (integrable_kappa _)

theorem integral_core : ∫ y, P.core y = 1 := by
  unfold core
  rw [integral_scale P.hη₁]
  unfold coreProfile
  rw [integral_add (f := fun y => P.l₁ * kappa 1 y + P.l₂ * kappa 2 y)
      (g := fun y => P.l₃ * kappa 3 y)
      (((integrable_kappa 1).const_mul _).add ((integrable_kappa 2).const_mul _))
      ((integrable_kappa 3).const_mul _),
    integral_add (f := fun y => P.l₁ * kappa 1 y) (g := fun y => P.l₂ * kappa 2 y)
      ((integrable_kappa 1).const_mul _) ((integrable_kappa 2).const_mul _),
    integral_const_mul, integral_const_mul, integral_const_mul, integral_kappa, integral_kappa,
    integral_kappa]
  linarith [P.hsum]

theorem integral_moll : ∫ y, P.moll y = 1 := by
  unfold moll; rw [integral_scale P.hδ, integral_kappa]

theorem core_eq_zero {y : ℝ} (hy : P.η₁ < |y|) : P.core y = 0 :=
  scale_eq_zero P.hη₁ (fun x hx => by simp [coreProfile, kappa_eq_zero hx]) hy

theorem moll_eq_zero {y : ℝ} (hy : P.δ < |y|) : P.moll y = 0 :=
  scale_eq_zero P.hδ (fun _ hx => kappa_eq_zero hx) hy

theorem K_nonneg (y : ℝ) : 0 ≤ P.K y := by
  unfold K
  rw [convolution_def]
  exact integral_nonneg fun t => by
    simpa [ContinuousLinearMap.lsmul_apply, smul_eq_mul] using
      mul_nonneg (P.core_nonneg t) (P.moll_nonneg (y - t))

theorem integrable_K : Integrable P.K :=
  P.integrable_core.integrable_convolution _ P.integrable_moll

/-- `∫ K = 1`. -/
theorem integral_K : ∫ y, P.K y = 1 := by
  unfold K
  rw [integral_convolution _ P.integrable_core P.integrable_moll, P.integral_core,
    P.integral_moll]
  simp

/-- `supp K ⊆ [-η₀, η₀]`. -/
theorem K_eq_zero {y : ℝ} (hy : P.η₀ < |y|) : P.K y = 0 := by
  by_contra hne
  have hmem : y ∈ Function.support P.K := hne
  have hsub := support_convolution_subset (L := ContinuousLinearMap.lsmul ℝ ℝ) (μ := volume)
    (f := P.core) (g := P.moll) hmem
  obtain ⟨a, ha, b, hb, rfl⟩ := hsub
  have ha' : |a| ≤ P.η₁ := not_lt.mp fun h => ha (P.core_eq_zero h)
  have hb' : |b| ≤ P.δ := not_lt.mp fun h => hb (P.moll_eq_zero h)
  have := abs_add_le a b
  unfold η₀ at hy
  linarith

theorem η₀_pos : 0 < P.η₀ := add_pos P.hη₁ P.hδ

/-- The window `[ω - η₀, ω + η₀]`. -/
def window (ω : ℝ) : Set ℝ := Icc (ω - P.η₀) (ω + P.η₀)

theorem K_shift_eq_zero {ω u : ℝ} (hu : u ∉ P.window ω) : P.K (u - ω) = 0 := by
  refine P.K_eq_zero ?_
  have := P.η₀_pos
  unfold window at hu
  rw [mem_Icc, not_and_or, not_le, not_le] at hu
  rcases hu with h | h
  · rw [abs_of_neg (by linarith)]; linarith
  · rw [abs_of_pos (by linarith)]; linarith

/-- `∫_{ω-η₀}^{ω+η₀} K(u-ω) du = 1`. -/
theorem integral_window (ω : ℝ) : ∫ u in P.window ω, P.K (u - ω) = 1 := by
  rw [setIntegral_eq_integral_of_forall_compl_eq_zero fun u hu => P.K_shift_eq_zero hu,
    integral_sub_right_eq_self (fun y => P.K y), P.integral_K]

theorem integrable_K_shift (ω : ℝ) : Integrable fun u => P.K (u - ω) :=
  P.integrable_K.comp_sub_right ω

end KParams

end SkewesVerify
