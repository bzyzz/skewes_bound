import SkewesVerify.Kernel

/-!
# The transform of the kernel (Lemma 3.1(c) and eq. (3.4))

Every certified zero sum in the paper evaluates `K̂(γ)` through the closed form
`κ̂_{2m}(ξ) = (sin ξ / ξ) ∏_{k=1}^m (kπ)² / ((kπ)² - ξ²)` (eq. (3.2)) and the factorisation
`K̂(γ) = κ̂_core(η₁ γ) κ̂_{2M}(δ γ)` (eq. (3.4)). This file proves both.

`I_k(ξ) = ∫_{-1}^1 cos^{2k}(πx/2) cos(ξx) dx` satisfies
`((kπ)² - ξ²) I_k = k (2k - 1) π² / 2 · I_{k-1}`, as in the paper's proof: with `c = cos(πx/2)`,
`s = sin(πx/2)`, the function `H = -kπ c^{2k-1} s cos(ξx) + ξ c^{2k} sin(ξx)` has derivative
`((ξ² - k²π²) c^{2k} + k(2k-1)π²/2 · c^{2k-2}) cos(ξx)` and vanishes at `±1`.
-/

open Real MeasureTheory Set

namespace SkewesVerify

/-- `I_k(ξ) = ∫_{-1}^1 cos^{2k}(πx/2) cos(ξx) dx`. -/
noncomputable def Icos (k : ℕ) (ξ : ℝ) : ℝ :=
  ∫ x in (-1 : ℝ)..1, cos (π * x / 2) ^ (2 * k) * cos (ξ * x)

theorem Icos_recurrence (j : ℕ) (ξ : ℝ) :
    (((j + 1 : ℝ) * π) ^ 2 - ξ ^ 2) * Icos (j + 1) ξ =
      (j + 1) * (2 * j + 1) * π ^ 2 / 2 * Icos j ξ := by
  set c := fun x : ℝ => cos (π * x / 2)
  set s := fun x : ℝ => sin (π * x / 2)
  let H : ℝ → ℝ := fun x =>
    -((j + 1 : ℝ) * π) * c x ^ (2 * j + 1) * s x * cos (ξ * x) + ξ * c x ^ (2 * j + 2) * sin (ξ * x)
  let h' : ℝ → ℝ := fun x =>
    (ξ ^ 2 - ((j + 1 : ℝ) * π) ^ 2) * (c x ^ (2 * (j + 1)) * cos (ξ * x)) +
      (j + 1) * (2 * j + 1) * π ^ 2 / 2 * (c x ^ (2 * j) * cos (ξ * x))
  have hderiv : ∀ x, HasDerivAt H (h' x) x := by
    intro x
    have hlin : HasDerivAt (fun x : ℝ => π * x / 2) (π / 2) x := by
      simpa using ((hasDerivAt_id x).const_mul π).div_const 2
    have hc : HasDerivAt c (-sin (π * x / 2) * (π / 2)) x := hlin.cos
    have hs : HasDerivAt s (cos (π * x / 2) * (π / 2)) x := hlin.sin
    have hξ : HasDerivAt (fun x : ℝ => ξ * x) ξ x := by simpa using (hasDerivAt_id x).const_mul ξ
    have h1 := ((((hc.pow (2 * j + 1)).const_mul (-((j + 1 : ℝ) * π))).mul hs).mul hξ.cos)
    have h2 := ((hc.pow (2 * j + 2)).const_mul ξ).mul hξ.sin
    convert h1.add h2 using 1
    have hsq : sin (π * x / 2) ^ 2 = 1 - cos (π * x / 2) ^ 2 := sin_sq _
    simp only [h', c, s, Pi.mul_apply, Pi.pow_apply, Nat.add_sub_cancel,
      show 2 * j + 2 - 1 = 2 * j + 1 by omega]
    push_cast
    linear_combination (-((j + 1 : ℝ) * π) * (2 * j + 1) * (π / 2)) *
      (cos (π * x / 2) ^ (2 * j) * cos (ξ * x)) * hsq
  have hc1 : c 1 = 0 := by simp [c]
  have hcm1 : c (-1) = 0 := by
    simp only [c]; rw [show π * -1 / 2 = -(π / 2) by ring, cos_neg, cos_pi_div_two]
  have hH1 : H 1 = 0 := by simp only [H, hc1]; simp
  have hHm1 : H (-1) = 0 := by simp only [H, hcm1]; simp
  have hint : ∫ x in (-1 : ℝ)..1, h' x = 0 := by
    rw [intervalIntegral.integral_eq_sub_of_hasDerivAt (fun x _ => hderiv x)
      (by apply Continuous.intervalIntegrable; fun_prop), hH1, hHm1, sub_zero]
  have hi1 : IntervalIntegrable (fun x => c x ^ (2 * (j + 1)) * cos (ξ * x)) volume (-1) 1 := by
    apply Continuous.intervalIntegrable; fun_prop
  have hi0 : IntervalIntegrable (fun x => c x ^ (2 * j) * cos (ξ * x)) volume (-1) 1 := by
    apply Continuous.intervalIntegrable; fun_prop
  simp only [h'] at hint
  rw [intervalIntegral.integral_add (hi1.const_mul _) (hi0.const_mul _),
    intervalIntegral.integral_const_mul,
    intervalIntegral.integral_const_mul] at hint
  unfold Icos
  linarith

/-- The closed form (3.2): `(sin ξ / ξ) ∏_{k=1}^m (kπ)² / ((kπ)² - ξ²)`. -/
noncomputable def kappaHat (m : ℕ) (ξ : ℝ) : ℝ :=
  sin ξ / ξ * ∏ k ∈ Finset.range m, ((k + 1 : ℝ) * π) ^ 2 / (((k + 1 : ℝ) * π) ^ 2 - ξ ^ 2)

/-- `ξ` is not a removable singularity of (3.2): `ξ ≠ 0` and `ξ² ≠ (kπ)²` for `1 ≤ k ≤ m`. -/
def Regular (m : ℕ) (ξ : ℝ) : Prop := ξ ≠ 0 ∧ ∀ k < m, ((k + 1 : ℝ) * π) ^ 2 - ξ ^ 2 ≠ 0

theorem Regular.mono {m n : ℕ} {ξ : ℝ} (h : Regular n ξ) (hmn : m ≤ n) : Regular m ξ :=
  ⟨h.1, fun k hk => h.2 k (lt_of_lt_of_le hk hmn)⟩

theorem Icos_zero {ξ : ℝ} (hξ : ξ ≠ 0) : Icos 0 ξ = 2 * sin ξ / ξ := by
  unfold Icos
  simp only [mul_zero, pow_zero, one_mul]
  rw [intervalIntegral.integral_comp_mul_left (fun x => cos x) hξ, integral_cos]
  simp only [mul_one, mul_neg, sin_neg, smul_eq_mul]
  field_simp
  ring

/-- `C_m I_m(ξ)` is the closed form (3.2). -/
theorem C_mul_Icos (m : ℕ) {ξ : ℝ} (h : Regular m ξ) : C m * Icos m ξ = kappaHat m ξ := by
  induction m with
  | zero => rw [Icos_zero h.1]; simp [C, kappaHat]; ring
  | succ m ih =>
    have hD := h.2 m (Nat.lt_succ_self m)
    have ih' := ih (h.mono (Nat.le_succ m))
    have hrec := Icos_recurrence m ξ
    have hI : Icos (m + 1) ξ =
        (m + 1) * (2 * m + 1) * π ^ 2 / 2 * Icos m ξ / (((m + 1 : ℝ) * π) ^ 2 - ξ ^ 2) := by
      rw [eq_div_iff hD]; linarith [hrec]
    have h21 : (2 * (m : ℝ) + 1) ≠ 0 := by positivity
    unfold kappaHat at ih'
    rw [C_succ, hI, kappaHat, Finset.prod_range_succ, ← mul_assoc (sin ξ / ξ), ← ih']
    field_simp

/-! ### From `I_m` to `κ_{2m}` -/

theorem kappa_mul_eq (m : ℕ) (g : ℝ → ℝ) (x : ℝ) :
    kappa m x * g x = indicator (Icc (-1) 1) (fun x => C m * cos (π * x / 2) ^ (2 * m) * g x) x := by
  by_cases hx : x ∈ Icc (-1 : ℝ) 1 <;> simp [kappa, hx]

theorem integral_kappa_mul (m : ℕ) (g : ℝ → ℝ) :
    ∫ x, kappa m x * g x = C m * ∫ x in (-1 : ℝ)..1, cos (π * x / 2) ^ (2 * m) * g x := by
  simp_rw [kappa_mul_eq]
  rw [integral_indicator measurableSet_Icc, integral_Icc_eq_integral_Ioc,
    ← intervalIntegral.integral_of_le (by norm_num), ← intervalIntegral.integral_const_mul]
  congr 1; funext x; ring

/-- The cosine transform of `κ_{2m}` is (3.2). -/
theorem integral_kappa_mul_cos (m : ℕ) {ξ : ℝ} (h : Regular m ξ) :
    ∫ x, kappa m x * cos (ξ * x) = kappaHat m ξ := by
  rw [integral_kappa_mul, ← C_mul_Icos m h, Icos]

/-- The sine transform of `κ_{2m}` vanishes (`κ_{2m}` is even). -/
theorem integral_kappa_mul_sin (m : ℕ) (ξ : ℝ) : ∫ x, kappa m x * sin (ξ * x) = 0 := by
  rw [integral_kappa_mul]
  set g : ℝ → ℝ := fun x => cos (π * x / 2) ^ (2 * m) * sin (ξ * x)
  have hodd : ∀ x, g (-x) = -g x := by
    intro x; simp only [g, mul_neg, neg_div, cos_neg, sin_neg]
  have h := intervalIntegral.integral_comp_neg (a := (-1 : ℝ)) (b := 1) g
  simp only [neg_neg, hodd, intervalIntegral.integral_neg] at h
  have : ∫ x in (-1 : ℝ)..1, g x = 0 := by linarith
  simp only [g] at this
  rw [this, mul_zero]

/-! ### The Fourier transform `f̂(γ) = ∫ f(y) e^{iγy} dy` -/

/-- `f̂(γ) = ∫ f(y) e^{iγy} dy`. -/
noncomputable def FT (f : ℝ → ℝ) (γ : ℝ) : ℂ :=
  ∫ y, (f y : ℂ) * Complex.exp (((γ * y : ℝ) : ℂ) * Complex.I)

theorem norm_exp_I (r : ℝ) : ‖Complex.exp ((r : ℂ) * Complex.I)‖ = 1 :=
  Complex.norm_exp_ofReal_mul_I r

theorem integrable_mul_exp {f : ℝ → ℝ} (hf : Integrable f) (γ : ℝ) :
    Integrable fun y => (f y : ℂ) * Complex.exp (((γ * y : ℝ) : ℂ) * Complex.I) := by
  have := hf.ofReal (𝕜 := ℂ) |>.bdd_mul (c := 1) (f := fun y =>
    Complex.exp (((γ * y : ℝ) : ℂ) * Complex.I)) (by fun_prop)
    (ae_of_all _ fun y => (norm_exp_I _).le)
  simpa [mul_comm] using this

theorem integrable_mul_trig {f : ℝ → ℝ} (hf : Integrable f) (t : ℝ → ℝ) (ht : Continuous t)
    (hb : ∀ y, |t y| ≤ 1) : Integrable fun y => f y * t y := by
  have := hf.bdd_mul (c := 1) ht.aestronglyMeasurable
    (ae_of_all _ fun y => by simpa [Real.norm_eq_abs] using hb y)
  simpa [mul_comm] using this

theorem FT_eq {f : ℝ → ℝ} (hf : Integrable f) (γ : ℝ) :
    FT f γ = ((∫ y, f y * cos (γ * y) : ℝ) : ℂ) + ((∫ y, f y * sin (γ * y) : ℝ) : ℂ) * Complex.I := by
  have hc := integrable_mul_trig hf (fun y => cos (γ * y)) (by fun_prop) fun y => abs_cos_le_one _
  have hs := integrable_mul_trig hf (fun y => sin (γ * y)) (by fun_prop) fun y => abs_sin_le_one _
  unfold FT
  have e : (fun y => (f y : ℂ) * Complex.exp (((γ * y : ℝ) : ℂ) * Complex.I)) =
      fun y => (((f y * cos (γ * y) : ℝ)) : ℂ) + (((f y * sin (γ * y) : ℝ)) : ℂ) * Complex.I := by
    funext y
    rw [Complex.exp_mul_I, ← Complex.ofReal_cos, ← Complex.ofReal_sin]
    push_cast; ring
  rw [e, integral_add hc.ofReal (hs.ofReal.mul_const _), integral_mul_const,
    integral_complex_ofReal, integral_complex_ofReal]

/-- **Lemma 3.1(c)**: `κ̂_{2m}(ξ) = (sin ξ / ξ) ∏_{k=1}^m (kπ)² / ((kπ)² - ξ²)`. -/
theorem FT_kappa (m : ℕ) {ξ : ℝ} (h : Regular m ξ) : FT (kappa m) ξ = (kappaHat m ξ : ℂ) := by
  rw [FT_eq (integrable_kappa m), integral_kappa_mul_cos m h, integral_kappa_mul_sin]
  simp

/-- `\widehat{f_η}(γ) = f̂(ηγ)` for `f_η(y) = η⁻¹ f(y/η)`. -/
theorem FT_scale {η : ℝ} (hη : 0 < η) (f : ℝ → ℝ) (γ : ℝ) : FT (scale η f) γ = FT f (η * γ) := by
  unfold FT scale
  set g : ℝ → ℂ := fun z => (f z : ℂ) * Complex.exp (((η * γ * z : ℝ) : ℂ) * Complex.I)
  have e : (fun y => (((η⁻¹ * f (y / η) : ℝ)) : ℂ) * Complex.exp (((γ * y : ℝ) : ℂ) * Complex.I)) =
      fun y => (η⁻¹ : ℂ) * g (y / η) := by
    funext y
    simp only [g]
    have : η * γ * (y / η) = γ * y := by field_simp
    rw [this]; push_cast; ring
  rw [e, integral_const_mul, Measure.integral_comp_div g η, abs_of_pos hη, Complex.real_smul,
    ← mul_assoc]
  rw [inv_mul_cancel₀ (by exact_mod_cast hη.ne'), one_mul]

/-- `\widehat{f ⋆ g} = f̂ ĝ`. -/
theorem FT_convolution {f g : ℝ → ℝ} (hf : Integrable f) (hg : Integrable g) (γ : ℝ) :
    FT (convolution f g (ContinuousLinearMap.lsmul ℝ ℝ) volume) γ = FT f γ * FT g γ := by
  set F : ℝ → ℂ := fun t => (f t : ℂ) * Complex.exp (((γ * t : ℝ) : ℂ) * Complex.I)
  set G : ℝ → ℂ := fun t => (g t : ℂ) * Complex.exp (((γ * t : ℝ) : ℂ) * Complex.I)
  have key : ∀ y, ((convolution f g (ContinuousLinearMap.lsmul ℝ ℝ) volume y : ℝ) : ℂ) *
      Complex.exp (((γ * y : ℝ) : ℂ) * Complex.I) =
      convolution F G (ContinuousLinearMap.mul ℂ ℂ) volume y := by
    intro y
    rw [convolution_def, convolution_def, ← integral_complex_ofReal, ← integral_mul_const]
    congr 1; funext t
    simp only [F, G, ContinuousLinearMap.lsmul_apply, ContinuousLinearMap.mul_apply',
      smul_eq_mul]
    rw [mul_mul_mul_comm, ← Complex.exp_add]
    push_cast
    ring_nf
  unfold FT
  simp_rw [key]
  rw [integral_convolution _ (integrable_mul_exp hf γ) (integrable_mul_exp hg γ)]
  rfl

theorem FT_add {f g : ℝ → ℝ} (hf : Integrable f) (hg : Integrable g) (γ : ℝ) :
    FT (fun y => f y + g y) γ = FT f γ + FT g γ := by
  unfold FT
  rw [← integral_add (integrable_mul_exp hf γ) (integrable_mul_exp hg γ)]
  congr 1; funext y; push_cast; ring

theorem FT_const_mul (c : ℝ) (f : ℝ → ℝ) (γ : ℝ) : FT (fun y => c * f y) γ = c * FT f γ := by
  unfold FT
  rw [← integral_const_mul]
  congr 1; funext y; push_cast; ring

/-- The transform of the core profile `λ₁κ₂ + λ₂κ₄ + λ₃κ₆`. -/
theorem FT_coreProfile (l₁ l₂ l₃ : ℝ) {ξ : ℝ} (h : Regular 3 ξ) :
    FT (coreProfile l₁ l₂ l₃) ξ =
      ((l₁ * kappaHat 1 ξ + l₂ * kappaHat 2 ξ + l₃ * kappaHat 3 ξ : ℝ) : ℂ) := by
  have hi : ∀ m (l : ℝ), Integrable fun y => l * kappa m y := fun m l =>
    (integrable_kappa m).const_mul l
  show FT (fun x => (l₁ * kappa 1 x + l₂ * kappa 2 x) + l₃ * kappa 3 x) ξ = _
  rw [FT_add (f := fun x => l₁ * kappa 1 x + l₂ * kappa 2 x) (g := fun x => l₃ * kappa 3 x)
      ((hi 1 l₁).add (hi 2 l₂)) (hi 3 l₃),
    FT_add (f := fun x => l₁ * kappa 1 x) (g := fun x => l₂ * kappa 2 x) (hi 1 l₁) (hi 2 l₂),
    FT_const_mul, FT_const_mul, FT_const_mul, FT_kappa 1 (h.mono (by norm_num)),
    FT_kappa 2 (h.mono (by norm_num)), FT_kappa 3 h]
  push_cast; ring

/-- **Eq. (3.4)**: `K̂(γ) = κ̂_core(η₁γ) κ̂_{2M}(δγ)`, with `κ̂_core = Σ λ_m κ̂_{2m}` and each
`κ̂_{2m}` given by the closed form (3.2). -/
theorem KParams.FT_K (P : KParams) (γ : ℝ) (hcore : Regular 3 (P.η₁ * γ))
    (hmoll : Regular P.M (P.δ * γ)) :
    FT P.K γ = (((P.l₁ * kappaHat 1 (P.η₁ * γ) + P.l₂ * kappaHat 2 (P.η₁ * γ) +
      P.l₃ * kappaHat 3 (P.η₁ * γ)) * kappaHat P.M (P.δ * γ) : ℝ) : ℂ) := by
  unfold KParams.K
  rw [FT_convolution P.integrable_core P.integrable_moll, KParams.core, KParams.moll,
    FT_scale P.hη₁, FT_scale P.hδ, FT_coreProfile _ _ _ hcore, FT_kappa _ hmoll]
  push_cast; ring

/-- **Lemma 3.1(c)**, first bound: `|κ̂_{2m}(ξ)| ≤ 1`. -/
theorem abs_kappaHat_le_one (m : ℕ) {ξ : ℝ} (h : Regular m ξ) : |kappaHat m ξ| ≤ 1 := by
  have hint := integrable_mul_trig (integrable_kappa m) (fun x => cos (ξ * x)) (by fun_prop)
    fun x => abs_cos_le_one _
  rw [← integral_kappa_mul_cos m h, ← integral_kappa m]
  refine (abs_integral_le_integral_abs).trans (integral_mono hint.abs (integrable_kappa m) ?_)
  intro x
  simp only
  rw [abs_mul, abs_of_nonneg (kappa_nonneg m x)]
  exact mul_le_of_le_one_right (kappa_nonneg m x) (abs_cos_le_one _)

/-- **Lemma 3.1(c)**, second bound: for `ξ > mπ`,
`|κ̂_{2m}(ξ)| ≤ Φ_{2m}(ξ) = ξ⁻¹ ∏_{k=1}^m (kπ)² / (ξ² - (kπ)²)`. -/
theorem abs_kappaHat_le_envelope (m : ℕ) {ξ : ℝ} (hξ : (m : ℝ) * π < ξ) :
    |kappaHat m ξ| ≤ 1 / ξ * ∏ k ∈ Finset.range m, ((k + 1 : ℝ) * π) ^ 2 / (ξ ^ 2 - ((k + 1 : ℝ) * π) ^ 2) := by
  have hξ0 : 0 < ξ := lt_of_le_of_lt (by positivity) hξ
  have hk : ∀ k ∈ Finset.range m, ((k + 1 : ℝ) * π) ^ 2 < ξ ^ 2 := by
    intro k hk
    have : (k + 1 : ℝ) ≤ m := by exact_mod_cast Finset.mem_range.mp hk
    have h1 : 0 < (k + 1 : ℝ) * π := by positivity
    have h2 : (k + 1 : ℝ) * π < ξ := lt_of_le_of_lt (by nlinarith [pi_pos]) hξ
    nlinarith
  unfold kappaHat
  rw [abs_mul, Finset.abs_prod, abs_div, abs_of_pos hξ0]
  have hprod : ∏ k ∈ Finset.range m, |((k + 1 : ℝ) * π) ^ 2 / (((k + 1 : ℝ) * π) ^ 2 - ξ ^ 2)| =
      ∏ k ∈ Finset.range m, ((k + 1 : ℝ) * π) ^ 2 / (ξ ^ 2 - ((k + 1 : ℝ) * π) ^ 2) := by
    refine Finset.prod_congr rfl fun k hkm => ?_
    have := hk k hkm
    rw [abs_div, abs_of_pos (by positivity), abs_of_neg (by linarith), neg_sub]
  rw [hprod]
  refine mul_le_mul_of_nonneg_right ?_ (Finset.prod_nonneg fun k hkm => ?_)
  · exact div_le_div_of_nonneg_right (abs_sin_le_one _) hξ0.le
  · have := hk k hkm
    exact div_nonneg (by positivity) (by linarith)

end SkewesVerify
