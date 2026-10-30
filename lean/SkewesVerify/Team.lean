import SkewesVerify.Averaging
import SkewesVerify.Decimal

/-!
# Theorem 1.3: residues against non-residues modulo 8 and 12

For `q ∈ {8, 12}`, the theorem says that `3 π(x; q, 1) > Σ_{b ≠ 1} π(x; q, b)` for some `x` in an
explicit window. Write `G(u) = u e^{-u/2} (3 π(e^u; q, 1) - Σ_{b ≠ 1} π(e^u; q, b))`.
The hypotheses are:

* `hExpl` **[cited + derived]** §4.5, "the explicit inequality":
  `G(u) ≥ -3 Φ(u) - Z(u) - ε₀` for almost every `u` in the window, where `Z(u)` is the zero sum
  `Σ_χ Σ_ρ u e^{-u/2} li(e^{ρu})` and `ε₀` bounds `ε₀^{(q)}(u)` on the window;
* `hZint`, `hZ` **[cited + derived]** §4.5, "the criterion" (the argument of Theorem 4.7):
  `∫ K(u - ω) Z(u) du ≤ S + R₃ + R₄ + R₅`;
* `hS`, `hR3`, `hR4`, `hR5`, `hε` **[Arb]** `race/team/team_cert_q8_T100000_w64.764612.json` and
  `race/team/team_cert_q12_T50000_w339.632308.json` (upper ends of the balls `sum_S`, `R3`, `R4`,
  `R5`, `eps0`).
-/

open Real MeasureTheory Set

namespace SkewesVerify.Team

/-- `G(u)` for modulus `q` with non-residue classes `b₁, b₂, b₃`. -/
noncomputable def G (q b₁ b₂ b₃ : ℕ) (u : ℝ) : ℝ :=
  u * exp (-(u / 2)) * (3 * piMod q 1 (exp u) -
    (piMod q b₁ (exp u) + piMod q b₂ (exp u) + piMod q b₃ (exp u)))

/-- The common step: if `∫ K (-3Φ - Z - ε₀) > 0`, then `G > 0` somewhere on the window. -/
theorem criterion (P : KParams) (ω : ℝ) (hlo : 0 < ω - P.η₀) (q b₁ b₂ b₃ : ℕ)
    (Z : ℝ → ℝ) (ε₀ B : ℝ)
    (hExpl : ∀ᵐ u, u ∈ P.window ω → -3 * Phi u - Z u - ε₀ ≤ G q b₁ b₂ b₃ u)
    (hZint : IntegrableOn (fun u => P.K (u - ω) * Z u) (P.window ω))
    (hZ : ∫ u in P.window ω, P.K (u - ω) * Z u ≤ B)
    (hpos : 0 < -3 * Phi (ω - P.η₀) - B - ε₀) :
    ∃ u ∈ P.window ω, 0 < G q b₁ b₂ b₃ u := by
  have hΦ := P.integrableOn_K_mul (P.continuousOn_Phi_window hlo)
  have h3 := P.integrableOn_K_const_mul 3 hΦ
  have hint : ∫ u in P.window ω, P.K (u - ω) * (-3 * Phi u - Z u - ε₀) =
      -((∫ u in P.window ω, P.K (u - ω) * (3 * Phi u + Z u + ε₀))) := by
    rw [← integral_neg]; congr 1; ext u; ring
  rw [P.integral_K_add (P.integrableOn_K_add h3 hZint) (P.integrableOn_K_const ω ε₀),
    P.integral_K_add h3 hZint, P.integral_K_const, P.integral_K_const_mul] at hint
  refine exists_pos_of_integral_pos (P.measurableSet_window ω) (fun u _ => P.K_shift_nonneg ω u)
    hExpl ?_
  rw [hint]
  linarith [P.integral_K_Phi_le hlo]

/-! ### Modulo 8 -/

/-- The kernel for `q = 8`: `κ₆` core, `η₁ = 1.4·10⁻⁴`, `δ = 2·10⁻⁵`, `M = 40`. -/
noncomputable def P8 : KParams where
  l₁ := 0
  l₂ := 0
  l₃ := 1
  η₁ := 14 / 10 ^ 5
  δ := 2 / 10 ^ 5
  M := 40
  hl₁ := le_rfl
  hl₂ := le_rfl
  hl₃ := zero_le_one
  hsum := by norm_num
  hη₁ := by norm_num
  hδ := by norm_num

/-- `ω = 64.764612`. -/
noncomputable def ω8 : ℝ := 64764612 / 10 ^ 6

theorem window8 : P8.window ω8 = Icc (64764452 / 10 ^ 6) (64764772 / 10 ^ 6) := by
  simp only [KParams.window, KParams.η₀, P8, ω8]; norm_num

/-- **Theorem 1.3(a).** -/
theorem team8 (Z : ℝ → ℝ) (ε₀ S R₃ R₄ R₅ : ℝ)
    (hExpl : ∀ᵐ u, u ∈ P8.window ω8 → -3 * Phi u - Z u - ε₀ ≤ G 8 3 5 7 u)
    (hZint : IntegrableOn (fun u => P8.K (u - ω8) * Z u) (P8.window ω8))
    (hZ : ∫ u in P8.window ω8, P8.K (u - ω8) * Z u ≤ S + R₃ + R₄ + R₅)
    (hS : S ≤ -3150 / 1000 + 404 / 10 ^ 6)
    (hR3 : R₃ ≤ 4570310812 / 10 ^ 14 + 925 / 10 ^ 18)
    (hR4 : R₄ ≤ 1020613101 / 10 ^ 12 + 281 / 10 ^ 15)
    (hR5 : R₅ ≤ 3518406877 / 10 ^ 17 + 474 / 10 ^ 20)
    (hε : ε₀ ≤ 1067077596 / 10 ^ 13 + 246 / 10 ^ 16) :
    ∃ x : ℝ, log x ∈ Icc (64764452 / 10 ^ 6) (64764772 / 10 ^ 6) ∧
      (piMod 8 3 x + piMod 8 5 x + piMod 8 7 x) / 3 < piMod 8 1 x ∧
      x < (10 : ℝ) ^ (281270 / 10000 : ℝ) := by
  have hlo : 0 < ω8 - P8.η₀ := by simp only [KParams.η₀, P8, ω8]; norm_num
  have hΦ : 3 * Phi (ω8 - P8.η₀) < 3099011998 / 10 ^ 9 := by
    simp only [KParams.η₀, P8, ω8, Phi]; norm_num
  obtain ⟨u, hu, hpos⟩ := criterion P8 ω8 hlo 8 3 5 7 Z ε₀ _ hExpl hZint hZ (by linarith)
  rw [window8] at hu
  refine ⟨exp u, by rw [log_exp]; exact hu, ?_, ?_⟩
  · have hu0 : 0 < u * exp (-(u / 2)) := mul_pos (lt_of_lt_of_le (by norm_num) hu.1) (exp_pos _)
    unfold G at hpos
    have := (pos_iff_pos_of_mul_pos hpos).mp hu0
    linarith
  · exact exp_lt_ten_rpow (by norm_num) (by unfold L10; linarith [hu.2])

/-! ### Modulo 12 -/

/-- The kernel for `q = 12`: `κ₆` core, `η₁ = 2.8·10⁻⁴`, `δ = 1.1·10⁻⁴`, `M = 60`. -/
noncomputable def P12 : KParams where
  l₁ := 0
  l₂ := 0
  l₃ := 1
  η₁ := 28 / 10 ^ 5
  δ := 11 / 10 ^ 5
  M := 60
  hl₁ := le_rfl
  hl₂ := le_rfl
  hl₃ := zero_le_one
  hsum := by norm_num
  hη₁ := by norm_num
  hδ := by norm_num

/-- `ω = 339.632308`. -/
noncomputable def ω12 : ℝ := 339632308 / 10 ^ 6

theorem window12 : P12.window ω12 = Icc (339631918 / 10 ^ 6) (339632698 / 10 ^ 6) := by
  simp only [KParams.window, KParams.η₀, P12, ω12]; norm_num

/-- **Theorem 1.3(b).** -/
theorem team12 (Z : ℝ → ℝ) (ε₀ S R₃ R₄ R₅ : ℝ)
    (hExpl : ∀ᵐ u, u ∈ P12.window ω12 → -3 * Phi u - Z u - ε₀ ≤ G 12 5 7 11 u)
    (hZint : IntegrableOn (fun u => P12.K (u - ω12) * Z u) (P12.window ω12))
    (hZ : ∫ u in P12.window ω12, P12.K (u - ω12) * Z u ≤ S + R₃ + R₄ + R₅)
    (hS : S ≤ -3147 / 1000 + 639 / 10 ^ 6)
    (hR3 : R₃ ≤ 1255923169 / 10 ^ 15 + 128 / 10 ^ 18)
    (hR4 : R₄ ≤ 9440472439 / 10 ^ 13 + 611 / 10 ^ 17)
    (hR5 : R₅ ≤ 4963374721 / 10 ^ 16 + 723 / 10 ^ 21)
    (hε : ε₀ ≤ 1357142235 / 10 ^ 33 + 479 / 10 ^ 36) :
    ∃ x : ℝ, log x ∈ Icc (339631918 / 10 ^ 6) (339632698 / 10 ^ 6) ∧
      (piMod 12 5 x + piMod 12 7 x + piMod 12 11 x) / 3 < piMod 12 1 x ∧
      x < (10 : ℝ) ^ (1475007 / 10000 : ℝ) := by
  have hlo : 0 < ω12 - P12.η₀ := by simp only [KParams.η₀, P12, ω12]; norm_num
  have hΦ : 3 * Phi (ω12 - P12.η₀) < 3017878732 / 10 ^ 9 := by
    simp only [KParams.η₀, P12, ω12, Phi]; norm_num
  obtain ⟨u, hu, hpos⟩ := criterion P12 ω12 hlo 12 5 7 11 Z ε₀ _ hExpl hZint hZ (by linarith)
  rw [window12] at hu
  refine ⟨exp u, by rw [log_exp]; exact hu, ?_, ?_⟩
  · have hu0 : 0 < u * exp (-(u / 2)) := mul_pos (lt_of_lt_of_le (by norm_num) hu.1) (exp_pos _)
    unfold G at hpos
    have := (pos_iff_pos_of_mul_pos hpos).mp hu0
    linarith
  · exact exp_lt_ten_rpow (by norm_num) (by unfold L10; linarith [hu.2])

end SkewesVerify.Team
