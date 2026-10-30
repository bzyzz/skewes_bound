import SkewesVerify.Averaging
import SkewesVerify.Decimal

/-!
# Theorem 1.1 and Corollary 1.2: the race `π(x; 24, 11)` against `π(x; 24, 1)`

Parameters (§4.4): pure `κ₆` core, `η₁ = 1.2·10⁻⁴`, `δ = 7.62·10⁻⁵`, `M = 40`, `ω = 111.54583`, so
`η₀ = 1.962·10⁻⁴` and the window is `[111.5456338, 111.5460262]`.

Write `E(u) = 8 u e^{-u/2} (π(e^u; 24, 11) - π(e^u; 24, 1))`. The hypotheses are:

* `hExpl` **[cited + derived]** Lemma 4.3 (Landau's explicit formula, orthogonality, Dusart):
  `E(u) ≤ 8 Φ(u) + Z(u) + ε₀` for almost every `u` in the window, where `Z(u)` is the zero sum
  `2 Σ_χ Σ_ρ u e^{-u/2} li(e^{ρu})` and `ε₀` bounds `ε₀(u)` on the window;
* `hZint`, `hZ` **[cited + derived]** the proof of Theorem 4.7 (term-by-term integration, GRH to
  height `A_χ` (Platt), BMOR zero counts, Lemma 3.2): `∫ K(u - ω) Z(u) du ≤ S + R₃ + R₄ + R₅`;
* `hS`, `hR3`, `hR4`, `hR5`, `hε` **[Arb]** `race/criterion/cert_T100000_w111.5458300.json`:
  `two_sum_S = [-8.165 ± 7.98e-4]`, `R3 = [5.779762829e-5 ± 2.81e-15]`,
  `R4 = [0.01130541581 ± 3.22e-12]`, `R5 = [9.129669155e-6 ± 2.62e-16]`,
  `eps0 = [8.771848710e-8 ± 4.55e-18]`. Each hypothesis is the upper end of the ball.

Lean proves the rest: the kernel facts, `∫ K Φ ≤ Φ(u_lo)`, the exact value of `8 Φ(u_lo)`, the
sign of `U`, the averaging step and `e^{111.5460262} < 10^{48.4439}`.
-/

open Real MeasureTheory Set

namespace SkewesVerify.Race

/-- The kernel of Theorem 1.1. -/
noncomputable def P : KParams where
  l₁ := 0
  l₂ := 0
  l₃ := 1
  η₁ := 12 / 10 ^ 5
  δ := 762 / 10 ^ 7
  M := 40
  hl₁ := le_rfl
  hl₂ := le_rfl
  hl₃ := zero_le_one
  hsum := by norm_num
  hη₁ := by norm_num
  hδ := by norm_num

/-- `ω = 111.54583`. -/
noncomputable def ω : ℝ := 11154583 / 10 ^ 5

theorem η₀_eq : P.η₀ = 1962 / 10 ^ 7 := by simp only [KParams.η₀, P]; norm_num

theorem window_eq : P.window ω = Icc (1115456338 / 10 ^ 7) (1115460262 / 10 ^ 7) := by
  simp only [KParams.window, η₀_eq, ω]; norm_num

/-- `E(u) = 8 u e^{-u/2} Δ_{24,11,1}(e^u)`. -/
noncomputable def E (u : ℝ) : ℝ := 8 * u * exp (-(u / 2)) * (piMod 24 11 (exp u) - piMod 24 1 (exp u))

theorem integral_le (Z : ℝ → ℝ) (ε₀ S R₃ R₄ R₅ : ℝ)
    (hZint : IntegrableOn (fun u => P.K (u - ω) * Z u) (P.window ω))
    (hZ : ∫ u in P.window ω, P.K (u - ω) * Z u ≤ S + R₃ + R₄ + R₅) :
    ∫ u in P.window ω, P.K (u - ω) * (8 * Phi u + Z u + ε₀) ≤
      8 * Phi (ω - P.η₀) + (S + R₃ + R₄ + R₅) + ε₀ := by
  have hlo : 0 < ω - P.η₀ := by rw [η₀_eq, ω]; norm_num
  have hΦ := P.integrableOn_K_mul (P.continuousOn_Phi_window hlo)
  have h8 := P.integrableOn_K_const_mul 8 hΦ
  rw [P.integral_K_add (P.integrableOn_K_add h8 hZint) (P.integrableOn_K_const ω ε₀),
    P.integral_K_add h8 hZint, P.integral_K_const, P.integral_K_const_mul]
  linarith [P.integral_K_Phi_le hlo]

/-- **Theorem 1.1.** -/
theorem race (Z : ℝ → ℝ) (ε₀ S R₃ R₄ R₅ : ℝ)
    (hExpl : ∀ᵐ u, u ∈ P.window ω → E u ≤ 8 * Phi u + Z u + ε₀)
    (hZint : IntegrableOn (fun u => P.K (u - ω) * Z u) (P.window ω))
    (hZ : ∫ u in P.window ω, P.K (u - ω) * Z u ≤ S + R₃ + R₄ + R₅)
    (hS : S ≤ -8165 / 1000 + 798 / 10 ^ 6)
    (hR3 : R₃ ≤ 5779762829 / 10 ^ 14 + 281 / 10 ^ 17)
    (hR4 : R₄ ≤ 1130541581 / 10 ^ 11 + 322 / 10 ^ 14)
    (hR5 : R₅ ≤ 9129669155 / 10 ^ 15 + 262 / 10 ^ 18)
    (hε : ε₀ ≤ 8771848710 / 10 ^ 17 + 455 / 10 ^ 20) :
    ∃ x : ℝ, log x ∈ Icc (1115456338 / 10 ^ 7) (1115460262 / 10 ^ 7) ∧
      piMod 24 11 x < piMod 24 1 x ∧ x < (10 : ℝ) ^ (484439 / 10000 : ℝ) := by
  have hU := integral_le Z ε₀ S R₃ R₄ R₅ hZint hZ
  have hΦ : 8 * Phi (ω - P.η₀) < 8148920290 / 10 ^ 9 := by
    rw [η₀_eq, ω, Phi]; norm_num
  obtain ⟨u, hu, hneg⟩ := exists_neg_of_integral_neg (P.measurableSet_window ω)
    (fun u _ => P.K_shift_nonneg ω u) hExpl (by linarith)
  rw [window_eq] at hu
  refine ⟨exp u, by rw [log_exp]; exact hu, ?_, ?_⟩
  · have hpos : 0 < 8 * u * exp (-(u / 2)) := by
      have : 0 < u := lt_of_lt_of_le (by norm_num) hu.1
      positivity
    unfold E at hneg
    by_contra h
    exact absurd hneg (not_lt.mpr (mul_nonneg hpos.le (by linarith)))
  · exact exp_lt_ten_rpow (by norm_num) (by unfold L10; linarith [hu.2])

/-- **Corollary 1.2.** For `b ≠ 11` the windows are those of Ford–Hudson [FH01, Table 3] (cited as
`hFH`); the largest ends at `161.22837`. -/
theorem race24 (hFH : ∀ b ∈ ({5, 7, 13, 17, 19, 23} : Finset ℕ),
      ∃ x : ℝ, log x ≤ 16122837 / 10 ^ 5 ∧ piMod 24 b x < piMod 24 1 x)
    (h11 : ∃ x : ℝ, log x ∈ Icc (1115456338 / 10 ^ 7) (1115460262 / 10 ^ 7) ∧
      piMod 24 11 x < piMod 24 1 x) :
    ∀ b ∈ ({5, 7, 11, 13, 17, 19, 23} : Finset ℕ),
      ∃ x : ℝ, x < (10 : ℝ) ^ (70021 / 1000 : ℝ) ∧ piMod 24 b x < piMod 24 1 x := by
  have key : ∀ x : ℝ, log x ≤ 16122837 / 10 ^ 5 → x < (10 : ℝ) ^ (70021 / 1000 : ℝ) := by
    intro x hx
    rcases le_or_gt x 0 with h | h
    · exact lt_of_le_of_lt h (by positivity)
    · rw [← exp_log h]; exact exp_lt_ten_rpow (by norm_num) (by unfold L10; linarith)
  intro b hb
  by_cases h : b = 11
  · subst h
    obtain ⟨x, hx, hlt⟩ := h11
    exact ⟨x, key x (by linarith [hx.2]), hlt⟩
  · have hb' : b ∈ ({5, 7, 13, 17, 19, 23} : Finset ℕ) := by
      simp only [Finset.mem_insert, Finset.mem_singleton] at hb ⊢; omega
    obtain ⟨x, hx, hlt⟩ := hFH b hb'
    exact ⟨x, key x hx, hlt⟩

end SkewesVerify.Race
