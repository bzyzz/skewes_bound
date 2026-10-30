import SkewesVerify.Remainders
import SkewesVerify.Decimal

/-!
# Theorems 1.4 and 1.5: Mertens' first theorem

Write `F(u) = e^{u/2} Δ₁(e^u)`, where `Δ₁(x) = Σ_{p ≤ x} log p / p - log x - E` (`Defs.lean`).
The hypotheses are:

* `hExpl` **[cited + derived]** Proposition 5.3 (Lemmas 5.1–5.2, the explicit formula,
  Platt–Trudgian's `|θ(y) - y| ≤ ε y`, Rosser–Schoenfeld):
  `F(u) ≤ 1 + 3ε + Z(u) + r(u)` for almost every `u` in the window, with `ε = 1.5423·10⁻⁹`,
  the zero sum `Z(u) = -Σ_ρ e^{(ρ-1/2)u} / (ρ - 1)` and `r = mertensR` (`Remainders.lean`);
* `hZint`, `hZ` **[cited + derived]** the proof of Theorem 5.5 (Platt–Trudgian's RH verification
  to `3·10¹²`, Lemmas 3.2–3.3, and Lemma 2.9): `∫ K(u - ω) Z(u) du ≤ S + R₄ + R₅`;
* `hS`, `hR4`, `hR5` **[Arb]** `mertens/lmfdb/first_cert_lmfdb_T1000000000.json` (Theorem 1.4)
  and `mertens/criterion/first_cert_T1300000.json` (Theorem 1.5). `hS` is the upper end of the
  ball `S_M`. `R4_stieltjes` and `R5` are printed to 10 and 6 digits, so the hypotheses raise them
  by one unit in the last printed digit.

Lean bounds `r(u)` on the window itself (`mertensR_le`), so the certificates' `R_small` is not
needed.
-/

open Real MeasureTheory Set

namespace SkewesVerify.Mertens

/-- `F(u) = e^{u/2} Δ₁(e^u)`. -/
noncomputable def F (u : ℝ) : ℝ := exp (u / 2) * Delta1 (exp u)

/-- `1 + 3ε` with `ε = 1.5423·10⁻⁹`. -/
noncomputable def bias : ℝ := 1 + 3 * (15423 / 10 ^ 13)

/-- The common step of Theorems 1.4 and 1.5. -/
theorem criterion (P : KParams) (ω : ℝ) (hlo : 0 ≤ ω - P.η₀) (Z : ℝ → ℝ) (B : ℝ)
    (hExpl : ∀ᵐ u, u ∈ P.window ω → F u ≤ bias + Z u + mertensR u)
    (hZint : IntegrableOn (fun u => P.K (u - ω) * Z u) (P.window ω))
    (hZ : ∫ u in P.window ω, P.K (u - ω) * Z u ≤ B)
    (hneg : bias + B + (5 + 136 / 100 * ((ω + P.η₀) / (69 / 100) + 2)) *
      (1 / (1 + (ω - P.η₀) / 120) ^ 20) < 0) :
    ∃ u ∈ P.window ω, F u < 0 := by
  set c := (5 + 136 / 100 * ((ω + P.η₀) / (69 / 100) + 2)) * (1 / (1 + (ω - P.η₀) / 120) ^ 20)
  have hr := P.integrableOn_K_mul (ω := ω) (continuousOn_mertensR (s := P.window ω))
  have hrle : ∫ u in P.window ω, P.K (u - ω) * mertensR u ≤ c :=
    P.integral_K_le_of_le continuousOn_mertensR fun u hu => mertensR_le hlo hu
  have hint : ∫ u in P.window ω, P.K (u - ω) * ((Z u + mertensR u) + bias) =
      (∫ u in P.window ω, P.K (u - ω) * Z u) + (∫ u in P.window ω, P.K (u - ω) * mertensR u) +
        bias := by
    rw [P.integral_K_add (P.integrableOn_K_add hZint hr) (P.integrableOn_K_const ω bias),
      P.integral_K_add hZint hr, P.integral_K_const]
  refine exists_neg_of_integral_neg (g := fun u => (Z u + mertensR u) + bias)
    (P.measurableSet_window ω) (fun u _ => P.K_shift_nonneg ω u)
    (hExpl.mono fun u h hu => by linarith [h hu]) ?_
  rw [hint]
  linarith

theorem Delta1_neg {u : ℝ} (h : F u < 0) : Delta1 (exp u) < 0 := by
  unfold F at h
  by_contra h'
  exact absurd h (not_lt.mpr (mul_nonneg (exp_pos _).le (not_lt.mp h')))

/-! ### Theorem 1.4 (`T₀ = 10⁹`) -/

/-- `η₁ = 1.4·10⁻⁸`, `δ = 5·10⁻¹⁰`, `M = 135`, `λ = (0.04, 0.52, 0.44)`. -/
noncomputable def P215 : KParams where
  l₁ := 4 / 100
  l₂ := 52 / 100
  l₃ := 44 / 100
  η₁ := 14 / 10 ^ 9
  δ := 5 / 10 ^ 10
  M := 135
  hl₁ := by norm_num
  hl₂ := by norm_num
  hl₃ := by norm_num
  hsum := by norm_num
  hη₁ := by norm_num
  hδ := by norm_num

/-- `ω = 495.7028331322`. -/
noncomputable def ω215 : ℝ := 4957028331322 / 10 ^ 10

theorem window215 :
    P215.window ω215 = Icc (4957028331177 / 10 ^ 10) (4957028331467 / 10 ^ 10) := by
  simp only [KParams.window, KParams.η₀, P215, ω215]; norm_num

/-- **Theorem 1.4.** -/
theorem mertens215 (Z : ℝ → ℝ) (S R₄ R₅ : ℝ)
    (hExpl : ∀ᵐ u, u ∈ P215.window ω215 → F u ≤ bias + Z u + mertensR u)
    (hZint : IntegrableOn (fun u => P215.K (u - ω215) * Z u) (P215.window ω215))
    (hZ : ∫ u in P215.window ω215, P215.K (u - ω215) * Z u ≤ S + R₄ + R₅)
    (hS : S ≤ -100024 / 10 ^ 5 + 573 / 10 ^ 8)
    (hR4 : R₄ ≤ 4786914752 / 10 ^ 14)
    (hR5 : R₅ ≤ 389001 / 10 ^ 43) :
    ∃ x : ℝ, log x ∈ Icc (4957028331177 / 10 ^ 10) (4957028331467 / 10 ^ 10) ∧
      (∑ p ∈ (Finset.range (⌊x⌋₊ + 1)).filter Nat.Prime, log (p : ℝ) / p) < log x + mertensE ∧
      x < (10 : ℝ) ^ (2152811 / 10000 : ℝ) := by
  have hlo : 0 ≤ ω215 - P215.η₀ := by simp only [KParams.η₀, P215, ω215]; norm_num
  obtain ⟨u, hu, hneg⟩ := criterion P215 ω215 hlo Z _ hExpl hZint hZ (by
    simp only [KParams.η₀, P215, ω215, bias]
    norm_num
    linarith)
  rw [window215] at hu
  refine ⟨exp u, by rw [log_exp]; exact hu, ?_, ?_⟩
  · have := Delta1_neg hneg; unfold Delta1 at this; linarith
  · exact exp_lt_ten_rpow (by norm_num) (by unfold L10; linarith [hu.2])

/-! ### Theorem 1.5 (`T₀ = 1.3·10⁶`) -/

/-- `η₁ = 9.24·10⁻⁶`, `δ = 2·10⁻⁹`, `M = 494`, `λ = (0.04, 0.52, 0.44)`. -/
noncomputable def P807 : KParams where
  l₁ := 4 / 100
  l₂ := 52 / 100
  l₃ := 44 / 100
  η₁ := 924 / 10 ^ 8
  δ := 2 / 10 ^ 9
  M := 494
  hl₁ := by norm_num
  hl₂ := by norm_num
  hl₃ := by norm_num
  hsum := by norm_num
  hη₁ := by norm_num
  hδ := by norm_num

/-- `ω = 1859.1290514`. -/
noncomputable def ω807 : ℝ := 18591290514 / 10 ^ 7

theorem window807 :
    P807.window ω807 = Icc (1859129042158 / 10 ^ 9) (1859129060642 / 10 ^ 9) := by
  simp only [KParams.window, KParams.η₀, P807, ω807]; norm_num

/-- **Theorem 1.5.** -/
theorem mertens807 (Z : ℝ → ℝ) (S R₄ R₅ : ℝ)
    (hExpl : ∀ᵐ u, u ∈ P807.window ω807 → F u ≤ bias + Z u + mertensR u)
    (hZint : IntegrableOn (fun u => P807.K (u - ω807) * Z u) (P807.window ω807))
    (hZ : ∫ u in P807.window ω807, P807.K (u - ω807) * Z u ≤ S + R₄ + R₅)
    (hS : S ≤ -101044645655603 / 10 ^ 14 + 822 / 10 ^ 17)
    (hR4 : R₄ ≤ 7495480210 / 10 ^ 14)
    (hR5 : R₅ ≤ 680017 / 10 ^ 179) :
    ∃ x : ℝ, log x ∈ Icc (1859129042158 / 10 ^ 9) (1859129060642 / 10 ^ 9) ∧ Delta1 x < 0 ∧
      x < (10 : ℝ) ^ (8074095 / 10000 : ℝ) := by
  have hlo : 0 ≤ ω807 - P807.η₀ := by simp only [KParams.η₀, P807, ω807]; norm_num
  obtain ⟨u, hu, hneg⟩ := criterion P807 ω807 hlo Z _ hExpl hZint hZ (by
    simp only [KParams.η₀, P807, ω807, bias]
    norm_num
    linarith)
  rw [window807] at hu
  exact ⟨exp u, by rw [log_exp]; exact hu, Delta1_neg hneg,
    exp_lt_ten_rpow (by norm_num) (by unfold L10; linarith [hu.2])⟩

end SkewesVerify.Mertens
