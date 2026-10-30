import SkewesVerify.Remainders
import SkewesVerify.Decimal

/-!
# Theorem 1.6: Skewes's number

Write `f(u) = u (π(e^u) - li(e^u)) / (e^{u/2} Φ(u))`, the integrand of eq. (6.1).
The hypotheses are:

* `hExpl` **[cited + derived]** Lemma 6.1 (the Riemann–von Mangoldt explicit formula, Dusart):
  `f(u) ≥ -1 - Z(u) - rem(u)` on the window, where `Z(u) = Σ_ρ u li(e^{ρu}) / (e^{u/2} Φ(u))`
  and `rem = skewesRem` (`Remainders.lean`);
* `hZint`, `hZ` **[cited + derived]** the proof of Theorem 6.2 (Platt–Trudgian's RH verification
  to `3·10¹²`, Lemmas 2.1, 2.9, 3.2, 3.3):
  `∫ K(u - ω) Z(u) du ≤ -S_K / Φ(ω) + R₃ + R₄ + R₅`;
* `hS`, `hR3`, `hR4`, `hR5` **[Arb]** `certification/s_certificate_arb.json`
  (`S_K = [1.00292 ± 2.27e-6]`, lower end) and `certification/error_bounds.py` at
  `ω = 727.951335444`, `T₀ = 8·10⁷` (`R₃ = 5.5408…e-9`, `R₄,env = 1.32497…e-4`,
  `R₅ = 8.5026…e-20`). The hypotheses use the values of Table 6 (`5.55e-9`, `1.3250e-4`,
  `8.51e-20`), which are larger.

Lean computes `Φ(ω)` exactly and bounds the remainder `rem` itself (`skewesRem_le`), so the
paper's `R₂` is not needed. Lean does not prove the `2.75·10¹⁵³` consecutive integers, which
rely on [SD10, Theorem 6.2].
-/

open Real MeasureTheory Set

namespace SkewesVerify.Skewes

/-- `f(u) = u (π(e^u) - li(e^u)) / (e^{u/2} Φ(u))`. -/
noncomputable def f (u : ℝ) : ℝ := u * (piR (exp u) - li (exp u)) / (exp (u / 2) * Phi u)

/-- The kernel: `η₁ = 1.5·10⁻⁷`, `δ = 2·10⁻⁹`, `M = 60`, `λ = (0.04, 0.52, 0.44)`. -/
noncomputable def P : KParams where
  l₁ := 4 / 100
  l₂ := 52 / 100
  l₃ := 44 / 100
  η₁ := 15 / 10 ^ 8
  δ := 2 / 10 ^ 9
  M := 60
  hl₁ := by norm_num
  hl₂ := by norm_num
  hl₃ := by norm_num
  hsum := by norm_num
  hη₁ := by norm_num
  hδ := by norm_num

/-- `ω = 727.951335444`. -/
noncomputable def ω : ℝ := 727951335444 / 10 ^ 9

theorem window_eq : P.window ω = Icc (727951335292 / 10 ^ 9) (727951335596 / 10 ^ 9) := by
  simp only [KParams.window, KParams.η₀, P, ω]; norm_num

/-- **Theorem 1.6.** -/
theorem skewes (Z : ℝ → ℝ) (S R₃ R₄ R₅ : ℝ)
    (hExpl : ∀ᵐ u, u ∈ P.window ω → -1 - Z u - skewesRem u ≤ f u)
    (hZint : IntegrableOn (fun u => P.K (u - ω) * Z u) (P.window ω))
    (hZ : ∫ u in P.window ω, P.K (u - ω) * Z u ≤ -S / Phi ω + R₃ + R₄ + R₅)
    (hS : 100292 / 10 ^ 5 - 227 / 10 ^ 8 ≤ S)
    (hR3 : R₃ ≤ 555 / 10 ^ 11)
    (hR4 : R₄ ≤ 13250 / 10 ^ 8)
    (hR5 : R₅ ≤ 851 / 10 ^ 22) :
    ∃ x : ℝ, log x ∈ Icc (727951335292 / 10 ^ 9) (727951335596 / 10 ^ 9) ∧ li x < piR x ∧
      2 ≤ x ∧ x < 1397166434 / 10 ^ 9 * 10 ^ 316 := by
  have hlo : 0 < ω - P.η₀ := by simp only [KParams.η₀, P, ω]; norm_num
  have hwin : ∀ u ∈ P.window ω, 0 < u := fun u hu => lt_of_lt_of_le hlo hu.1
  set c := 4 * (ω + P.η₀) * (1 / (1 + (ω - P.η₀) / 120) ^ 20)
  have hrem := P.integrableOn_K_mul (ω := ω) (continuousOn_skewesRem hwin)
  have hremle : ∫ u in P.window ω, P.K (u - ω) * skewesRem u ≤ c :=
    P.integral_K_le_of_le (continuousOn_skewesRem hwin) fun u hu => skewesRem_le hlo hu
  have hint : ∫ u in P.window ω, P.K (u - ω) * (-1 - Z u - skewesRem u) =
      -((∫ u in P.window ω, P.K (u - ω) * Z u) +
        (∫ u in P.window ω, P.K (u - ω) * skewesRem u) + 1) := by
    have e : (fun u => P.K (u - ω) * (-1 - Z u - skewesRem u)) =
        fun u => -(P.K (u - ω) * ((Z u + skewesRem u) + 1)) := by funext u; ring
    rw [e, integral_neg,
      P.integral_K_add (P.integrableOn_K_add hZint hrem) (P.integrableOn_K_const ω 1),
      P.integral_K_add hZint hrem, P.integral_K_const]
  -- `S / Φ(ω) ≥ S_lo / Φ(ω)`, and the numerical margin
  have hΦ : 0 < Phi ω := lt_trans one_pos (one_lt_Phi (by unfold ω; norm_num))
  have hSΦ : (100292 / 10 ^ 5 - 227 / 10 ^ 8) / Phi ω ≤ S / Phi ω :=
    div_le_div_of_nonneg_right hS hΦ.le
  have hmargin : 1 + 555 / 10 ^ 11 + 13250 / 10 ^ 8 + 851 / 10 ^ 22 + c <
      (100292 / 10 ^ 5 - 227 / 10 ^ 8) / Phi ω := by
    simp only [c, KParams.η₀, P, ω, Phi]; norm_num
  obtain ⟨u, hu, hpos⟩ := exists_pos_of_integral_pos (P.measurableSet_window ω)
    (fun u _ => P.K_shift_nonneg ω u) hExpl (by rw [hint]; linarith [neg_div (Phi ω) S])
  rw [window_eq] at hu
  have hu0 : 0 < u := lt_of_lt_of_le (by norm_num) hu.1
  refine ⟨exp u, by rw [log_exp]; exact hu, ?_, ?_, ?_⟩
  · unfold f at hpos
    have hden : 0 < exp (u / 2) * Phi u := mul_pos (exp_pos _) (lt_trans one_pos (one_lt_Phi hu0))
    have h1 : 0 < u * (piR (exp u) - li (exp u)) := by
      have := mul_pos hpos hden
      rwa [div_mul_cancel₀ _ hden.ne'] at this
    linarith [(pos_iff_pos_of_mul_pos h1).mp hu0]
  · calc (2 : ℝ) ≤ 1 + u := by linarith [hu.1]
      _ ≤ exp u := by linarith [add_one_le_exp u]
  · exact exp_lt_skewes_decimal (by unfold L1397 L10; linarith [hu.2])

end SkewesVerify.Skewes
