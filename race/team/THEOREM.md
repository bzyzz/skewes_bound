# Residue against non-residue races modulo 8 and 12

Status: draft. Certificates: `team_cert_q8_T100000_w64.764612.json`, `team_cert_q12_T50000_w339.632308.json` (2026-10-26).

Independent review (2026-10-26): no error affecting L > 0.
- It re-derived the orthogonality identity and the lower-bound directions, and checked the coefficient-1 rescaling of R₃, R₄, R₅ and ε₀.
- It confirmed the characters, conductors and parities, and Platt's heights.
- Its independent float recomputation gives L ≈ 0.1390 (mod 8) and 0.1286 (mod 12).
- It spot-checked the new zeros with its own evaluation of Z.
- One fix was made: the decimal bounds in the statement had been rounded down; they now read 10^{68.8976} and 10^{147.5007}.
- Novelty check done: Bays–Ford–Hudson–Rubinstein (2001) located the mod-8 region heuristically, without a certificate (see below). Dusart: our main terms also follow from Axler's open-access Theorem 2 for log x in [60, 800] (see `../THEOREM.md` (D)).

## Statement

For q ∈ {8, 12}, the squares in (ℤ/qℤ)^× form the subgroup {1}. The "team" race asks when the residue class 1 overtakes the average of the non-residue classes. Equivalently, it asks when

G_q(x) := φ(q)·π(x;q,1) − Σ_{(b,q)=1} π(x;q,b) > 0.

The first sign change is CPNT problem list (PIMS 2024) #9 for these moduli. No explicit value appears to be known. Bailleul–Hayani–Untrau (arXiv:2603.20093) give only asymptotic bounds, conditional on GRH and an effective linear-independence hypothesis.

**Theorem.**
1. π(x;8,1) > ⅓(π(x;8,3) + π(x;8,5) + π(x;8,7)) for some x with log x ∈ [64.764452, 64.764772], so for some x < e^{64.764772} < 10^{28.1270}.
   - This is the region located heuristically by Bays–Ford–Hudson–Rubinstein (2001, Fig. 7).
   - Certificate `team_cert_q8_T100000_w64.764612.json`, with T₀ = 10⁵ and u_lo ≥ 45.
   - The earlier, weaker certificate for log x ∈ [158.64181, 158.64244] (x < 10^{68.8976}) is kept as `team_cert_q8_T50000_w158.642125.json`.
2. π(x;12,1) > ⅓(π(x;12,5) + π(x;12,7) + π(x;12,11)) for some x with log x ∈ [339.631918, 339.632698], so for some x < e^{339.632698} < 10^{147.5007}.

## Derivation [derived]

By orthogonality:

φπ(x;q,1) − Σ_b π(x;q,b) = Σ_{χ≠χ₀} S(χ;x) − φA₁(x) + A(x),

where A₁ is the sum of 1/m over prime powers p^m ≤ x with m ≥ 2 and p^m ≡ 1, and A is the same sum over all p ∤ q.

- **Squares (m = 2).** Every p ∤ q has p² ≡ 1 (mod q). So the m = 2 terms give (1 − φ)·½(π(√x) − ω(q)).
- **Higher powers (m ≥ 3).** Together they are at most φ·1.3x^{1/3}/log x (Ford–Hudson (2.12)).
- **The characters.** The three nonprincipal characters are real and are induced by:
  - q = 8: χ₋₄ (4.3), χ₋₈ (8.3), χ₈ (8.5). There are no imprimitive corrections, since 2 divides each conductor.
  - q = 12: χ₋₃ (3.2), χ₋₄ (4.3), χ₁₂ (12.11). There are corrections at p = 2 for χ₋₃ and at p = 3 for χ₋₄.
- **The explicit formula.** For each primitive character use the corrected explicit formula (`../THEOREM.md` §2): S(χ*;x) = −Σ_ρ li(x^ρ) − (1−δ) log log x + C_χ + ϑ(1/x).

Multiply by u e^{−u/2}, bound ½u e^{−u/2}π(e^{u/2}) ≤ Φ_D(u) (Dusart), and drop the positive term (φ−1)ω(q)/2. Here, explicitly:
- ε₀(u) = u e^{−u/2}·[Σ_χ(e^{−u} + (1−δ_χ) log u + |C_χ|) + Σ_χ n_χ(log u + 1.37)] + φ·1.3·e^{−u/6};
- n_χ is the number of primes dividing q but not the conductor (q = 12: one each for 3.2 and 4.3);
- 1.37 bounds 1 − log log p for p = 2 (it is 0.91 for p = 3);
- the m ≥ 3 bound needs x ≥ e³⁰ (Ford–Hudson (2.12));
- φ·1.3 is used where (φ−1)·1.3 would suffice.

This gives

u e^{−u/2}G(eᵘ) ≥ −(φ−1)Φ_D(u) − Σ_χ Σ_ρ u e^{−u/2} li(e^{ρu}) − ε₀(u).

Integrate against the two-scale kernel exactly as in `../THEOREM.md` §4, with each character's coefficient 1 instead of 2. This gives

I(ω) ≥ L := −(φ−1)Φ_D(u_lo) − Σ_χ S_χ(ω,T₀) − R₃ − R₄ − R₅ − ε₀.

If L > 0, then G_q(eᵘ) > 0 for some u ∈ [ω−η₀, ω+η₀]. The code is in `team_criterion.py` (error terms) and `team_certify.py` (Arb), and it reuses the reviewed race code.

**Validation.**
- At ω = 16.5, 17.3 and 18.1, the smoothed u e^{−u/2}G(eᵘ) computed directly from the primes agrees with the bias plus the zero sum to within 0.01, while the values themselves range from −2.4 to −4.8 (`validate_team.py`).
- The independent review included all lower-order terms exactly: π(√x), the m ≥ 3 prime powers, the constants C_χ, the imprimitive corrections, and the 2/(ρ³u²) term. With them, the agreement is within 2.3·10⁻⁵ at ω ∈ {15, 16, 17, 18, 18.9} for both moduli.
- This confirms the bias (φ−1) = 3, the coefficient 1 per character, and the sign.
- A direct count also gives G(p) < 0 at every prime p ≤ 2.18·10⁸ for both moduli.

## Inputs

- **Zeros.** All zeros of the five L-functions involved, certified in Arb, up to T₀ = 10⁵ for mod 8 and up to T₀ = 5·10⁴ for mod 12. Turing's method proves the lists complete and the zeros simple and on the critical line.
  - The new characters 8.3 and 12.11 have 80,066 and 83,292 zeros, from `team/run_zeros_team.sh`.
  - The characters 3.2, 4.3 and 8.5 reuse `../data`.
- **Cited:** Platt's GRH heights A_χ: 8.3 and 8.5 at 1.25·10⁷, 12.11 at 8.3·10⁶, 3.2 at 3.3·10⁷, 4.3 at 2.5·10⁷. Also BMOR (2021), Dusart (2018), and Ford–Hudson (2.11)–(2.12), as in `../THEOREM.md`.

## Certified numbers

| | mod 8 | mod 12 |
|---|---|---|
| ω; η₁ (η₁T₀); δ; M; core | 158.642125; 2.8·10⁻⁴ (14); 3.5·10⁻⁵; 60; κ₆ | 339.632308; 2.8·10⁻⁴ (14); 1.1·10⁻⁴; 60; κ₆ |
| S_χ | 4.3: −0.648; 8.3: −0.7029; 8.5: −1.8277 | 3.2: −0.3566; 4.3: −0.6484; 12.11: −2.142 |
| Σ S_χ | −3.179 ± 4.8·10⁻⁴ | −3.147 ± 6.4·10⁻⁴ |
| bias 3Φ_D(u_lo) | 3.03882 | 3.01788 |
| R₃; R₄; R₅; ε₀ | 7.6·10⁻⁶; 9.61·10⁻⁴; 1.1·10⁻⁷; 1.7·10⁻¹¹ | 1.3·10⁻⁶; 9.44·10⁻⁴; 5.0·10⁻⁷; 1.4·10⁻²⁴ |
| **L** | **0.139 ± 2.6·10⁻⁴ > 0** | **0.129 ± 5.4·10⁻⁴ > 0** |

## Heuristics and prior art

- **Prior art for mod 8.** P_{8,N,R}(x) in Bays–Ford–Hudson–Rubinstein (J. Number Theory 2001) is exactly −G₈ up to normalization. Their Figure 7 plots "the first negative region for P_{8,N,R}(x; 10000)", near log₁₀ x ≈ 28.127, and calls it "possibly the first region where P_{8,N,R}(x) becomes negative". They give no certificate. Our full-weight sums confirm that region:

  | u | log₁₀ x | full weight, zeros ≤ 5·10⁴ | width |
  |---|---|---|---|
  | 64.76460 | 28.1269 | −0.061 | 1.8·10⁻⁴ |
  | 94.45984 | 41.0234 | −0.018 | 1.0·10⁻⁴ |
  | 158.64255 | 68.8976 | −0.187 | 1.6·10⁻³ |

  An earlier version of this note wrongly claimed that 10^{68.90} is the first region; the smoothed (JP) scan had missed the two narrower regions. The region near 10^{28.127} is now certified (Theorem 1), and the weaker bound x < 10^{68.8976} is kept only as a second certificate.
  
  It was certified after lowering the criterion's requirement to u_lo ≥ 45 and extending the χ₋₈ zeros to 10⁵ (see the section below).
- **Mod 12.** A full-weight rescan (zeros ≤ 5·10⁴; candidates where the JP(10⁴) value is below 0.35) finds no negative region below u = 339.63. So 10^{147.50} is very probably the first sign change. No prior art for mod 12 was found.
- **Mod 24.** For q = 24 (7 non-residue classes) no sign change occurs up to 10^{1300} heuristically. The closest approach is +0.2 at u ≈ 1537.

## Mod 8 near 10^{28.127} (2026-10-26)

- **Change:** the team criterion now requires only u_lo ≥ 45. Dusart needs u ≥ 44.22; Ford–Hudson (2.12) needs x ≥ e³⁰; (2.10) needs x > 10; the R₅ disk argument and the li expansions hold for any u_lo > 0.
- **Zeros:** L(s,χ₋₈) was extended to 10⁵: 171,163 zeros, Turing passing at every multiple of 5000, L > 0 on [½, 2]. The 5·10⁴ run is kept in `../data/arb_8_3_T5e4`.
- **Parameters:** ω = 64.764612, η₁ = 1.4·10⁻⁴ (η₁T₀ = 14), δ = 2·10⁻⁵, M = 40, κ₆ core, T₀ = 10⁵.
- **Zero sums:** S_χ is −0.3266 (4.3), −1.7005 (8.3) and −1.1225 (8.5); Σ = −3.150 ± 4.0·10⁻⁴.
- **Error terms:** bias 3Φ_D(u_lo) = 3.09901; R₃ = 4.6·10⁻⁵; R₄ = 1.02·10⁻³; R₅ = 3.5·10⁻⁸; ε₀ = 1.07·10⁻⁴ (larger at this smaller u).
- **Result:** **L = 0.0495 ± 8.9·10⁻⁵ > 0.**
- **Follow-up review (2026-10-26), no error.**
  - No dependence on u ≥ 100 remains anywhere in the chain. The only u-thresholds are Dusart (u ≥ 44.22), Ford–Hudson (2.12) (x ≥ e³⁰; the constant 1.3 was re-checked directly: it is ≤ 1.112 for u ≥ 45) and (2.10) (x > 10).
  - An independent float recomputation gives L ≈ 0.0494.
  - The explicit formula matches the primes to 1.2·10⁻⁶ at ω ∈ [15, 19].
  - The pointwise formula is positive on u ∈ [64.76453, 64.76470], peaking at +0.065, inside the certified window.
  - The 8.3 zeros were spot-checked with an independent Z evaluation.
