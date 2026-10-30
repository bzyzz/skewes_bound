# Census of first sign changes in prime races: feasibility (2026-10-27)

Everything here is floating point or literature unless marked otherwise. Sources: `LITERATURE.md` (literature sweep), `count_races.c` (direct counts), `team_scan.py` (explicit-formula scans), `ld_team.py` (random-phase model).

## 1. Two-way races π(x;q,a) vs π(x;q,b), q ≤ 30

- **Direct count** (`count_races 1e10 30`): every ordered pair with q ≤ 30 has crossed by 10⁹, except 3;1,2, 8;1,3, 8;1,7, 12;1,{5,7,11} and 24;1,b. All of these involve q | 24, b = 1.
- **The exceptions are covered by the literature** (`LITERATURE.md` §1–2):
  - exact first crossings (Bays–Hudson; A. Martin; Shchebetov, exhaustive to 10¹⁵; Chaffin, OEIS) for 3;2, 8;7, 12;5, 12;7, 24;13, 24;17, 24;19, 24;23;
  - Ford–Hudson 2001 certificates for 8;3 (≈ 6.7·10¹⁸), 12;11 (10^81.44), 24;5 (10^70.0) and 24;7 (10^40.18), plus our 24;11 (10^48.44).
- **Our heuristic scans** (`../race/heuristics/targets_m*.log`) put the first region of 8;3, 12;11, 24;5 and 24;7 exactly in Ford–Hudson's certified windows. Re-certifying them would only narrow windows.
- **Conclusion: for q ≤ 30 the two-way census is essentially complete in the literature.** A census of two-way races would only matter beyond q = 30 (Fiorilli–Martin moduli 40–168; complex characters, `float_zeros_any.py` is ready and validated against Arb's `hardy_z`).

## 2. Residue vs non-residue ("team") races — Fiorilli's x_q

x_q = first x with avg_{R} π(x;q,·) > avg_{NR} π(x;q,·), where R are the squares of units and ρ = [G : G²]. Distinct races are determined by the set of real characters mod q.

| q | ρ | x_q | how |
|---|---|---|---|
| 5 (=25) | 2 | 2,082,927,221 | BFHR (literature); reproduced by our count and by our scan (u = 21.458) |
| 3, 4, 7, 11, 13, 17, 19, 23, 29 | 2 | known / small | count (3: Bays–Hudson 6.1·10¹¹) |
| 8 (=16) | 4 | < 10^28.1270 | **certified** (our paper, Thm 1.3a) |
| 12 | 4 | < 10^147.5007 | **certified** (our paper, Thm 1.3b) |
| 15 | 4 | heuristic first region u ≈ 192.641 (10^83.66), wide (0.012), deep (+0.46) | scan, zeros ≤ 10⁴ |
| 20 | 4 | heuristic u ≈ 41.436 (≈ 10^18.0), narrow (7.8·10⁻⁴), +0.09 | scan |
| 21 | 4 | heuristic u ≈ 29.077 (≈ 4.2·10¹²), very narrow, +0.008 | scan; within reach of a direct count |
| 28 | 4 | heuristic u ≈ 40.217 (≈ 10^17.47), wide (0.013), deep (+0.40) | scan |
| 24 | 8 | heuristic first region **u ≈ 13619.1292 (x ≈ 10^5914.7)**: +0.30 / +0.36 / +0.40 with zeros ≤ 5·10³ / 10⁴ / 2·10⁴, width ≈ 0.002. None with u < 13619 | scan to u = 30000 (`scans/team_q24_U30000.json`); the model predicted u ≈ 1.3·10⁴ beforehand |

No first crossing or certificate is in the literature for ρ ≥ 4 (LITERATURE.md §3), other than our mod 8 and 12 certificates.

**Checks of the tools.**
- `count_races` reproduces Leech's 26,861 and BFHR's 2,082,927,221.
- `team_scan` reproduces 2,082,927,221 (mod 5), the certified mod 8 region (u = 64.7646) and the mod 12 region (u = 339.632).

## 3. Fiorilli's conjecture and the random-phase model

- **Fiorilli (ANT 2014, Conj. 1.19):** log log x_q ≍ ρ(q)/log rad(q). Bailleul–Hayani–Untrau disprove it in general. No numerical test exists.
- **Data so far, all ρ = 4:**

  | q | q′ = rad q | log log x_q |
  |---|---|---|
  | 8 | 2 | 4.17 |
  | 12 | 6 | 5.83 |
  | 15 | 15 | 5.26 |
  | 20 | 10 | 3.72 |
  | 21 | 21 | 3.37 |
  | 28 | 14 | 3.69 |

  These do not follow ρ/log q′ monotonically (12 and 15 are large); more moduli are needed.
- **Random-phase model** (`ld_team.py`): the Rice rate of h = Σ w_j cos(γ_j u + φ_j) ≥ ρ − 1 under LI, with zeros truncated at T_eff.
  - With T_eff ≈ 300 it is unbiased on the seven known races, with a scatter of a factor ≈ 2.3 in u. That is about the spread of an exponential waiting time, so the model is consistent with the data.
  - For q = 24 it predicts u ≈ 1.3·10⁴ (x ≈ 10^5600).

## 4. Scope available

Team races with q ≤ 200 (q ≢ 2 mod 4), distinct by character set:
- ρ = 2: 46 races (Chebyshev's residue races mod p; mostly small);
- ρ = 4: 47 races, needing zeros of 49 more real characters (≈ 1.5 h);
- ρ = 8: 16 races;
- ρ = 16: 2 races (q = 120, 168).

## 5. Certification notes

- **u ≥ 45** (q = 8, 12, 15): the team criterion in `../race/team/` applies as is, given Arb zeros and Turing for the new characters (15.14; 5.4 is in `../mertens/zeros` to 2·10⁴).
- **40 ≤ u < 45** (q = 20, 28): needs an explicit bound for π(y) at y ≈ e^20, below Dusart's range. Small-x bounds exist, or π(√x) can be computed exactly, since √x ≈ 5·10⁸.
- **u ≈ 29** (q = 21): direct count to ≈ 5·10¹² (parallel two-pass sieve, a few hours), or an exact-data criterion.
- **u ~ 10⁴** (q = 24, ρ = 8): the mollifier floor δ ≳ eπu/(4A_χ) ≈ 7·10⁻³ (A = 4.1·10⁶ for conductor 24). This is certifiable only if the positive region found is wider than that. R₅ needs M ≈ u/4, which is feasible.

## 6. Mod 24 (ρ = 8) after the deep scan

- **First region.** The scan to u = 30000 finds the first region at u ≈ 13619.129, matching the random-phase prediction made before the scan. The next regions are at u ≈ 13721.42 (tiny), 16971.90, 17191.44, 22873.52 (+0.45, width 0.0029) and 28288.34.
- **Unconditional certificate: not feasible.** With Platt's GRH heights (A ≈ 4.1·10⁶ for conductor 24), the mollifier needs δ ≳ eπu/(4A) ≈ 7·10⁻³. That is wider than the positive region (≈ 2·10⁻³). It would need GRH verified to ≈ 10⁸ for seven L-functions, about 2·10⁹ zeros.
- **Under GRH: feasible.** R₅ = 0, so δ can be tiny. With η₁ ≈ 5·10⁻⁴ we need T₀ ≈ 2.4·10⁴. Our Arb zeros reach 10⁵ for six of the seven characters; 24.5 still needs them.
