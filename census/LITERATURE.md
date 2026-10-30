# Literature census: first sign changes of prime number races, q ≤ 30

Compiled 2026-10-27 for the planned census paper (certified first crossings or upper bounds
for all biased races with q ≤ ~30). The goal is to know what is already known, so it can be
cited rather than claimed.

## Conventions

* Δ_{q;a,b}(x) = π(x;q,a) − π(x;q,b). In a "biased race", a is a quadratic **non**-residue
  (the usual leader) and b is a residue. "First crossing" means the least x with
  Δ_{q;a,b}(x) < 0, i.e. π(x;q,b) > π(x;q,a). This x is always a prime ≡ b (mod q).
* Most sources state it as "first x with Δ = −1" (Bays–Hudson, Shchebetov/OEIS) or
  "first time b takes the lead" (OEIS A275939). These are the same.
* Team race (R vs N): R = squares in (Z/q)^×, N = non-squares. Fiorilli's normalisation is
  π(x;q,N) vs (ρ(q)−1)·π(x;q,R) with ρ(q) = [G:G²], which equals the number of square roots
  of 1. BFHR write it as P_{q;N,R} = Σ_N π − c(q)·Σ_R π with c(q) = #N/#R. These two
  normalisations agree.
* **Caveat on "first crossing" for weakly biased races.** For many races the residue class
  leads at very small x, before the bias sets in. Examples: mod 7 R vs N at x = 2; 1 vs −1
  mod 5 at x = 11. The census must fix a definition, either "least x" or "least x after a
  stated threshold". Literature values below use "least x" unless noted.

### Verification tags

| Tag | Meaning |
|---|---|
| **[P]** | Read in the primary source (PDF text or page image). |
| **[S]** | Read in a reliable secondary source (named), not in the original. |
| **[O]** | Read in OEIS. This is data, not peer reviewed. |
| **[C]** | Recomputed here by a direct sieve to 2·10⁷ (scratchpad). Agrees. |
| **UNVERIFIED** | Claimed somewhere but not confirmed. Where it is claimed is stated. |

### Sources read in full or in the relevant part

* Ford–Hudson, Acta Arith. 100 (2001) 297–314.
  [PDF](https://ford126.web.illinois.edu/wwwpapers/lehman.pdf). Pages 1–3 and 10–16 read as
  images.
* Bays–Ford–Hudson–Rubinstein, J. Number Theory 87 (2001) 54–76.
  [PDF](https://ford126.web.illinois.edu/wwwpapers/bfhr.pdf).
* Deléglise–Dusart–Roblot, Math. Comp. 73 (2004) 1565–1575, DOI 10.1090/S0025-5718-04-01649-7.
  Only the author preprint was read, via a Wayback snapshot of
  `math.univ-lyon1.fr/~roblot/resources/pikl.pdf` (2017-08-09). The AMS PDF was not
  accessible.
* Fiorilli, Algebra Number Theory 8 (2014) 1733–1767, DOI 10.2140/ant.2014.8.1733.
  Published PDF from msp.org; arXiv:1210.6946v2 also read.
* Fiorilli–Martin, J. reine angew. Math. 676 (2013) 121–212. arXiv:0912.4908v3 read.
* Granville–Martin, "Prime number races", Amer. Math. Monthly 113 (2006). arXiv:math/0408319
  read.
* Bailleul–Hayani–Untrau, arXiv:2603.20093v2 (24 Apr 2026).
* Comparative prime number theory problem list, arXiv:2407.03530 (local copy `cpntpl.txt`).
* Martin et al., annotated bibliography, arXiv:2309.08729 (local copy `bib.txt`).
* Schlage-Puchta, arXiv:1105.1619 (= Acta Math. Hungar. 102 (2004) 305–320).
* Sneed, "Lead changes in the prime number race" (preprint, local copy).
* A. S. Shchebetov, EUCYS 2018 presentation
  ["Testing Chebyshev's bias for prime numbers up to 10^15"](https://math101.guru/wp-content/uploads/2018/09/01-A3-Presentation-v7.3EN-no.pdf),
  30 slides. Also his ISEF abstract
  ["Testing Chebyshev's bias for prime numbers up to 5·10^15"](https://abstracts.societyforscience.org/Home/PrintPdf/16842).
* OEIS entries, read in `fmt=text` on 2026-10-27: A275939, A007350, A007352, A051025,
  A297006, A297448, A295354, A297355, A297357, A295356, A297450, A298821, A306499, A306500,
  A321857, A321859, A321860, A321858, A071838, A320857, A329242, A392284, A326615, A329224,
  A398234–A398236, A398377, A399143, A399174, A401243–A401245.

### Sources not read directly

These are cited through [S]:

* Leech 1957
* Bays–Hudson 1977, 1978a, 1978b, 1978c, 1979, 1983, 2000
* Rubinstein–Sarnak 1994
* Feuerverger–Martin 2000

---

## 1. Known exact first crossings (direct computation)

### 1a. Two-way races with b = 1 and q | 24 (the 15 "most unfair" races)

Each Fiorilli–Martin density in column 4 is δ(q;a,1), the logarithmic density of
{π(x;q,a) > π(x;q,1)}. The rank is its position in their top-120 list (§4).

| q | a (leader) | b | δ(q;a,1) [FM rank] | First x with π(x;q,1) > π(x;q,a) | Search limit (exhaustive unless noted) | Source(s) / status |
|---|---|---|---|---|---|---|
| 3 | 2 | 1 | 0.999063 [10] | **608,981,813,029** | Exhaustive to ~6.1·10¹¹ (BH 1978). Shchebetov: a 2nd zone starts at 6,148,171,711,663 and no other zone ≤ 10¹⁵. | Bays–Hudson, Math. Comp. 32 (1978) 571–576, MR0476616 [S: annotated bib [137]; FH p.2; BFHR §4; G–M p.3]. OEIS A297006/A007352 [O]. 2nd zone: Shchebetov 2018 [P, slide 19]; A297006 [O], "first two zones up to 10¹³". BFHR predicted the 2nd zone at ≈6.150·10¹² and a 3rd near 3.96555843·10¹⁹ (heuristic) [P]. |
| 4 | 3 | 1 | 0.995928 [15] | **26,861** | Leech: to 3·10⁶. Shchebetov: all 10 zones ≤ 10¹⁵ (see list below). | Leech, J. London Math. Soc. 32 (1957) 56–58 [S: bib [48]; BFHR §4 table; FH p.11] [C]. |
| 8 | 3 | 1 | 0.999569 [7] | **none known** | None ≤ 10¹² (Bays–Hudson 1979; FH p.11 (3.1)). None ≤ 10¹⁵ (Shchebetov 2018). | [P: FH p.11; Shchebetov slides 22/27]. FH certify a negative value in x ∈ [10^18.825, 10^18.842] (§2). No exhaustive search ≥ 10¹⁶ found. |
| 8 | 5 | 1 | 0.997395 [14] | **588,067,889** | Bays–Hudson 1979 to 10¹²: every prime < 10⁹, then every 10⁷ with refinement (not strictly exhaustive above 10⁹). Shchebetov: exhaustive to 10¹⁵, 7 zones. | Bays–Hudson, Internat. J. Math. Math. Sci. 2 (1979) 111–119, MR529694 [S: BFHR §4 "computed Δ_{8,5,1}(x) = −1 for x = 588,067,889"; G–M p.29; bib [141]]. OEIS A297448 [O]. 2nd zone starts 35,615,130,497 (BH 1979, per Shchebetov slide 10). |
| 8 | 7 | 1 | 0.998939 [11] | **192,252,423,729,713** | None ≤ 10¹² (BH 1979). Exhaustive to 10¹⁵ (Shchebetov); only 1 zone. | A. Martin, OEIS A275939 a(8), Aug 2016 [O]; Shchebetov 2017, A295354 [O, P slide 22]. FH p.16 / BFHR §4: Δ_{8,7,1}(1.9282·10¹⁴) = −105 by Ford's Meissel program. That is a point inside the zone, not the first point [P]. |
| 12 | 5 | 1 | 0.999206 [8] | **25,726,067,172,577** | OEIS: "full sequence checked up to 10¹⁴". Project range 10¹⁵. | Shchebetov 2017, A297355 [O, P slide 24]. |
| 12 | 7 | 1 | 0.998606 [13] | **27,489,101,529,529** | As above. | Shchebetov 2017, A297357 [O, P slide 24]. |
| 12 | 11 | 1 | 0.999977 [3] | **none known** | ≤ 1.1·10¹⁴ (A. Martin 2016), ≤ 10¹⁵ (Shchebetov), **≤ 10¹⁹ (B. Chaffin, 2026)**. | A275939 comments [O]: "for n=12, residue 1 does not lead −1 up to 10^19" (Chaffin, Jun 11 2026, corr. Jul 08 2026). Not peer reviewed. FH certify a crossing near 10^81.44 (§2). Ford (private communication in A275939) expects the first crossing ≈ e^187.536 ≈ 2.79·10⁸¹. |
| 24 | 5 | 1 | 0.999988 [1] | **none known** | ≤ 10¹⁵ (Shchebetov 2018) | [P slide 25]. FH window ≈ 10^70.0. |
| 24 | 7 | 1 | 0.999834 [5] | **none known** | ≤ 10¹⁵ | [P slide 25]. FH window ≈ 10^40.17. |
| 24 | 11 | 1 | 0.999983 [2] | **none known** | ≤ 10¹⁵ | [P slide 25]. FH window ≈ 10^352.84–10^353.01. FH Table 3 has an uncertified K* candidate at ω = 111.54595 (≈10^48.44). |
| 24 | 13 | 1 | 0.998722 [12] | **978,412,359,121** | Region near 10¹² found by Bays–Hudson 1978, exact start not given. Exhaustive to 10¹⁵ (6 zones), and a 7th zone in (10¹⁵, 10¹⁶). | Bays–Hudson, J. reine angew. Math. 299/300 (1978) 234–237, MR0472726 [S: bib [139]; FH p.14 quotes Δ_{24,13,1}(9.866·10¹¹) = −6091]. Exact first point: Shchebetov, A295356 [O, P slide 25]. A Jianing Song comment (Oct 2026) in A295356 says the 1st zone ends at 1,096,157,101,043. |
| 24 | 17 | 1 | 0.999125 [9] | **617,139,273,158,713** | Exhaustive to 10¹⁵ (3 zones) | Shchebetov 2018, A297450 [O, P slide 25]. |
| 24 | 19 | 1 | 0.999719 [6] | **706,866,045,116,113** | Exhaustive to 10¹⁵ (5 zones) | Shchebetov 2018, A298821 [O, P slide 25]. |
| 24 | 23 | 1 | 0.999889 [4] | **7,390,203,183,507,483,529** (≈ 7.39·10¹⁸ = e^43.447) | Implicitly exhaustive to that x (A275939 asks for the first lead). Shchebetov: none ≤ 10¹⁵. | **OEIS-only:** A275939 comment by Benjamin Chaffin (Jun 11 2026, corrected Jul 08 2026) [O]. No paper found. The number is prime and ≡ 1 (mod 24) (checked). It agrees with FH Table 3, which has an uncertified candidate at ω = 43.45318 (K* = −0.0204), and with Ford's prediction e^43.453 = 7.437·10¹⁸. Treat as **UNVERIFIED-by-publication**; the method and exhaustiveness are not documented. |

#### Search-limit notes

* **Shchebetov project limit.** The EUCYS presentation says "up to 10¹⁵ … extended to 10¹⁶ in
  May 2018" (slide 18). The later ISEF abstract says "25 new sign-changing zones for 9 prime
  number races … between 10¹² and 5·10¹⁵". So the exhaustive limit for the "not found" races
  is at least **10¹⁵ (documented)**. It is possibly 5·10¹⁵ or 10¹⁶, but that is UNVERIFIED
  per race.
* **Shchebetov method.** primesieve plus prime-by-prime counting, with at least two passes per
  range. A Mathematics of Computation article was announced as "under work" (slide 30). No
  publication was found.
* **Bays–Hudson ≤ 10¹².** The q = 8 checks above 10⁹ were sampled every 10⁷ with refinement
  near 0 (bib [141]). FH p.11 nevertheless state "Δ_{q,b,1}(x) ≥ 0 (x ≤ 10¹²) for q = 8,
  b ∈ {3,7} and q = 24, b ∈ {5,7,11,17,19,23}" [P].

#### All known zones for q = 4

Region starts for Δ_{4,3,1} < 0, i.e. x where π(x;4,1) > π(x;4,3):

1. 26,861 (Leech 1957)
2. 616,841 (Leech)
3. 12,306,137 (Lehmer 1969)
4. 951,784,481 (Lehmer 1969)
5. 6,309,280,709 (Bays–Hudson 1979)
6. 18,465,126,293 (Bays–Hudson 1979)
7. 1,488,478,427,089 (Bays–Hudson 1996)
8. 9,103,362,505,801
9. 64,083,080,712,569
10. 715,725,135,905,981

Zones 1–7 are from the BFHR §4 table [P]. Zones 8–10 are from Shchebetov (slide 20 [P];
A051025 [O]), exhaustive to 10¹⁵. Further notes:

* BFHR reported "the 8th region … near x = 9.318·10¹²", verified with Ford's program. That
  is a point inside zone 8.
* Granville–Martin (p.2) give 6,309,280,697 and 18,465,126,217 for zones 5 and 6. These
  disagree with BFHR/Shchebetov (…709, …293); the cause is unresolved and does not affect
  first crossings.
* Deléglise–Dusart–Roblot found a further region starting ≈ 9.97·10¹⁷ and extending to
  ≈ 1.005·10¹⁸, so it surrounds 10¹⁸ but δ(10¹⁸) > 0 [P, preprint Table 1]. Their search is
  explicitly **not exhaustive**: they sampled δ(1000·1.004ⁿ) up to x = 1,088,537,721,123,564,252
  and refined near small values. They computed π(x;4,1) and π(x;4,3) at x = d·10^j
  (1 ≤ d ≤ 9, 10 ≤ j ≤ 19) and at 10²⁰ (Table 2).
* Because zones 9 and 10 lie below the DDR region, the DDR region is at least the 11th zone.
  DDR call it "region 9" because their list (Table 1) had 8 earlier ones.
* DDR give **only mod 4**. There is no table of first sign changes for other moduli.

### 1b. Two-way races "1 vs −1" for all q ≤ 30 (OEIS A275939)

Definition: "Consider the prime race mod n between n·k+1 and n·k−1. Terms are numbers n·k+1
where n·k+1 first takes lead over n·k−1." The sequence is by Andy Martin (Aug 2016) [O].
Values for q ≤ 30 except 3, 4, 6, 8, 12, 24 were recomputed here [C].

| q | first x | q | first x | q | first x |
|---|---|---|---|---|---|
| 3 | 608,981,813,029 | 13 | 53 | 22 | 23 |
| 4 | 26,861 | 14 | 71 | 23 | 47 |
| 5 | 11 | 15 | 331 | 24 | 7,390,203,183,507,483,529 (Chaffin 2026) |
| 6 | 608,981,813,017 (tie point of mod 3, since 2 ∉ class 5 mod 6) | 16 | 17 | 25 | 101 |
| 7 | 71 | 17 | 239 | 26 | 53 |
| 8 | 192,252,423,729,713 | 18 | 37 | 27 | 379 |
| 9 | 37 | 19 | 213,827 | 28 | 29 |
| 10 | 11 | 20 | 1,381 | 29 | 59 |
| 11 | 23 | 21 | 673 | 30 | 331 |
| 12 | unknown; none ≤ 10¹⁹ | | | | |

* A275939 comments: "The longest n·k+1 versus n·k−1 races up to n = 999 are for
  n = 3, 6, 8, 12, 24 and 168. When n = 168 the race ends at prime 273084304417." The mod 12
  and mod 24 races were checked to 1.1·10¹⁴ (Martin 2016).
* **No source tabulates the other two-way pairs (a, b) with b ≠ 1** (b a residue,
  a a non-residue) for q ∈ {5, 7, 9, 11, 13, 15, 16, 20, 21, …}. This is GAP 1. These
  crossings are mostly tiny and need no literature citation; compute them.

### 1c. Team races (R vs N, and coset races)

| Race | First x with residue side ahead | Search limit | Source / status |
|---|---|---|---|
| q=3, q=4 | as in 1a | | |
| q=5: N={2,3} vs R={1,4} | **2,082,927,221** | (BFHR plots to 10¹²) | BFHR p.10: "the first sign change in P_{5;N,R}(x) occurs at x = 2082927221" [P]. Also OEIS A321857 (J. Song 2019), A306499 (D=5), A392295 [O]. |
| q=7: N={3,5,6} vs R={1,2,4} | R leads at x=2 (2 is a QR mod 7). The next first lead is at **48,673**. | | A321859 [O] [C] |
| q=8, χ₈: {3,5} vs {1,7} | **11,100,143** | | A071838 (J. Song 2019); A306499 (D=8) [O] [C] |
| q=8, χ₋₈: {5,7} vs {1,3} | trivial (x=3) | | A320857, A306500 (D=−8) [O] |
| q=12, χ₁₂: {5,7} vs {1,11} | **61,463** | | A321858; A306499 (D=12) [O] [C] |
| q=13 (D=13) | **2,083** | | A306499 [O] [C] |
| q=21 (D=21) | **1,217** | | A306499 [O] [C] |
| q=29 (D=29) | **719** | | A306499 [O] [C] |
| D = −7, −11, −15, −19, −20, −23, −24, 17, 24, 28 | 2, 5, 2, 11, 3, 2, 5, 2, 5, 3 (trivial) | | A306500 / A306499 [O] |
| **q=8 full N vs R**: Σ_{b=3,5,7} π(x;8,b) − 3π(x;8,1) (Fiorilli's x_8, BFHR P_{8;N,R}) | **none known** | Bias over [1,10¹²] "is exactly 1.000000" (BFHR) | BFHR Fig. 7 plots the truncated explicit formula P_{8;N,R}(x;10⁴) dipping below 0 near log₁₀ x ≈ 28.127 and call it "possibly the first region" [P]. This is **heuristic, not certified**. A329242 [O] repeats it ("negative for some n ~ 10^28.127"). Granville–Martin p.29 say π(x;8,1) first leads the 4-way race "beyond 10^28", also heuristic. |
| **q=12, 24, 15, 16, 20, 21, 28 full N vs R** (ρ(q) ≥ 4) | **none known** | | No literature found (GAP 3). |
| q=24, index-4 coset races (2 classes vs 2 classes), J. Song, OEIS, Jul–Oct 2026 [O] | See list below | | Not peer reviewed. |

Index-4 coset races mod 24 (OEIS, Jianing Song). For each: first x with the race below 0,
and the limit to which no later x is known.

| H vs S | First x | No further terms up to | OEIS |
|---|---|---|---|
| H={1,19} vs S={5,23} | 5,279,549,641 | 1.131·10¹² | A398234 |
| H={1,19} vs S={7,13} | 491,923 | 9.29·10¹¹ | A398235 |
| H={1,19} vs S={11,17} | 52,551,817 | 2.30·10¹² | A398236 |
| H={1,5} vs S={7,11} | 5 (then 67,137,387,817) | 1.28·10¹² | A398377 |
| H={1,5} vs S={13,17} | 5 (then 6073) | 1.26·10¹² | A399143 |
| H={1,5} vs S={19,23} | 5 | 1.03·10¹² | A399174 |
| H={1,23} vs S={5,19} | 46,698,887 | 1.44·10¹² | A401243 |
| H={1,23} vs S={7,17} | 77,017 | 1.37·10¹² | A401244 |
| H={1,23} vs S={11,13} | 81,119 | 1.36·10¹² | A401245 |

The mod 8 and mod 12 two-way races are the remaining index-4 cosets: H={1,17} is 1 mod 8,
and H={1,13} is 1 mod 12.

---

## 2. Known certified upper bounds (Lehman-type explicit-formula certificates)

### 2a. Ford–Hudson, Acta Arith. 100 (2001) 297–314 (MR1862054)

This is the only Lehman-type certificate for Dirichlet races found in the literature. All
details were read from page images [P].

**Theorem 2** (p.13). Suppose:

* ω−η ≥ 32.3 and 0 < η ≤ ω/30;
* q | 24, (b,q) = 1, 1 < b < q;
* for each χ with χ(b) = −1, all zeros with |Im| ≤ A_χ lie on Re s = ½;
* 150 ≤ T_χ ≤ A_χ and 2A_χ/η ≤ α ≤ A_χ².

Then the Gaussian-kernel average of u·e^{−u/2}(π_{q,b}(e^u) − π_{q,1}(e^u)) over
[ω−η, ω+η] is at most an explicit sum over |γ| ≤ T_χ plus error terms.

**Theorem 3.** For each row with a K value,
min_{ω−η ≤ u ≤ ω+η} u·e^{−u/2}(π_{q,b}(e^u) − π_{q,1}(e^u)) ≤ K < 0.

Setup:

* Zeros are Rumely's, with T_χ = 10⁴ for every χ.
* A_χ = 10⁵ in Table 2 (conductor < 10) and A_χ = 10⁴ in Table 3.
* The result is unconditional, because Rumely verified GRH to height A for these
  L-functions.
* Double precision, results "rounded up in the 4th decimal place".
* Only x ≥ e^32.3 ≈ 10¹⁴ was examined. Only races **b vs 1** with **q | 24** are covered.

Certified windows:

| q | b | ω | η | α | K | Certified x-window [e^{ω−η}, e^{ω+η}] |
|---|---|---|---|---|---|---|
| 3 | 2 | 45.12686 | 0.02 | 10⁷ | −0.0650 | [10^19.590, 10^19.607] |
| 3 | 2 | 58.36855 | 0.02 | 10⁷ | −0.1525 | [10^25.341, 10^25.358] |
| 4 | 3 | 2179.77584 | 0.05 | 4·10⁶ | −0.7761 | [10^946.64, 10^946.69] (this is a π_{4,1} − π_{4,3} > √x/log x type region) |
| 4 | 3 | 78683.67818 | 2.00 | 1.2·10⁵ | −0.8372 | ≈10^34171–10^34173 |
| 8 | 3 | 43.36630 | 0.02 | 10⁷ | −0.0013 | **[10^18.825, 10^18.842] = [6.68, 6.96]·10¹⁸** |
| 8 | 3 | 54.94255 | 0.02 | 10⁷ | −0.0280 | [10^23.853, 10^23.870] |
| 8 | 5 | 32.89388 | 0.02 | 10⁷ | −0.0503 | [1.89, 1.97]·10¹⁴ |
| 8 | 5 | 57.48058 | 0.02 | 10⁷ | −0.1915 | ≈10^24.96 |
| 8 | 7 | 45.34991 | 0.02 | 10⁷ | −0.0508 | [4.86, 5.06]·10¹⁹ |
| 8 | 7 | 48.79950 | 0.02 | 10⁷ | −0.1724 | ≈10^21.19 |
| 12 | 11 | 187.53674 | 0.02 | 10⁷ | −0.0191 | **[10^81.437, 10^81.455]** |
| 12 | 11 | 191.89007 | 0.02 | 10⁷ | −0.0182 | [10^83.328, 10^83.346] |
| 12 | 5 | 103.73736 | 0.03 | 7.5·10⁵ | −0.0445 | **[10^45.040, 10^45.066]** |
| 12 | 7 | 39.12144 | 0.02 | 1.55·10⁶ | −0.1410 | [9.58, 9.98]·10¹⁶ |
| 12 | 7 | 45.87795 | 0.02 | 1.4·10⁶ | −0.0871 | ≈10^19.92 |
| 24 | 5 | 161.18837 | 0.04 | 5.25·10⁵ | −0.0920 | **[10^69.986, 10^70.021]** |
| 24 | 7 | 92.49622 | 0.03 | 8.3·10⁵ | −0.0530 | **[10^40.158, 10^40.184]** |
| 24 | 11 | 812.63677 | 0.20 | 1.18·10⁵ | −0.0104 | **[10^352.837, 10^353.011]** |
| 24 | 13 | 34.14425 | 0.02 | 1.7·10⁶ | −0.3521 | ≈[6.61, 6.88]·10¹⁴ |
| 24 | 17 | 34.19749 | 0.02 | 1.65·10⁶ | −0.0110 | ≈[6.97, 7.25]·10¹⁴ |
| 24 | 19 | 34.20322 | 0.02 | 1.65·10⁶ | −0.1362 | ≈[7.01, 7.29]·10¹⁴ |
| 24 | 23 | 94.46170 | 0.03 | 8·10⁵ | −0.0113 | **[10^41.011, 10^41.037]** |

Notes on the table:

* Window endpoints were converted here as e^{ω±η}. FH give only ω and η.
* Uncertified K* < 0 candidates (no K) also appear in Tables 2–3:
  * (8,5) at 34.46826 and (8,7) at 32.89284;
  * (12,5) at 39.12815, 69.00554, 73.93306, 88.98310 and 102.08460;
  * (24,11) at 111.54595;
  * (24,17) at 34.05708;
  * (24,23) at 43.45318.
* FH p.16: "For each pair (q,b) except (4,3), the first few likely regions of negative values
  are listed." So for (12,11), (24,5), (24,7) and (24,11), FH's K* heuristic found **no
  candidate region** below the first listed ω.

**Corollary 4** (p.15) [P]:

* π_{8,b}(x) < π_{8,1}(x) for some x < 5·10¹⁹ (b = 3, 5, 7).
* π_{12,b}(x) < π_{12,1}(x) for some x < 10⁸⁴ (b = 5, 7, 11).
* π_{24,b}(x) < π_{24,1}(x) for some x < 10³⁵³ (b = 5, …, 23).
* If the zeros of L(s,χ₋₄) up to height 630000 are on the critical line, then
  π_{4,1}(x) − π_{4,3}(x) > √x/log x for some x near e^78683.7.

Corrections and cautions:

* **Rounding in Cor. 4.** The rigorous window ends are e^{45.36991} ≈ 5.06·10¹⁹ for (8,7)
  and e^{812.83677} ≈ 1.025·10³⁵³ for (24,11). These are slightly above the stated "5·10¹⁹"
  and "10³⁵³". When citing, prefer the window ends.
* **Misprint in the annotated bibliography.** Entry [230] of arXiv:2309.08729 says 10²⁵³;
  the paper says 10³⁵³ (p.15) [P].

**Best rigorous upper bound per still-open race (from FH):**

| Race | Upper bound |
|---|---|
| (8,3) | 6.96·10¹⁸ |
| (12,11) | 10^81.455 |
| (24,5) | 10^70.021 |
| (24,7) | 10^40.184 |
| (24,11) | 10^353.011 |

Other notes:

* For (24,23), FH's rigorous bound 10^41.037 is superseded by the OEIS exact value if one
  accepts it.
* For (12,5), (12,7), (24,13), (24,17), (24,19) and (8,7), exact values are now known (§1a).
* The user's own certificate (24,11) < 10^48.4439 (github.com/bzyzz/skewes_bound) shows up in
  web searches. It is not prior literature.

### 2b. Bays–Ford–Hudson–Rubinstein, J. Number Theory 87 (2001) 54–76 (MR1816036)

**No certified bounds** [P]. The paper:

* computes heuristic biases from truncated explicit formulas (Table I);
* plots P_{q;N,R}(x) for x ≤ 10¹² and P_{q;N,R}(x;T) (truncated, T = 10⁴) beyond;
* lists the known region starts for q = 3 and 4 (§4);
* gives P_{5;N,R} first sign change = 2,082,927,221;
* gives the heuristic predictions quoted in §1;
* shows the mod 8 N vs R dip at ≈10^28.127 (Fig. 7), "possibly the first region". This is
  not certified.

For rigorous sign changes it defers to [11] = Ford–Hudson. The "certified" statements in the
paper are exact Meissel-type evaluations, for example Δ_{8,7,1}(1.9282·10¹⁴) = −105 and the
mod 4 8th region near 9.318·10¹².

### 2c. Other upper-bound results (none numerically competitive)

* **Schlage-Puchta**, Acta Math. Hungar. 102 (2004) 305–320 (arXiv:1105.1619) [P]. Assume
  GRH mod q, let f(q) = #{x : x² ≡ 1 (mod q)} and q⁺ = max(q, e^1260). Then there is an
  x < e₂((q⁺)^170 + e^{18 f(q)}) (iterated exponential) such that π(x;q,1) > π(x;q,a) for all
  a ≢ 1.
* **Bailleul–Hayani–Untrau 2026**, Theorems 1.5–1.6, assuming GRH plus an effective LI
  (ELI_A):
  * log log log x_{q;a,b} ≪_A log φ(q);
  * for R vs N, log log log x_{q;R,NR} ≪_A ϱ(q)/log rad(q) + log ϱ(q).

  The constants are not explicit. See §5.
* **Knapowski–Turán** (Comparative prime-number theory I–VIII, 1962–63, incl. VIII for
  k = 8): tower-type bounds for the first sign change (FH p.11: "towering exponentials"). No
  numbers usable for the census.
* **Not found:** any Lehman-type certificate for races mod 5, 7, 9, 11, 13, 15, 16, 20, 21,
  … or for any race with b ≠ 1, or for full team races (N vs R) with ρ(q) ≥ 4.
  * The annotated bibliography lists none.
  * Its Lehman citations are to π−li work: Lehman, te Riele, Bays–Hudson 2000,
    Chao–Plymen, Saouter–Demichel, Saouter–Trudgian–Demichel, Revers 2025.

---

## 3. Fiorilli's conjecture on Skewes numbers of R vs N races

Source: D. Fiorilli, *Highly biased prime number races*, Algebra Number Theory 8 (2014),
no. 7, 1733–1767, DOI 10.2140/ant.2014.8.1733. Read from the **published** PDF, pp.1736 and
1742 (page images) [P].

**Setup** (p.1736):

* A = NR := {a ∈ (Z/qZ)^× : a ≢ □ mod q}, the non-squares. B = R := the squares.
* |B| = φ(q)/ρ(q) and |A| = φ(q)(1 − 1/ρ(q)), where ρ(q) := [G : G²] for G = (Z/qZ)^×:

  | Condition on q | ρ(q) |
  |---|---|
  | 2 ∤ q | 2^{ω(q)} |
  | 2 ‖ q | 2^{ω(q)−1} |
  | 4 ‖ q | 2^{ω(q)} |
  | 8 \| q | 2^{ω(q)+1} |

* The race is π(x;q,NR) > (ρ(q)−1)·π(x;q,R) (§2, p.1742). Lemma 2.1 under GRH:
  E_q(x) = [π(x;q,NR) − (ρ(q)−1)π(x;q,R)] / (√x/log x)
  = ρ(q) − 1 + Σ_{χ²=χ₀, χ≠χ₀} Σ_γ x^{iγ}/ρ_χ + o(1).

**Definitions** (p.1742, verbatim content):

* Generalized Skewes number: x_{q;a,b} := inf{x : π(x;q,a) < π(x;q,b)}. "For large q we
  expect this generalized Skewes number to be small and uninteresting."
* "**x_q := inf{x : π(x;q,NR) < (ρ(q) − 1)π(x;q,R)}** to tend to infinity as ρ(q)/log q′
  tends to infinity (q′ is the radical of q)."
* Theorem 1.6 defines q′ := ∏_{p|q} p.

**Conjecture 1.19.** "As ρ(q)/log q′ tends to infinity we have
**log log x_q ≍ ρ(q)/log q′**."

* There are **no explicit constants**: ≍ means bounded above and below by absolute positive
  constants. In the abstract form, x_q = exp(exp(Θ(ρ(q)/log rad q))).
* The problem list's paraphrase, "roughly exp(exp(ρ(q)/log rad(q)))" (arXiv:2407.03530,
  Problem 9, posed by A. Bailleul), is therefore correct up to the unspecified constants.
* Motivation: "similar arguments to those of [Montgomery 1980; Ng 2000]", i.e. large
  deviations. Theorem 1.6 (GRH+LI):
  exp(−a₁ρ(q)/log q′) ≤ 1 − δ(q;NR,R) ≤ exp(−a₂ρ(q)/log q′).

**arXiv v2 differs** (arXiv:1210.6946v2, 28 Oct 2012) [P]:

* It has typos: the set labels are swapped ("NR := {a ≡ □}").
* x_q is printed as inf{x : (ρ(q)−1)π(x;q,NR) < π(x;q,R)}.
* The published version is the correct one to cite.

**Restatement by Bailleul–Hayani–Untrau** (arXiv:2603.20093, Conj. 1.2) [P]:

* They use averages: (1/|R|)Σ_R π(x;q,a) > (1/|NR|)Σ_NR π(x;q,b) defines x_{q;R,NR}.
* They take ϱ(q) = number of square roots of 1 mod q, which equals [G:G²].
* Hypothesis: ϱ(q_n)/log rad(q_n) → ∞.
* Conclusion: log log x_{q_n;R,NR} ≍ ϱ(q_n)/log rad(q_n).

**Status:**

* **Disproved unconditionally in general** by BHU Theorem 1.3. For any increasing h → ∞ and
  any L > 0 there are square-free sequences q_n, q′_n with ϱ/log q → ∞ such that
  h(x_{q′_n;R,NR}) > L·ϱ(q′_n)/log q′_n and h(x_{q_n;NR,R}) > L·ϱ(q_n)/log q_n.
* The construction uses Chowla's least-quadratic-nonresidue method: all of 2, p₁, …, p_n are
  residues mod q_n.
* BHU "expect the large deviation heuristics used by Fiorilli … should still hold for some
  sequences … for example the primorials".
* **No numerical test** of Conjecture 1.19 was found. Fiorilli's Table 1 (published p.1738)
  gives only densities δ(q;NR,R) for half-primorials:

  | q | ω(q) | ρ(q)/log q′ | δ(q;NR,R) |
  |---|---|---|---|
  | 3 | 1 | 1.82 | 0.999063 |
  | 15 | 2 | 1.47 | 0.999907 |
  | 105 | 3 | 1.71 | 0.999928 |
  | 1155 | 4 | 2.26 | 0.999877 |
  | 15015 | 5 | 3.33 | 0.999950 |
  | 255255 | 6 | 5.14 | 0.9999946 |
  | 4849845 | 7 | 8.31 | 0.999999928 |
  | 111546435 | 8 | 13.81 | 0.999999999954 |

* **Relevance to q ≤ 30.** The conjecture is asymptotic and says nothing quantitative for
  q ≤ 30. For orientation, ρ(q)/log rad(q) (computed here):

  | q | ρ(q)/log rad(q) | Notes |
  |---|---|---|
  | 8 | 5.77 | δ(16;NR,R) = δ(8;NR,R) by Fiorilli Remark 1.4, though q = 16 also gives 5.77 |
  | 24 | 4.47 | |
  | 4 | 2.89 | |
  | 12 | 2.23 | |
  | 3, 9, 27 | 1.82 | |
  | 20 | 1.74 | |
  | 28 | 1.52 | |
  | 15 | 1.48 | |
  | 21 | 1.31 | |
  | 5, 25 | 1.24 | |

* For prime q (ρ = 2), x_q is just the R vs N team race. For q ∈ {8, 12, 15, 16, 20, 21, 28}
  (ρ = 4) and q = 24 (ρ = 8), x_q is the "full team race" of §1c, for which **no first
  crossing and no certificate is known** (except BFHR's heuristic 10^28.127 for q = 8).
* Notably, δ(15;NR,R) = 0.999907 exceeds δ(3;2,1).

---

## 4. Fiorilli–Martin list of the most biased two-way races

D. Fiorilli and G. Martin, *Inequities in the Shanks–Rényi prime number race: an asymptotic
formula for the densities*, J. reine angew. Math. 676 (2013) 121–212, MR3028758.
arXiv:0912.4908v3 (21 Jul 2011) was read [P]. Table numbering in the journal version was
not checked.

**Content.** Theorem 1.11 (GRH+LI): Table 1 gives the ten largest δ(q;a,b). Table 8 gives
"the top 120 most unfair prime number races"; "there are precisely 117 distinct densities
above 9/10" (abstract and §5.4).

**Reading the table:**

* Each entry (q, a, a⁻¹, δ(q;a,1)) represents the whole class of equivalent races, since
  δ(q;a,1) = δ(q;a⁻¹,1) and δ(q;a,1) = δ(q;ab,b) for any non-square a and square b.
* Moduli q ≡ 2 (mod 4) do not appear. Races mod 2m (m odd) coincide with races mod m.

**Moduli in the top 120** (number of entries):

| q | entries | q | entries | q | entries |
|---|---|---|---|---|---|
| 3 | 1 | 20 | 4 | 56 | 3 |
| 4 | 1 | 21 | 2 | 60 | 10 |
| 5 | 1 | 24 | 7 | 72 | 2 |
| 8 | 3 | 28 | 1 (rank 119, below 0.9) | 84 | 6 |
| 12 | 3 | 40 | 10 | 120 | 22 |
| 15 | 4 | 48 | 10 | 168 | 26 |
| 16 | 4 | | | | |

No race mod 7, 9, 11, 13, 17, 19, 23, 25, 27 or 29 reaches 0.9. For comparison:

* BFHR/R–S give team densities δ(7;N,R) = 0.978258, δ(11;N,R) = 0.916795,
  δ(13;N,R) = 0.944319.
* Two-way races are less biased than the corresponding team race. For example, FM's
  δ(5;2,1) = 0.952175, against the team race δ(5;N,R) = 0.995422.

**All 120 entries, ranked** (q, a, a⁻¹, δ(q;a,1)). Ranks 1–117 are above 0.9.

| # | q | a | a⁻¹ | δ | # | q | a | a⁻¹ | δ | # | q | a | a⁻¹ | δ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 24 | 5 | 5 | 0.999988 | 41 | 120 | 91 | 91 | 0.975051 | 81 | 20 | 11 | 11 | 0.931367 |
| 2 | 24 | 11 | 11 | 0.999983 | 42 | 120 | 83 | 107 | 0.975001 | 82 | 168 | 139 | 139 | 0.931362 |
| 3 | 12 | 11 | 11 | 0.999977 | 43 | 120 | 29 | 29 | 0.974634 | 83 | 168 | 55 | 55 | 0.931346 |
| 4 | 24 | 23 | 23 | 0.999889 | 44 | 120 | 19 | 19 | 0.974408 | 84 | 48 | 47 | 47 | 0.929478 |
| 5 | 24 | 7 | 7 | 0.999834 | 45 | 120 | 11 | 11 | 0.971988 | 85 | 168 | 67 | 163 | 0.928944 |
| 6 | 24 | 19 | 19 | 0.999719 | 46 | 48 | 31 | 31 | 0.970470 | 86 | 84 | 71 | 71 | 0.928657 |
| 7 | 8 | 3 | 3 | 0.999569 | 47 | 40 | 7 | 23 | 0.969427 | 87 | 168 | 41 | 41 | 0.927933 |
| 8 | 12 | 5 | 5 | 0.999206 | 48 | 40 | 13 | 37 | 0.969114 | 88 | 84 | 55 | 55 | 0.927755 |
| 9 | 24 | 17 | 17 | 0.999125 | 49 | 120 | 73 | 97 | 0.967355 | 89 | 168 | 71 | 71 | 0.927349 |
| 10 | 3 | 2 | 2 | 0.999063 | 50 | 20 | 19 | 19 | 0.966662 | 90 | 16 | 15 | 15 | 0.926101 |
| 11 | 8 | 7 | 7 | 0.998939 | 51 | 15 | 7 | 13 | 0.964719 | 91 | 168 | 65 | 137 | 0.923960 |
| 12 | 24 | 13 | 13 | 0.998722 | 52 | 120 | 31 | 31 | 0.963190 | 92 | 168 | 53 | 149 | 0.923937 |
| 13 | 12 | 7 | 7 | 0.998606 | 53 | 60 | 13 | 37 | 0.963058 | 93 | 168 | 83 | 83 | 0.923868 |
| 14 | 8 | 5 | 5 | 0.997395 | 54 | 60 | 59 | 59 | 0.962016 | 94 | 21 | 5 | 17 | 0.923779 |
| 15 | 4 | 3 | 3 | 0.995928 | 55 | 40 | 31 | 31 | 0.960718 | 95 | 168 | 79 | 151 | 0.922597 |
| 16 | 120 | 71 | 71 | 0.988747 | 56 | 48 | 5 | 29 | 0.960195 | 96 | 40 | 21 | 21 | 0.922567 |
| 17 | 120 | 59 | 59 | 0.988477 | 57 | 40 | 3 | 27 | 0.960099 | 97 | 168 | 37 | 109 | 0.922359 |
| 18 | 60 | 11 | 11 | 0.987917 | 58 | 16 | 7 | 7 | 0.959790 | 98 | 168 | 17 | 89 | 0.920542 |
| 19 | 60 | 29 | 29 | 0.986855 | 59 | 48 | 11 | 35 | 0.959245 | 99 | 48 | 17 | 17 | 0.918910 |
| 20 | 120 | 109 | 109 | 0.986835 | 60 | 120 | 119 | 119 | 0.957182 | 100 | 56 | 27 | 27 | 0.918015 |
| 21 | 60 | 19 | 19 | 0.986459 | 61 | 15 | 11 | 11 | 0.955226 | 101 | 168 | 59 | 131 | 0.917874 |
| 22 | 120 | 89 | 89 | 0.986364 | 62 | 120 | 41 | 41 | 0.955189 | 102 | 168 | 23 | 95 | 0.917718 |
| 23 | 120 | 79 | 79 | 0.986309 | 63 | 48 | 19 | 43 | 0.952194 | 103 | 168 | 31 | 103 | 0.917278 |
| 24 | 120 | 101 | 101 | 0.984792 | 64 | 5 | 2 | 3 | 0.952175 | 104 | 168 | 29 | 29 | 0.915514 |
| 25 | 15 | 2 | 8 | 0.983853 | 65 | 20 | 13 | 17 | 0.948637 | 105 | 72 | 53 | 53 | 0.913533 |
| 26 | 120 | 13 | 37 | 0.980673 | 66 | 120 | 61 | 61 | 0.948586 | 106 | 21 | 2 | 11 | 0.911872 |
| 27 | 40 | 19 | 19 | 0.980455 | 67 | 60 | 41 | 41 | 0.947870 | 107 | 168 | 19 | 115 | 0.911412 |
| 28 | 60 | 7 | 43 | 0.979323 | 68 | 16 | 3 | 11 | 0.947721 | 108 | 168 | 11 | 107 | 0.909850 |
| 29 | 120 | 23 | 47 | 0.979142 | 69 | 48 | 13 | 37 | 0.946479 | 109 | 168 | 73 | 145 | 0.908239 |
| 30 | 15 | 14 | 14 | 0.979043 | 70 | 40 | 17 | 33 | 0.946002 | 110 | 168 | 5 | 101 | 0.908206 |
| 31 | 120 | 17 | 113 | 0.978762 | 71 | 40 | 11 | 11 | 0.945757 | 111 | 56 | 31 | 47 | 0.906135 |
| 32 | 120 | 7 | 103 | 0.978247 | 72 | 40 | 39 | 39 | 0.942554 | 112 | 84 | 67 | 79 | 0.905578 |
| 33 | 48 | 23 | 23 | 0.978096 | 73 | 60 | 31 | 31 | 0.941802 | 113 | 168 | 13 | 13 | 0.904525 |
| 34 | 120 | 43 | 67 | 0.978013 | 74 | 48 | 7 | 7 | 0.939000 | 114 | 168 | 97 | 97 | 0.904162 |
| 35 | 60 | 17 | 53 | 0.977433 | 75 | 16 | 5 | 13 | 0.938369 | 115 | 72 | 35 | 35 | 0.903755 |
| 36 | 48 | 41 | 41 | 0.977183 | 76 | 168 | 125 | 125 | 0.936773 | 116 | 84 | 47 | 59 | 0.902413 |
| 37 | 40 | 29 | 29 | 0.977161 | 77 | 168 | 155 | 155 | 0.935843 | 117 | 56 | 37 | 53 | 0.900863 |
| 38 | 20 | 3 | 7 | 0.976713 | 78 | 168 | 47 | 143 | 0.932099 | 118 | 84 | 53 | 65 | 0.899063 |
| 39 | 120 | 53 | 77 | 0.976527 | 79 | 168 | 61 | 157 | 0.931981 | 119 | 28 | 11 | 23 | 0.898807 |
| 40 | 60 | 23 | 47 | 0.975216 | 80 | 84 | 41 | 41 | 0.931702 | 120 | 168 | 127 | 127 | 0.898647 |

Parsed from the arXiv v3 Table 8 text and sorted by δ. Ties between equal printed values are
broken by text order.

**q ≤ 30 subset** (rank: q; a/a⁻¹; δ):

* 1–15: the q | 24 races and (3;2), (4;3), as in §1a.
* 25: 15; 2/8; 0.983853
* 30: 15; 14; 0.979043
* 38: 20; 3/7; 0.976713
* 50: 20; 19; 0.966662
* 51: 15; 7/13; 0.964719
* 58: 16; 7; 0.959790
* 61: 15; 11; 0.955226
* 64: 5; 2/3; 0.952175
* 65: 20; 13/17; 0.948637
* 68: 16; 3/11; 0.947721
* 75: 16; 5/13; 0.938369
* 81: 20; 11; 0.931367
* 90: 16; 15; 0.926101
* 94: 21; 5/17; 0.923779
* 106: 21; 2/11; 0.911872
* (119: 28; 11/23; 0.898807, below 0.9)

Shchebetov (2018) selected exactly FM ranks 1–15 for his 10¹⁵ search (slides 14 and 27) [P].
For FM ranks 16+ with q ≤ 30 (moduli 5, 15, 16, 20, 21), **no first-crossing literature was
found**. These are presumably small; the census can compute them directly.

Related density sources (not first crossings):

* Rubinstein–Sarnak, Experiment. Math. 3 (1994) 173–197. δ for q = 3, 4, 5, 7, 11, 13
  (N vs R) [S: BFHR Table I col. B5; bib [215]].
* Feuerverger–Martin, Experiment. Math. 9 (2000) 535–570. All race densities for q ∈ {5, 8, 12};
  two- and three-way for q ∈ {7, 9}; two-way for q = 11 [S: bib [225]].
* Neither paper tabulates first crossings (per the annotated bibliography summaries; the
  originals were not read).

---

## 5. Recent work (2018–2026) and number fields

### 5a. A. Bailleul, M. Hayani, T. Untrau

*A Wasserstein metric approach to generalized Skewes' numbers. I. Prime number races*,
arXiv:2603.20093 (v1 20 Mar 2026, v2 24 Apr 2026; MPIM-Bonn-2026; 32 pp.) [P].

* **Definition.** x_{q;a,b} := inf{x ≥ 2 : π(x;q,a) > π(x;q,b)} with a a residue and b a
  non-residue (note the order), plus x_{q;R,NR} for the averaged R vs N race.
* **Thm 1.3.** Unconditional disproof of Fiorilli's Conj. 1.19 (see §3).
* **Thm 1.4.** Quantitative rate of convergence to δ_{q;a,b} under GRH + ELI_A(X_q), uniform
  for X ≥ (φ(q) log q)^{(Lφ(q) log q)^A}.
* **Thm 1.5.** log log log x_{q;a,b} ≪_A log φ(q).
* **Thm 1.6.** log log log x_{q;R,NR} ≪_A ϱ(q)/log rad(q) + log ϱ(q).
* There are **no numerical values and no explicit constants**.
* A function-field follow-up is announced ("will be detailed in a follow-up"). **No Part II
  was found** as of 2026-10-27.
* Related talks: Hayani (MPIM, Mar 2026); Bailleul (Apr 2026).
* Bailleul–Hayani, "Prime ideal races with several competitors", arXiv:2508.04087. Not read;
  appeared in a search.

### 5b. A. S. Shchebetov and S. D. Shchebetov (2017–2019)

Exhaustive prime-by-prime search to 10¹⁵ (and 5·10¹⁵ or 10¹⁶ for some races) for FM ranks
1–15.

* 25 new zones in 9 races.
* First zones found for (12,5), (12,7), (24,17), (24,19), and the exact first point of (24,13).
* None found for (8,3), (12,11), (24,5), (24,7), (24,11), (24,23).
* The data are in 18 OEIS sequences (A051024/5, A295353–6, A297005/6, A297354–7,
  A297447–50, A298820/1). Some names were corrected by J. Song in Oct 2026.
* Grey literature only: the EUCYS 2018 presentation and the ISEF abstract. A Math. Comp.
  paper was announced but not found.

### 5c. B. Chaffin (OEIS A275939 comment, Jun/Jul 2026)

* (24,23): first crossing at 7,390,203,183,507,483,529.
* (12,11): none ≤ 10¹⁹.
* This is the only known search to 10¹⁹ for mod 12/24 races. It is undocumented apart from the
  OEIS comment. It may or may not cover (8,3), (24,5), (24,7), (24,11); **UNVERIFIED** whether
  those were also searched.

### 5d. J. Song (OEIS 2018–2026)

* Real-character team races: A321857–A321865, A071838, A320857, A306499, A306500, A392284,
  A326615 (Kronecker-symbol sums), A329224.
* The nine index-4 coset races mod 24 (A398234–A401245), searched to ≈1–2.3·10¹².

### 5e. Earlier explicit and conditional work referenced above

* Schlage-Puchta 2004 (GRH, iterated-exponential bounds).
* Sneed (preprint, c. 2008–2010; local copy): proves infinitely many lead changes for all
  races with q ≤ 100 via Ingham-type numerics with lcalc/Rumely zeros. It gives **no
  location** of first crossings.
* Ford–Konyagin survey (2002; arXiv:1910.08983) poses the first-sign-change problem.
* Problem list arXiv:2407.03530:
  * Problem 7: bound the first sign change of π(x;q,a) − π(x;q,b) in terms of q.
  * Problem 9: Skewes numbers of highly biased races; Fiorilli's conjecture.
  * Problem 8: localize the first negative value of Δ_K(x) for Mertens' third theorem over
    Q(√5), Q(√13); Büthe gives < 2·10²¹⁵ for Q.

### 5f. Number fields and Chebotarev races

**No explicit first-sign-change computation or certificate was found.** Related work is all
densities or biases:

* Fiorilli–Jouve, J. Inst. Math. Jussieu 23 (2024), and J. Éc. polytech. 9 (2022);
* Bailleul, Algebra Number Theory 15 (2021);
* Aoki–Koyama, JNT 245 (2023);
* Hathi–Lee, arXiv:2112.02166, for Mertens' third theorem (densities and Cramér's
  inequality).

Problem 8 of the problem list states the localization is open. The annotated bibliography
(which covers literature to ~2023) lists no such computation.

---

## 6. Gaps and cautions for the census paper

1. **Two-way races with a residue b ≠ 1, or q ∤ 24.** Examples: (5;2,4), (7;3,2), all of
   q = 15, 16, 20, 21. No tabulation exists except "1 vs −1" (A275939). These are presumably
   tiny and need no citation, but the census should compute and state them.
2. **Open q | 24 races.**
   * (8,3): ≤ 10¹⁵ exhaustive, nothing found. The FH certified window is 6.68–6.96·10¹⁸.
     An exhaustive or Meissel search to 10¹⁹ would settle it, and none has been reported.
   * (12,11): none ≤ 10¹⁹ (OEIS only). FH certified ≈ 10^81.45.
   * (24,5), (24,7), (24,11): ≤ 10¹⁵ only. FH certified 10^70.0, 10^40.2, 10^353.0.
   * Only FH 2001 gives rigorous bounds. Any certificate for b ≠ 1, q ∤ 24, or team races
     would be new.
3. **Full team races with ρ(q) ≥ 4** (q = 8, 12, 15, 16, 20, 21, 24, 28 in Fiorilli's
   normalisation). No crossing or certificate is known. BFHR give only a heuristic plot for
   q = 8 near 10^28.127, explicitly "possibly the first region".
4. **(24,23) = 7,390,203,183,507,483,529** rests on a single OEIS comment. If cited, cite it as
   "B. Chaffin, OEIS A275939 (2026)". The census should re-verify independently, for example
   with a Meissel/LMO computation of π(x;24,1) and π(x;24,23) at that x and an argument that
   no earlier x exists. Exhaustiveness to 7.39·10¹⁸ is not documented.
5. **The Shchebetov data are grey literature** (OEIS b-files and slides). The 10¹⁵ limit is
   documented. Higher limits (5·10¹⁵, 10¹⁶) are only partially documented.
6. **Deléglise–Dusart–Roblot.** Only the author preprint was read (Table 1 there). The
   published Math. Comp. version was not accessible. Their mod 4 region (9.97·10¹⁷ to
   1.005·10¹⁸) was not exhaustively located, and they give no other moduli.
7. **Not read in the original** (taken via the annotated bibliography, BFHR, FH,
   Granville–Martin, Shchebetov): Leech 1957, Bays–Hudson 1977/78/79/83/2000,
   Rubinstein–Sarnak 1994, Feuerverger–Martin 2000.
8. **FH Cor. 4 rounding.** 5·10¹⁹ and 10³⁵³ are slightly below the rigorous window ends
   (5.06·10¹⁹ and 1.025·10³⁵³). The annotated bibliography misprints 10³⁵³ as 10²⁵³.
9. **Granville–Martin vs BFHR/Shchebetov** disagree on the mod 4 zone 5 and 6 starts
   (…697/…217 vs …709/…293). This is not relevant to first crossings.
