# Statements, inputs and certified numbers

This file summarises what is proved, what each proof relies on, and where every number comes from. Each input is marked **[cited]** (published source), **[derived]** (proved in the accompanying manuscript, in preparation; the proofs are outlined here) or **[computed]** (by a script in [`certification/`](certification/)). The prime-race result has its own file, [`race/THEOREM.md`](race/THEOREM.md).

Notation: $\Xi := \min\{x \ge 2 : \pi(x) > \mathrm{li}(x)\}$; nontrivial zeros of $\zeta(s)$ are $\rho = \beta + i\gamma$; $A = 3\cdot 10^{12}$.

## Status

| Result | Status |
|---|---|
| Main theorem: $\Xi < e^{727.951335596}$ ($T_0 = 8\cdot 10^7$) | **certified in Arb**: all 195,565,267 zeros isolated with Arb (`certification/arb_recert.py`, `arb_combine.py`, `s_certificate_arb.json`) |
| $\Xi < e^{727.951337382}$ ($T_0 = 6.8\cdot 10^6$) | certified in Arb |
| $\theta(x) > x$ for some $x \in [e^{727.9513323698}, e^{727.9513328538}]$ ($T_0 = 5\cdot 10^7$) | certified in Arb (zeros above $6.8\cdot 10^6$ from the Arb re-isolation) |
| Smith refinement: $\Xi < e^{727.951335402030}$ | conditional on one value tabulated in Smith (2016), Table 4 |

## 1. The two-scale kernel [derived]

$$K = K_{\mathrm{core}} * \psi_\delta, \qquad K_{\mathrm{core}}(y) = \frac{1}{\eta_1}\sum_{m=1}^{3}\lambda_m \kappa_{2m}\!\left(\frac{y}{\eta_1}\right), \qquad \psi_\delta(y) = \frac{1}{\delta}\kappa_{2M}\!\left(\frac{y}{\delta}\right),$$

with $\kappa_{2m}(x) = C_m\cos^{2m}(\pi x/2)\mathbf{1}_{[-1,1]}(x)$, $C_m = 4^m(m!)^2/(2(2m)!)$. Then $K \ge 0$, $\mathrm{supp} K = [-\eta_0, \eta_0]$ with $\eta_0 = \eta_1 + \delta$, $\int K = 1$, and $\hat K(\gamma) = \hat\kappa_{\mathrm{core}}(\eta_1\gamma)\,\hat\kappa_{2M}(\delta\gamma)$ with
$\hat\kappa_{2m}(\xi) = \frac{\sin\xi}{\xi}\prod_{k=1}^m \frac{(k\pi)^2}{(k\pi)^2 - \xi^2}$.

- Compact support with unit mass: no truncation error and no exterior tail.
- $K$ vanishes to order $2M$ at $\pm\eta_0$, so $2M$-fold integration by parts bounds the zeros with $|\gamma| > A$ (possibly off the critical line) by a constant times $e^{\omega/2}(M\pi/(A\delta))^{2M}$. This uses only $A \gg T_0$.
- The core only has to control the on-line zeros $T_0 < |\gamma| \le A$; this needs $\eta_1 T_0 \approx 10$–$12$.
- Limit: beating $e^{-\omega/2}$ needs $A\delta \gtrsim 1600$, so $\eta_0 \gtrsim 10/T_0 + 1600/A$.

For comparison, the truncated Gaussian of Lehman's theorem and its refinements needs $2A/\alpha \le \eta$ and $A^2/(2\alpha) > \omega/2 + O(\log A)$, so $\eta\sqrt\alpha \gtrsim 55$ before trimming. In practice $\eta_0 T_0 = 114$–$136$ (Saouter–Demichel 2010; Saouter–Trudgian–Demichel 2015). Logan-type compact kernels (Büthe 2016, 2018) control the high zeros by $1/\sinh c$ with $c \gtrsim \omega/2$, giving $\eta_0 T_0 \approx 70$ at our window.

## 2. Criterion for $\pi(x) - \mathrm{li}(x)$ [derived]

With $\Phi(u) = 1 + 2/u + 8/u^2 + 58.56/u^3$,

$$I(\omega, \eta_0) := \int_{\omega-\eta_0}^{\omega+\eta_0} K(u-\omega)\,\frac{u\{\pi(e^u) - \mathrm{li}(e^u)\}}{e^{u/2}\Phi(u)}\,\mathrm{d}u \;\ge\; J_K(\omega, T_0) - (R_2 + R_3 + R_4 + R_5),$$

$$J_K = -1 + \frac{S_K}{\Phi(\omega)}, \qquad S_K(\omega, T_0) = -\sum_{0 < |\gamma| \le T_0}\hat K(\gamma)\,e^{i\gamma\omega}\left(\frac{1}{\rho} + \frac{1}{\rho^2\omega}\right).$$

The error terms are implemented in `certification/error_bounds.py`:
- $R_2$: prime powers and $\log 2$;
- $R_3$: $u$-variation and the $\mathrm{li}$ remainder;
- $R_4$: on-line zeros $T_0 < |\gamma| \le A$, either the direct envelope from Lehman's Lemma 1 or the Stieltjes refinement;
- $R_5$: zeros with $|\gamma| > A$.

Inputs:
- [cited] Riemann–von Mangoldt explicit formula.
- [cited] Dusart (2018), Theorem 5.1: $\pi(y) \le \frac{y}{\log y}(1 + \frac{1}{\log y} + \frac{2}{\log^2 y} + \frac{7.32}{\log^3 y})$ for $y \ge 4\cdot 10^9$.
  - **Independent of Dusart's constant:** the original Ramanujan J. article was not accessible, but Axler, Integers 18 (2018) #A52, Theorem 2 (open access, read in the original) gives π(x) < x/(log x − 1 − 1/log x − 3.15/log²x − 12.85/log³x − 71.3/log⁴x − 463.2275/log⁵x − 4585/log⁶x) for x ≥ 49. `certification/check_axler_vs_dusart.py` certifies in Arb that this implies the Dusart-form bound for 60 ≤ log x ≤ 800, which covers every main term in this repository. (At log x below about 54 Axler's bound is the weaker one; there we rely on Dusart.)
- [cited] Platt–Trudgian (2021): all zeros with $0 < |\gamma| \le 3\cdot 10^{12}$ are simple and on the critical line.
- [cited] Lehman (1966), Lemmas 1–2; Saouter–Trudgian–Demichel (2015), Lemmas 2.1, 2.7–2.10.

## 3. Criterion for $\theta(x) - x$ [derived]

$$I_\theta := \int K(u-\omega)\,e^{-u/2}\{\theta(e^u) - e^u\}\,\mathrm{d}u \;\ge\; -1 - 1.5423\cdot 10^{-9} + S_\theta(\omega, T_0) - R_3^\theta - R_4 - R_5^\theta,$$

with $S_\theta = -\sum_{0 < |\gamma| \le T_0}\hat K(\gamma)e^{i\gamma\omega}/\rho$, $R_3^\theta = e^{-u_{\mathrm{lo}}/2}\log 2\pi + 3e^{-u_{\mathrm{lo}}/6}$ and $R_5^\theta = 4C_M\log A\,e^{u_{\mathrm{hi}}/2}(M\pi/(A\delta))^{2M}$.

Input: [cited] Platt–Trudgian (2016), (6) and (9), the latter being $|\theta(x) - x| \le 1.5423\cdot 10^{-9}x$ for $x \ge e^{200}$. Implemented in `certification/theta_certify.py`.

## 4. Conditional refinement of Smith (2016) [derived, conditional]

Assume Platt's Gaussian sum $S_1 + S_2 = 1.00288049783641096204 \pm 10^{-20}$ at $\omega = 727.951335401$, $\alpha = 4.5\cdot 10^{19}$, $T = 30{,}610{,}046{,}000$ (Smith 2016, Table 4). With the asymmetric kernel $G_\alpha\cdot(\mathbf{1}_{[-\eta_L, \eta_R]} * \psi_\delta)$, $\eta_L = 1.00\cdot 10^{-9}$, $\eta_R = 0.725\cdot 10^{-9}$, $\delta = 0.305\cdot 10^{-9}$, $M = 190$, the RH-dependent steps of Smith's Theorem 4.1 are replaced as follows:
- the zeros with $T < |\gamma| \le A$ are on the line by Platt–Trudgian (2021);
- $-\frac12\pi(x^{1/2})$ is bounded through the explicit formula for $\Pi(x^{1/2})$;
- the zeros with $|\gamma| > A$ are bounded by 380-fold integration by parts.

This gives $\Xi < e^{727.951335402030}$ unconditionally and $\Xi < e^{727.951335401708}$ under RH. Implemented in `certification/certify_smith_unconditional.py`.

## 5. Computed numbers [computed]

| Quantity | $T_0 = 5\cdot 10^7$ | $T_0 = 8\cdot 10^7$ (env.) | $T_0 = 8\cdot 10^7$ (St.) | $\theta$, $T_0 = 6.8\cdot 10^6$ | Smith, conditional |
|---|---|---|---|---|---|
| Certificate | `s_cert_T5e7.json` | `s_certificate.json` | `s_cert_T8e7_st.json` | `theta_cert_T6800000.json` | `smith_unconditional_certificate.json` |
| $N(T_0)$ | 118,488,123 | 195,565,267 | 195,565,267 | 13,955,181 | $1.04\cdot 10^{11}$ |
| $\omega$ | 727.951335479 | 727.951335444 | 727.951335441 | 727.9513329060 | 727.951335401 |
| $\eta_1$; $\delta$ ($M$) | $2.39\cdot 10^{-7}$; $2\cdot 10^{-9}$ (60) | $1.50\cdot 10^{-7}$; $2\cdot 10^{-9}$ (60) | $1.50\cdot 10^{-7}$; $2\cdot 10^{-9}$ (60) | $1.618\cdot 10^{-6}$; $2\cdot 10^{-9}$ (60) | see §4 |
| Margin, Arb-only zeros (`*_arb.json`) | $3.369\cdot 10^{-5}$ | $2.322\cdot 10^{-5}$ | $2.021\cdot 10^{-5}$ (St.) | $5.11\cdot 10^{-6}$ (St.) | $1.249\cdot 10^{-5}$ |
| Margin, float zeros above $6.8\cdot10^6$ (comparison) | $3.273\cdot 10^{-5}$ | $2.178\cdot 10^{-5}$ | $1.877\cdot 10^{-5}$ | — | — |
| Upper end $u_{\mathrm{hi}}$ | 727.951335720 | **727.951335596** | 727.951335593 | 727.951334526 | 727.951335402030 |

The Arb-only certificates (`certification/s_cert_T5e7_arb.json`, `s_certificate_arb.json`, `s_cert_T8e7_st_arb.json`) replace S_K by the sum over Arb-isolated zeros. `arb_combine.py` checks that the batches cover n = 1, …, N(8·10⁷) contiguously and in increasing order, and certifies N(T₀) with `arb_zeta_nzeros` together with γ_N < T₀ < γ_{N+1}. `check_theorem.py <file>` then verifies each certificate.

## 6. Prime number race modulo 24

See [`race/THEOREM.md`](race/THEOREM.md): π(x;24,11) < π(x;24,1) for some x with log x ∈ [111.5456338, 111.5460262], so the first sign change is below 10^48.4439 (Ford–Hudson 2001: 10^353).
