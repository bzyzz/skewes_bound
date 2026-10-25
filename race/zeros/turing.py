#!/usr/bin/env python3
"""Turing's method for real primitive chi: certify N(T, chi) = F(T) at checkpoints T.

Setting. chi is real primitive mod Q, so W_chi = 1. theta(t) = (t/2) log(Q/pi) + Im log Gamma((1/2+delta+it)/2),
and Z(t) = e^{i theta(t)} L(1/2+it, chi) is real. N(T) counts zeros rho = beta+i gamma, 0 < gamma <= T, with
multiplicity. S(t) = (1/pi) arg L(1/2+it), by continuous variation along 2 -> 2+it -> 1/2+it.

(1) If L(sigma, chi) != 0 for sigma in [1/2, 2], then N(T) = theta(T)/pi + S(T) for T not an ordinate.
    Proof: apply the argument principle to xi(s, chi) = (Q/pi)^{s/2} Gamma((s+delta)/2) L(s, chi) on
    [-1, 2] x [0, T]. Since xi(s) = xi(1-s) and xi(conj s) = conj xi(s), the change in argument is twice
    the change along 1/2 -> 2 -> 2+iT -> 1/2+iT. The real segment contributes 0 because xi is real and
    nonzero there. Hypothesis (1) is checked here in Arb: L > 0 on a cover of [1/2, 2].
(2) Trudgian (Math. Comp. 80 (2011), Thm 3.3): |int_{t1}^{t2} S(t,chi) dt| <= 1.975 + 0.084 log(Q t2/(2 pi))
    for t2 > t1 > 50.
(3) Let F(t) = #{certified disjoint sign-change intervals of Z lying in (0, t]}, so N(t) >= F(t). Suppose
    N(T) >= F(T) + 1. Then for t in [T, T+h], N(t) >= F(t) + 1, hence S(t) >= F(t) + 1 - theta(t)/pi.
    If L := int_T^{T+h} (F + 1 - theta/pi) > bound (2), this is a contradiction. So N(T) = F(T), and all
    zeros with 0 < gamma <= T are simple and on the critical line.
    The integral is bounded below rigorously: int F >= sum_i (T+h - max(hi_i, T)) over the intervals with
    hi_i <= T+h; int theta is bounded above by an upper Riemann sum, since theta' > 0 on [T, T+h]
    (checked in Arb on balls).

usage: turing.py certified.npz h T1,T2,...   -> prints the results and writes turing_<q>_<n>.json next to the input
"""
import json
import os
import sys

import numpy as np
import flint
from flint import acb, arb, dirichlet_char

flint.ctx.dps = 30


def setup(q, n):
    chi = dirichlet_char(q, n)
    delta = int(chi.parity())
    lq = (arb(q)/arb.pi()).log()
    sh = (arb(1)/2 + delta)/2
    theta = lambda t: t/2*lq + acb(sh, t/2).lgamma().imag
    dtheta = lambda t: lq/2 + acb(sh, t/2).digamma().real/2     # valid for ball t
    return chi, theta, dtheta


def L_positive_on_half_to_two(chi, pieces=3000, r=arb("0.01"), R=arb("0.1"), arcs=400):
    """Check L(sigma, chi) > 0 for all sigma in [1/2, 2].

    Arb evaluates L through Hurwitz zeta functions, which have a pole at s = 1, so balls must avoid 1.
    [1/2, 1-r] and [1+r, 2] are covered by balls. For |sigma - 1| <= r we use
    |L(sigma) - L(1)| <= r * R * max_{|z-1|=R} |L(z)| / (R - r)^2 (Cauchy's formula for L'),
    with the circle covered by balls that stay at distance > R/2 from 1.
    """
    def positive(a, b, depth=0):
        mid = (a + b)/2
        ball = arb(mid.mid(), ((b - a)/2).upper()*1.0000001)
        if acb.dirichlet_l(acb(ball), chi).real > 0:
            return True
        if depth > 30:
            return False
        return positive(a, mid, depth + 1) and positive(mid, b, depth + 1)

    for a, b in ((arb(1)/2, 1 - r), (1 + r, arb(2))):
        w = (b - a)/(pieces//2)
        for k in range(pieces//2):
            if not positive(a + k*w, a + (k + 1)*w):
                return False, float((a + (k + arb(1)/2)*w).mid())
    rad = (arb.pi()*R/arcs).upper()*1.05
    mx = arb(0)
    for k in range(arcs):
        z = 1 + R*acb(0, 2*arb.pi()*k/arcs).exp()
        zb = acb(arb(z.real.mid(), rad), arb(z.imag.mid(), rad))
        mx = max(mx, abs(acb.dirichlet_l(zb, chi)).upper())
    near = acb.dirichlet_l(acb(1), chi).real - r*R*mx/(R - r)**2
    if not (near > 0):
        return False, 1.0
    return True, None


def turing_at(theta, dtheta, lo, hi, q, T, h, step=0.05):
    T = arb(T); Th = T + h
    Tf, Thf = float(T.mid()), float(Th.mid())
    inside = np.nonzero((lo < Tf) & (hi > Tf))[0]
    if len(inside):
        return dict(ok=False, why=f"T inside a certified interval {inside}")
    F_T = int(np.sum(hi <= Tf))
    # int_T^{T+h} F dt (lower bound)
    sel = np.nonzero(hi <= Thf)[0]
    intF = arb(0)
    for i in sel:
        intF += Th - max(arb(float(hi[i])), T)
    # theta' > 0 on [T, T+h]
    nb = int(h/1.0) + 1
    for k in range(nb):
        c = T + (k + arb(1)/2)*h/nb
        if not (dtheta(arb(c.mid(), (arb(h)/nb/2).upper()*1.0000001)) > 0):
            return dict(ok=False, why="theta' not certified positive")
    m = int(round(h/step))
    up = arb(0)
    for k in range(1, m + 1):
        up += theta(T + k*arb(h)/m)
    int_theta_up = (up*arb(h)/m).upper()
    lower = intF + h - int_theta_up/arb.pi()
    # Trudgian (2011) Thm 3.3, and Rumely's (1993) original constants as recomputed by Trudgian (2011) §3.3;
    # take the larger, so the result does not rest on either alone.
    ell = (q*Th/(2*arb.pi())).log()
    bound = max(arb("1.975") + arb("0.084")*ell, arb("1.8397") + arb("0.1242")*ell)
    ok = bool(lower > bound)
    return dict(ok=ok, T=Tf, h=float(h), F_T=F_T, lower=float(lower.lower()), bound=float(bound.upper()),
                theta_over_pi_T=float((theta(T)/arb.pi()).mid()))


def main():
    path, h, Ts = sys.argv[1], float(sys.argv[2]), [float(x) for x in sys.argv[3].split(",")]
    d = np.load(path)
    q, n = int(d["q"]), int(d["n"])
    if len(d["fail"]):
        print(f"WARNING: {len(d['fail'])} uncertified located zeros; they are not counted in F")
    lo, hi = d["lo"], d["hi"]
    good = ~np.isnan(lo)
    lo, hi = lo[good], hi[good]
    assert np.all(hi[:-1] < lo[1:]) and np.all(lo < hi)
    chi, theta, dtheta = setup(q, n)
    okL, where = L_positive_on_half_to_two(chi)
    print(f"chi = {q}.{n}: L(sigma) > 0 on [1/2, 2]: {okL}" + ("" if okL else f" (fails near {where})"))
    res = dict(q=q, n=n, L_positive_half_to_two=okL, eps=float(d["eps"]), checkpoints=[])
    for T in Ts:
        r = turing_at(theta, dtheta, lo, hi, q, T, h)
        res["checkpoints"].append(r)
        print(f"  T={T:.0f}: " + (f"N(T) = F(T) = {r['F_T']} (theta/pi = {r['theta_over_pi_T']:.2f}); "
              f"Turing lower {r['lower']:.3f} > bound {r['bound']:.3f}" if r.get("ok") else f"FAILED {r}"), flush=True)
    res["all_ok"] = okL and all(r.get("ok") for r in res["checkpoints"])
    out = os.path.join(os.path.dirname(path), f"turing_{q}_{n}.json")
    json.dump(res, open(out, "w"), indent=1)
    print("all_ok =", res["all_ok"], "->", out)


if __name__ == "__main__":
    main()
