#!/usr/bin/env python3
"""
Rigorous Arb certification of a negative value of the error term in Mertens' first theorem,
    Delta_1(x) = sum_{p <= x} log p / p - log x - E,   E = -gamma - sum_p log p/(p(p-1)) = -1.33258...,
with the two-scale kernel K = K_core * psi_delta (THEOREM.md). With x = e^u, F(u) := sqrt(x) Delta_1(x) satisfies
    I(omega) := int K(u - omega) F(u) du  <=  U := 1 + 3 eps + S_M(omega, T0) + R4 + R5 + R_small,
    S_M(omega, T0) = -2 sum_{0 < gamma <= T0} hat{K}(gamma) Re(e^{i gamma omega} / (rho - 1)),
where eps = 1.5423e-9 (Platt-Trudgian 2016, (9): |theta(y) - y| <= eps y for y >= e^200),
    R4 = bound_R4_envelope / bound_R4_stieltjes (error_bounds.py; their majorant bounds 2 sum |hat kappa_core(eta1 gamma)|/gamma,
         and |rho - 1| >= |gamma| for zeros on the critical line),
    R5 = 4 C_M log(A) e^{u_hi/2} (M pi/(A delta))^{2M}   (2M-fold integration by parts for |gamma| > A = 3e12),
    R_small = log(2 pi) e^{-u_lo/2} + e^{-5 u_lo/2} + 1.53 e^{-u_lo/6} + 1.36 (u_hi/log 2 + 2) e^{-u_lo/4}.
If U < 0, then Delta_1(x) < 0 for some x with log x in [omega - eta0, omega + eta0].

Zeros: Arb balls from the Mertens / Skewes zero caches (n <= 13,955,181 text caches; beyond, z_{n0:09d}.npz batches).
Usage: certify_first.py OMEGA ETA1 DELTA M LAM T0 [WORKERS]     e.g. LAM = 4,52,44 (percent weights)
"""
import glob
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import flint
from flint import acb, arb

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "certification"))
from certio import CertificationError, exact_int, fail, provenance, require
from error_bounds import A_PT, bound_R4_envelope, bound_R4_stieltjes

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))   # skewes_bound/
CERT_DIR = os.path.join(REPO, "certification")
MERTENS_ZCACHE = os.environ.get("MERTENS_ZERO_CACHE", os.path.join(REPO, "..", "mertens_bound", "certification", "zero_cache"))
EXT_ZCACHE = os.path.join(REPO, "zero_cache_ext")
ARB_ZERO_DIR = os.environ.get("ARB_ZERO_DIR", os.path.join(REPO, "zero_cache_arb"))
CODE = [os.path.abspath(__file__), os.path.join(CERT_DIR, "certio.py"), os.path.join(CERT_DIR, "error_bounds.py")]
N_LOW = 13955181    # N(6.8e6): last zero in the Arb text caches
BATCH = 2000        # batch size of the Arb zero files beyond N_LOW
DPS = 30


def khat(xi, m):
    """hat{kappa}_{2m}(xi) = sin(xi)/xi * prod_k (k pi)^2 / ((k pi)^2 - xi^2), with the removable singularities."""
    pi = arb.pi()
    for k in range(1, m + 1):
        if (xi - k * pi).contains(0):
            u = (xi - k * pi) / pi
            res = (-1) ** (k + 1) * u.sinc_pi() * (k * pi) ** 2 / (xi * (k * pi + xi))
            for j in range(1, m + 1):
                if j != k:
                    res *= (j * pi) ** 2 / ((j * pi) ** 2 - xi * xi)
            return res
    res = xi.sin() / xi
    for k in range(1, m + 1):
        res *= (k * pi) ** 2 / ((k * pi) ** 2 - xi * xi)
    return res


def partial_sum(gs, P):
    flint.ctx.dps = DPS
    om, e1, de, T0 = arb(P["omega"]), arb(P["eta1"]), arb(P["delta"]), arb(P["T0"])
    l = [arb(x) / 100 for x in P["lam"]]
    tot = arb(0)
    for g in gs:
        require(g < T0, "a loaded zero ball is not strictly below T0")   # every summed zero must lie below T0
        w = (l[0] * khat(e1 * g, 1) + l[1] * khat(e1 * g, 2) + l[2] * khat(e1 * g, 3)) * khat(de * g, P["M"])
        tot += -2 * w * (acb(0, g * om).exp() / acb(arb("-0.5"), g)).real      # rho - 1 = -1/2 + i gamma
    return tot


def summarise(gs, P):
    require(all(a < b for a, b in zip(gs, gs[1:])), "zeros not strictly increasing within a batch")
    return (partial_sum(gs, P).str(30, radius=True), gs[0].str(30, radius=True), gs[-1].str(30, radius=True), len(gs))


def low_job(args):
    path, k, P = args
    flint.ctx.dps = DPS
    gs = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i >= k:
                break
            gs.append(arb(line.split("\t")[0].strip()))
    require(len(gs) == k, f"{path}: expected {k} zeros")
    return summarise(gs, P)


def high_job(args):
    path, P = args
    flint.ctx.dps = DPS
    z = np.load(path)
    gs = [arb(float(h)) + arb(float(l)) + arb(0, float(r)) for h, l, r in zip(z["hi"], z["lo"], z["rad"])]
    return summarise(gs, P)


def low_tasks(n_end):
    best = {}
    for d in (MERTENS_ZCACHE, EXT_ZCACHE):
        for p in glob.glob(os.path.join(d, "*.txt")):
            n0, k, dps = map(int, os.path.basename(p)[:-4].split("_"))
            if n0 not in best or (k, dps) > best[n0][:2]:
                best[n0] = (k, dps, p)
    out, cur = [], 1
    while cur <= n_end:
        require(cur in best, f"missing Arb cache batch starting at n={cur}")
        k, _, p = best[cur]
        k = min(k, n_end + 1 - cur)
        out.append((p, k))
        cur += k
    return out


def high_tasks(n_end):
    out, cur = [], N_LOW + 1
    while cur <= n_end:
        p = os.path.join(ARB_ZERO_DIR, f"z_{cur:09d}.npz")
        require(os.path.exists(p), f"missing Arb zero batch {p}")
        out.append(p)
        cur += BATCH
    return out


def main():
    if len(sys.argv) < 7:
        fail(__doc__)
    omega_s, eta1_s, delta_s, M_s, lam_s, T0_s = sys.argv[1:7]
    workers = int(sys.argv[7]) if len(sys.argv) > 7 else 48
    lam = [int(x) for x in lam_s.split(",")]
    require(len(lam) == 3 and all(x >= 0 for x in lam) and sum(lam) == 100, "LAM must be three percentages summing to 100")
    M = int(M_s)
    P = {"omega": omega_s, "eta1": eta1_s, "delta": delta_s, "M": M, "lam": lam, "T0": T0_s}

    flint.ctx.dps = DPS
    T0 = arb(T0_s)
    N = exact_int(T0.zeta_nzeros(), "N(T0)")
    nb = acb.zeta_zeros(N, 2)
    require(nb[0].imag < T0 < nb[1].imag, "T0 not strictly between gamma_N and gamma_{N+1}")

    om, e1, de = arb(omega_s), arb(eta1_s), arb(delta_s)
    eta0 = e1 + de
    u_lo, u_hi = om - eta0, om + eta0
    require(u_lo >= 400, "omega - eta0 >= 400 is required for (9) of Platt-Trudgian (2016) at sqrt(x)")
    require(e1 * T0 > 3 * arb.pi(), "eta1 * T0 > 3 pi is required by the R4 bounds")
    require(de < e1, "mollifier must be narrower than the core")

    t0 = time.time()
    lt = low_tasks(min(N, N_LOW))
    ht = high_tasks(N) if N > N_LOW else []
    with Pool(workers) as pool:
        res = pool.map(low_job, [(p, k, P) for p, k in lt], chunksize=4)
        res += pool.map(high_job, [(p, P) for p in ht], chunksize=4)

    flint.ctx.dps = DPS
    for a, b in zip(res, res[1:]):
        require(arb(a[2]) < arb(b[1]), "zero batches not increasing across a boundary")
    require(sum(r[3] for r in res) == N, "number of zeros summed differs from N(T0)")
    S = arb(0)
    for r in res:
        S += arb(r[0])

    l = [arb(x) / 100 for x in lam]
    r4_env = bound_R4_envelope(om, e1, T0, *l)
    r4_st = bound_R4_stieltjes(om, e1, T0, *l)
    r4 = r4_env if r4_env < r4_st else r4_st
    C_M = arb(4 ** M * math.factorial(M) ** 2) / arb(2 * math.factorial(2 * M))
    rate = (M * arb.pi()) / (A_PT * de)
    require(rate < 1, "R5 rate must be < 1")
    r5 = (4 * C_M * A_PT.log() * (u_hi / 2).exp() * rate ** (2 * M)).upper()
    r_small = ((2 * arb.pi()).log() * (-u_lo / 2).exp() + (-5 * u_lo / 2).exp() + arb("1.53") * (-u_lo / 6).exp()
               + arb("1.36") * (u_hi / arb(2).log() + 2) * (-u_lo / 4).exp()).upper()
    eps = arb("1.5423e-9")
    U = 1 + 3 * eps + S + r4 + r5 + r_small

    cert = provenance(P, CODE)
    cert.update({
        "omega": omega_s, "eta1": eta1_s, "delta": delta_s, "M_moll": M, "lam": lam_s, "T0": T0_s, "N_T0": N,
        "S_M": S.str(20, radius=True), "R4_envelope": r4_env.str(10, radius=False), "R4_stieltjes": r4_st.str(10, radius=False),
        "R5": r5.str(6, radius=False), "R_small": r_small.str(6, radius=False),
        "U": U.str(15, radius=True), "U_upper": U.upper().str(10, radius=False),
        "u_lo": u_lo.lower().str(16, radius=False), "u_hi": u_hi.upper().str(16, radius=False),
        "log10_x_upper": (u_hi / arb(10).log()).upper().str(12, radius=False),
        "elapsed_s": round(time.time() - t0, 1),
    })
    print(json.dumps(cert, indent=2))
    require(U < 0, "U is not certified negative")
    print(f"CERTIFIED: Delta_1(x) < 0 for some x in [exp({cert['u_lo']}), exp({cert['u_hi']})], so X < 10^{cert['log10_x_upper']}")
    with open(os.path.join(HERE, f"first_cert_T{T0_s}.json"), "w") as f:
        json.dump(cert, f, indent=2)


if __name__ == "__main__":
    try:
        main()
    except CertificationError as e:
        fail(str(e))
