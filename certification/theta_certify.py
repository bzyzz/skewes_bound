#!/usr/bin/env python3
"""
Rigorous Arb ball-arithmetic certification of theta(x) > x with the two-scale kernel
K = K_core * psi_delta (paper, Theorem "theta criterion"):

    I_theta(omega, eta0) >= -1 - 1.5423e-9 + S_theta(omega, T0) - R3_theta - R4 - R5_theta,
    S_theta(omega, T0)    = -2 * sum_{0 < gamma <= T0} hat{K}(gamma) * Re(e^{i gamma omega} / rho),

using Platt-Trudgian (2016), (6) and (9), for the theta(x^{1/2}), log(2 pi) and x^{1/3} terms, Platt-Trudgian (2021)
for 0 < |gamma| <= A = 3e12, and 2M-fold integration by parts for |gamma| > A:
    R3_theta = e^{-u_lo/2} log(2 pi) + 3 e^{-u_lo/6},
    R5_theta = 4 C_M log(A) e^{u_hi/2} (M pi / (A delta))^{2M}.
R4 is bound_R4_envelope / bound_R4_stieltjes from error_bounds.py; their majorant carries an extra 1/(t^2 omega)
term, so they also bound the theta tail 2 sum_{T0 < gamma <= A} |hat{kappa}_core(eta1 gamma)| / gamma.

Zeros: n <= 13,955,181 from the Arb caches used by certify_S.py; n > 13,955,181 from Arb batches
z_{n0:09d}.npz (fields hi, lo, rad: the ball hi + lo +/- rad encloses gamma_n) produced by acb.zeta_zeros,
in the directory $ARB_ZERO_DIR (default: ../zero_cache_arb).

Usage: theta_certify.py OMEGA ETA1 DELTA M LAM T0 [WORKERS]     e.g. LAM = 4,52,44 (percent weights)
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

from certio import CertificationError, ball, exact_int, fail, provenance, require
from error_bounds import A_PT, N_SUB, bound_R4_envelope, bound_R4_stieltjes

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MERTENS_ZCACHE = os.environ.get("MERTENS_ZERO_CACHE", os.path.join(ROOT, "..", "mertens_bound", "certification", "zero_cache"))
EXT_ZCACHE = os.path.join(ROOT, "zero_cache_ext")
ARB_ZERO_DIR = os.environ.get("ARB_ZERO_DIR", os.path.join(ROOT, "zero_cache_arb"))
CODE = [os.path.abspath(__file__), os.path.join(HERE, "certio.py"), os.path.join(HERE, "error_bounds.py")]
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
    used = 0
    for g in gs:
        if g > T0:
            break
        require(not g.overlaps(T0), "a zero ball overlaps T0")
        w = (l[0] * khat(e1 * g, 1) + l[1] * khat(e1 * g, 2) + l[2] * khat(e1 * g, 3)) * khat(de * g, P["M"])
        tot += -2 * w * (acb(0, g * om).exp() / acb(arb("0.5"), g)).real
        used += 1
    # the zeros are increasing, so everything after the first ball above T0 must also be above T0
    require(all(g > T0 for g in gs[used:]), "a zero ball above the cut-off is not certified above T0")
    return tot, used


def summarise(gs, P):
    """(sum over the zeros <= T0, first ball, last ball, number loaded, number summed, last summed ball)."""
    require(len(gs) > 0 and all(g.is_finite() and g > 0 for g in gs), "empty batch or non-finite / non-positive zero ball")
    require(all(a < b for a, b in zip(gs, gs[1:])), "zeros not strictly increasing within a batch")
    tot, used = partial_sum(gs, P)
    return (ball(tot, 40), ball(gs[0], 30), ball(gs[-1], 30), len(gs), used, ball(gs[used - 1], 30) if used else None)


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
    path, need, P = args
    flint.ctx.dps = DPS
    with np.load(path, allow_pickle=False) as z:
        hi, lo, rad = z["hi"].copy(), z["lo"].copy(), z["rad"].copy()
    require(hi.ndim == 1 and hi.shape == lo.shape == rad.shape, f"{path}: malformed zero batch")
    require(bool(np.all(np.isfinite(hi)) and np.all(np.isfinite(lo)) and np.all(np.isfinite(rad)) and np.all(rad >= 0)),
            f"{path}: non-finite entries or negative radii")
    # the file z_{n0}.npz must hold the zeros n0, ..., n0 + BATCH - 1 (the last file needed may hold fewer)
    require(len(hi) == BATCH if need is None else need <= len(hi) <= BATCH, f"{path}: unexpected number of zeros {len(hi)}")
    gs = [arb(float(h)) + arb(float(l)) + arb(0, float(r)) for h, l, r in zip(hi, lo, rad)]
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
        # every file but the last must be full; the last must reach the zero n_end
        out.append((p, None if cur + BATCH <= n_end else n_end + 1 - cur))
        cur += BATCH
    return out


def main():
    if len(sys.argv) < 7:
        fail(__doc__)
    omega_s, eta1_s, delta_s, M_s, lam_s, T0_s = sys.argv[1:7]
    workers = int(sys.argv[7]) if len(sys.argv) > 7 else 16
    require(workers >= 1, "workers must be positive")
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
    require(u_lo >= 400, "omega - eta0 >= 400 is required for (9) of Platt-Trudgian (2016)")
    require(e1 * T0 > 3 * arb.pi(), "eta1 * T0 > 3 pi is required by the R4 bounds")
    require(de < e1, "mollifier must be narrower than the core")

    t0 = time.time()
    lt = low_tasks(min(N, N_LOW))
    ht = high_tasks(N) if N > N_LOW else []
    with Pool(workers) as pool:
        res = pool.map(low_job, [(p, k, P) for p, k in lt], chunksize=4)
        res += pool.map(high_job, [(p, need, P) for p, need in ht], chunksize=4)

    flint.ctx.dps = DPS
    for a, b in zip(res, res[1:]):
        require(arb(a[2]) < arb(b[1]), "zero batches not increasing across a boundary")
    # Completeness (consumer side). The batches are consecutive and increasing, so the summed zeros are the first
    # sum(r[4]) of them. This must equal N(T0) exactly, and the last summed ball must be the zero gamma_N that Arb
    # isolates independently here; a missing, truncated or repeated batch fails one of these checks.
    require(sum(r[3] for r in res) >= N, "fewer zeros loaded than N(T0)")
    n_used = sum(r[4] for r in res)
    require(n_used == N, f"summed {n_used} zeros, but N(T0) = {N}")
    k_last = max(i for i, r in enumerate(res) if r[4] > 0)
    require(all(r[4] == r[3] for r in res[:k_last]) and all(r[4] == 0 for r in res[k_last + 1:]),
            "the summed zeros are not an initial segment of the batches")
    require(arb(res[k_last][5]).overlaps(nb[0].imag), "the last summed zero is not gamma_N")
    require(arb(res[0][1]).overlaps(acb.zeta_zeros(1, 1)[0].imag), "the first zero is not gamma_1")
    S = arb(0)
    for r in res:
        S += arb(r[0])

    l = [arb(x) / 100 for x in lam]
    r4_env = bound_R4_envelope(om, e1, T0, *l, n_sub=N_SUB)
    r4_st = bound_R4_stieltjes(om, e1, T0, *l, n_sub=N_SUB)
    C_M = arb(4 ** M * math.factorial(M) ** 2) / arb(2 * math.factorial(2 * M))
    rate = (M * arb.pi()) / (A_PT * de)
    require(rate < 1, "R5 rate must be < 1")
    r5 = (4 * C_M * A_PT.log() * (u_hi / 2).exp() * rate ** (2 * M)).upper()
    r3 = ((-u_lo / 2).exp() * (2 * arb.pi()).log() + 3 * (-u_lo / 6).exp()).upper()
    # Serialise the budget first and take both margins from the strings that are written, so that a reader of the
    # certificate who recomputes the margin from its entries gets the same decision.
    saved = {"S_theta": ball(S, 30), "R3_theta": ball(r3, 30), "R4_envelope": ball(r4_env, 30),
             "R4_stieltjes": ball(r4_st, 30), "R5_theta": ball(r5, 30)}
    base = -1 - arb("1.5423e-9") + arb(saved["S_theta"]) - arb(saved["R3_theta"]) - arb(saved["R5_theta"])
    I_env, I_st = (base - arb(saved["R4_envelope"])).lower(), (base - arb(saved["R4_stieltjes"])).lower()

    cert = provenance(P, CODE)
    cert.update({
        "omega": omega_s, "eta1": eta1_s, "delta": delta_s, "M_moll": M, "lam": lam_s, "T0": T0_s, "N_T0": N,
        "N_summed": n_used, "tail_subintervals": N_SUB, **saved,
        "I_lower_envelope": ball(I_env, 30), "I_lower_stieltjes": ball(I_st, 30),
        "u_lo": ball(u_lo.lower(), 30), "u_hi": ball(u_hi.upper(), 30),
        "elapsed_s": round(time.time() - t0, 1),
    })
    print(json.dumps(cert, indent=2))
    best = I_env if I_env > I_st else I_st
    require(best > 0, "theta margin is not certified positive")
    print(f"CERTIFIED: theta(x) > x for some x in [exp({cert['u_lo']}), exp({cert['u_hi']})]")
    with open(os.path.join(HERE, f"theta_cert_T{T0_s}.json"), "w") as f:
        json.dump(cert, f, indent=2)


if __name__ == "__main__":
    try:
        main()
    except CertificationError as e:
        fail(str(e))
