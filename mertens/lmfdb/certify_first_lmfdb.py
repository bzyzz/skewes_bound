#!/usr/bin/env python3
"""
Arb certification of a negative value of Delta_1 (Mertens' first theorem) near u = 495.7028331, using
  - our Arb zero balls for 1 <= n <= N_ARB = 195,565,267 (gamma < 8e7), as in ../criterion/certify_first.py;
  - D. J. Platt's zeros (LMFDB files) for N_ARB < n <= N(T0): gamma = t0 + Z 2^-101, stored to +-2^-102.
    Platt's isolation also proves that these zeros are simple and on the critical line, and his Turing-method
    check proves the list complete (Platt, Math. Comp. 86 (2017)); Platt-Trudgian (2021) independently gives
    beta = 1/2 for all of them.
The criterion is THEOREM.md section 4:
  U = 1 + 3 eps + S_M(omega, T0) + R4 + R5 + R_small < 0  implies  Delta_1(e^u) < 0 for some u in [omega - eta0, omega + eta0].
Checks: index contiguity across all sources, strictly increasing gamma, every summed gamma < T0, and
N(T0) = number summed (the next zero, gamma_{N+1}, must exceed T0).
usage: certify_first_lmfdb.py OMEGA ETA1 DELTA M LAM T0 FILELIST [WORKERS]
"""
import json, math, os, sys, time
from multiprocessing import Pool
import numpy as np
import flint
from flint import acb, arb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "criterion"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "certification"))
DATA = os.environ.get("LMFDB_ZETA_DIR", os.path.join(HERE, "data"))   # LMFDB/Platt files (download.sh)
import platt
from kernel_fast import Mollifier
from certio import provenance, require, fail, exact_int, CertificationError
from error_bounds import A_PT, bound_R4_envelope, bound_R4_stieltjes
import certify_first as cf            # zero-cache readers for n <= N_ARB

N_ARB = 195565267
DPS = 35
CODE = [os.path.abspath(__file__), os.path.join(HERE, "platt.py"), os.path.join(HERE, "platt.c"), os.path.join(HERE, "libplatt.so"),
        os.path.join(HERE, "kernel_fast.py"), os.path.join(HERE, "..", "criterion", "certify_first.py"),
        os.path.join(HERE, "..", "..", "certification", "error_bounds.py"), os.path.join(HERE, "..", "..", "certification", "certio.py")]
_moll = {}


def term_sum(gs, P):
    """sum over gamma in gs of -2 hatK(gamma) Re(e^{i gamma omega}/(rho - 1)), rho - 1 = -1/2 + i gamma."""
    flint.ctx.dps = DPS
    if P["M"] not in _moll:
        _moll[P["M"]] = Mollifier(P["M"])
    moll = _moll[P["M"]]
    om, e1, de, T0 = arb(P["omega"]), arb(P["eta1"]), arb(P["delta"]), arb(P["T0"])
    l = [arb(x)/100 for x in P["lam"]]
    tot = arb(0); prev = None
    for g in gs:
        require(g < T0, "a summed zero is not strictly below T0")
        if prev is not None:
            require(prev < g, "zeros not strictly increasing")
        prev = g
        core = sum((l[m - 1]*cf.khat(e1*g, m) for m in (1, 2, 3) if P["lam"][m - 1]), arb(0))
        w = core*moll(de*g)
        tot += -2*w*(acb(0, g*om).exp()/acb(arb("-0.5"), g)).real
    return tot


def job_low(args):
    path, k, P = args
    flint.ctx.dps = DPS
    gs = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i >= k:
                break
            gs.append(arb(line.split("\t")[0].strip()))
    require(len(gs) == k, f"{path}: expected {k} zeros")
    return term_sum(gs, P).str(40), gs[0].str(30), gs[-1].str(30), len(gs), None, None


def job_high(args):
    path, P = args
    flint.ctx.dps = DPS
    z = np.load(path)
    gs = [arb(float(h)) + arb(float(lo)) + arb(0, float(r)) for h, lo, r in zip(z["hi"], z["lo"], z["rad"])]
    return term_sum(gs, P).str(40), gs[0].str(30), gs[-1].str(30), len(gs), None, None


def job_lmfdb(args):
    fname, P = args
    flint.ctx.dps = DPS
    d = platt.decode(os.path.join(DATA, fname))
    idx = d["idx"].astype(np.int64)
    require(np.all(np.diff(idx) == 1), f"{fname}: indices not contiguous")
    two = arb(2)**platt.EPS_EXP; rad = arb(0, arb(2)**platt.RAD_EXP)
    gs, first_idx, next_after, seam = [], None, None, None
    for t0, hi, lo, n in zip(d["t0"], d["hi"], d["lo"], idx):
        if n == N_ARB:
            seam = (arb(float(t0)) + arb((int(hi) << 64) | int(lo))*two + rad).str(40)
        if n <= N_ARB:
            continue
        g = arb(float(t0)) + arb((int(hi) << 64) | int(lo))*two + rad
        if g < arb(P["T0"]):
            if first_idx is None:
                first_idx = int(n)
            gs.append(g)
        else:
            require(g > arb(P["T0"]), f"{fname}: zero ball overlaps T0")
            next_after = (int(n), g.str(30))
            break
    if not gs:
        return arb(0).str(40), None, None, 0, None, next_after, seam
    return term_sum(gs, P).str(40), gs[0].str(30), gs[-1].str(30), len(gs), first_idx, next_after, seam


def main():
    if len(sys.argv) < 8:
        fail(__doc__)
    omega_s, eta1_s, delta_s, M_s, lam_s, T0_s, flist = sys.argv[1:8]
    workers = int(sys.argv[8]) if len(sys.argv) > 8 else 24
    lam = [int(x) for x in lam_s.split(",")]
    require(len(lam) == 3 and sum(lam) == 100 and min(lam) >= 0, "LAM must be three percentages summing to 100")
    M = int(M_s)
    P = {"omega": omega_s, "eta1": eta1_s, "delta": delta_s, "M": M, "lam": lam, "T0": T0_s}
    flint.ctx.dps = DPS
    om, e1, de, T0 = arb(omega_s), arb(eta1_s), arb(delta_s), arb(T0_s)
    eta0 = e1 + de; u_lo, u_hi = om - eta0, om + eta0
    require(u_lo >= 400, "omega - eta0 >= 400 is required")
    require(e1*T0 > 3*arb.pi(), "eta1 T0 > 3 pi is required by the R4 bounds")
    require(de*T0 < arb.pi()/2, "mollifier series needs delta*T0 < pi/2")
    files = [l.strip() for l in open(flist) if l.strip()]
    t0 = time.time()
    lt = cf.low_tasks(cf.N_LOW)
    ht = cf.high_tasks(N_ARB)
    with Pool(workers) as pool:
        resA = pool.map(job_low, [(p, k, P) for p, k in lt], chunksize=4)
        resA += pool.map(job_high, [(p, P) for p in ht], chunksize=4)
        resB = pool.map(job_lmfdb, [(f, P) for f in files], chunksize=1)
    flint.ctx.dps = DPS
    nA = sum(r[3] for r in resA)
    require(nA == N_ARB, f"our zeros: summed {nA}, expected {N_ARB}")
    for a, b in zip(resA, resA[1:]):
        require(arb(a[2]) < arb(b[1]), "our zero batches not increasing across a boundary")
    resB = [r for r in resB if r[3] > 0 or r[5] is not None or r[6] is not None]
    # contiguity of LMFDB indices: first summed index N_ARB + 1, then consecutive across files
    cur = N_ARB + 1; lastg = arb(resA[-1][2])
    nxt = None
    for r in sorted([r for r in resB if r[3] > 0], key=lambda r: r[4]):
        require(r[4] == cur, f"LMFDB index gap: expected {cur}, got {r[4]}")
        require(lastg < arb(r[1]), "LMFDB zeros not above the previous batch")
        cur += r[3]; lastg = arb(r[2])
    for r in resB:
        if r[5] is not None:
            nxt = r[5] if nxt is None or r[5][0] < nxt[0] else nxt
    require(nxt is not None and nxt[0] == cur, "the first zero above T0 was not found right after the last summed zero")
    N = cur - 1
    # zero counts certified independently of Platt's data (Turing's method in Arb)
    require(exact_int(arb(80000000).zeta_nzeros(), "N(8e7)") == N_ARB, "N(8e7) != N_ARB")
    zz = acb.zeta_zeros(N_ARB, 2)
    require(zz[0].imag < arb(80000000) < zz[1].imag, "8e7 not strictly between gamma_N_ARB and gamma_N_ARB+1")
    require(exact_int(T0.zeta_nzeros(), "N(T0)") == N, f"N(T0) from Arb differs from the number summed ({N})")
    zz = acb.zeta_zeros(N, 2)
    require(zz[0].imag < T0 < zz[1].imag, "T0 not strictly between gamma_N and gamma_N+1")
    # two-sided seam: Platt's zero #N_ARB must coincide with ours (so no zero is dropped or doubled at the seam)
    seams = [r[6] for r in resB if r[6] is not None]
    require(len(seams) == 1, "LMFDB zero #N_ARB not found exactly once")
    require(arb(seams[0]).overlaps(arb(resA[-1][2])), "LMFDB zero #N_ARB does not match our zero #N_ARB")
    S = sum((arb(r[0]) for r in resA), arb(0)) + sum((arb(r[0]) for r in resB), arb(0))
    l = [arb(x)/100 for x in lam]
    r4_env = bound_R4_envelope(om, e1, T0, *l); r4_st = bound_R4_stieltjes(om, e1, T0, *l)
    r4 = r4_env if r4_env < r4_st else r4_st
    C_M = arb(4**M*math.factorial(M)**2)/arb(2*math.factorial(2*M))
    rate = (M*arb.pi())/(A_PT*de)
    require(rate < 1, "R5 rate must be < 1")
    r5 = (4*C_M*A_PT.log()*(u_hi/2).exp()*rate**(2*M)).upper()
    r_small = ((2*arb.pi()).log()*(-u_lo/2).exp() + (-5*u_lo/2).exp() + arb("1.53")*(-u_lo/6).exp()
               + arb("1.36")*(u_hi/arb(2).log() + 2)*(-u_lo/4).exp()).upper()
    eps = arb("1.5423e-9")
    U = 1 + 3*eps + S + r4 + r5 + r_small
    cert = provenance(P, CODE)
    cert.update({"omega": omega_s, "eta1": eta1_s, "delta": delta_s, "M_moll": M, "lam": lam_s, "T0": T0_s,
                 "N_T0": N, "N_ours": N_ARB, "N_lmfdb": N - N_ARB, "next_zero_above_T0": nxt[1],
                 "S_M": S.str(20, radius=True), "R4_envelope": r4_env.str(10, radius=False),
                 "R4_stieltjes": r4_st.str(10, radius=False), "R5": r5.str(6, radius=False),
                 "R_small": r_small.str(6, radius=False), "U": U.str(15, radius=True),
                 "U_upper": U.upper().str(10, radius=False),
                 "u_lo": u_lo.lower().str(16, radius=False), "u_hi": u_hi.upper().str(16, radius=False),
                 "log10_x_upper": (u_hi/arb(10).log()).upper().str(12, radius=False),
                 "elapsed_s": round(time.time() - t0, 1)})
    print(json.dumps(cert, indent=2))
    json.dump(cert, open(os.path.join(HERE, f"first_cert_lmfdb_T{T0_s}.json"), "w"), indent=2)
    require(U < 0, "U is not certified negative")
    print(f"CERTIFIED: Delta_1(x) < 0 for some x in [exp({cert['u_lo']}), exp({cert['u_hi']})], so X < 10^{cert['log10_x_upper']}")


if __name__ == "__main__":
    try:
        main()
    except CertificationError as e:
        fail(str(e))
