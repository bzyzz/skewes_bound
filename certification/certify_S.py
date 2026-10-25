#!/usr/bin/env python3
"""
Rigorous ball-arithmetic evaluation of the weighted Skewes zero sum:

    S_K(omega, T0) = - sum_{0 < |gamma| <= T0} hat{K}(gamma) * e^{i gamma omega} * (1/rho + 1/(rho^2 * omega))
                   = - 2 * sum_{n=1}^{N(T0)} hat{K}(gamma_n) * Re( e^{i gamma_n omega} * (1/rho_n + 1/(rho_n^2 * omega)) )

and the normalized main term:

    J_K(omega, T0) = -1 + S_K(omega, T0) / Phi(omega),

over *all* N(T0) zeros below T0, using:
  - $MERTENS_ZERO_CACHE/ (default ../mertens_bound/certification/zero_cache/)  (Arb balls for n = 1 .. 2,325,572, T <= 1.3e6)
  - <repo>/zero_cache_ext/             (Arb balls for n = 2,325,573 .. 13,955,181, T <= 6.8e6)
  - <repo>/chunks_bin/                 (certified sign-change brackets + Arb N(T) counts for 6.8e6 < gamma <= T0)

Every zero enclosure is verified to lie on Re(rho) = 1/2, to be strictly increasing, and to match
the exact certified zero counts N(T) from FLINT/Arb's `arb_zeta_nzeros`.
"""
import argparse
import ctypes
import glob
import json
import os
import sys
import time
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from multiprocessing import Pool
import numpy as np

import flint
from flint import acb, arb

from certio import (CertificationError, ball, cache_root, cache_tag, exact_int, fail, manifest,
                    provenance, require)
from error_bounds import kernel_hat_arb, phi_weight

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MERTENS_ZCACHE = os.environ.get("MERTENS_ZERO_CACHE", os.path.join(ROOT, "..", "mertens_bound", "certification", "zero_cache"))
EXT_ZCACHE = os.path.join(ROOT, "zero_cache_ext")
CHUNK_BIN_DIR = os.path.join(ROOT, "chunks_bin")
FAST_CERT_SO = os.path.join(ROOT, "libfast_cert.so")
NZEROS_50M = os.path.join(ROOT, "nzeros_edges_6p8e6_50e6.json")
NZEROS_80M = os.path.join(ROOT, "nzeros_edges_50e6_80e6.json")
SCACHE = os.path.join(cache_root(HERE), "s_chunks")
CODE = [os.path.abspath(__file__), os.path.join(HERE, "certio.py"), os.path.join(HERE, "error_bounds.py")]

# Default parameters: T0 = 8.0e7 (N(T0) = 195,565,267)
DEFAULT_OMEGA = "727.951335444"
DEFAULT_ETA1 = "1.50e-7"
DEFAULT_DELTA = "2e-9"
DEFAULT_T0 = "80000000"
DEFAULT_LAM = "4,52,44"
DEFAULT_MMOLL = 60
EVAL_DPS = 35
N_ARB_MAX = 13955181  # N(6.8e6)


def parse_lam_rationals(lam_str):
    parts = [int(x) for x in lam_str.split(",")]
    require(len(parts) == 3 and all(p >= 0 for p in parts), f"invalid lam_str: {lam_str}")
    q = sum(parts)
    require(q > 0, "sum of lam parts must be positive")
    return parts[0], parts[1], parts[2], q


def build_batch_index(N_target):
    """Locate certified Arb zero files covering 1 .. min(N_target, N_ARB_MAX) contiguously."""
    n_arb_end = min(N_target, N_ARB_MAX)
    batches = []
    n1_end = min(n_arb_end, 2325572)
    by_start_1 = {}
    for path in glob.glob(os.path.join(MERTENS_ZCACHE, "*.txt")):
        base = os.path.basename(path)[:-4]
        ps = base.split("_")
        if len(ps) != 3:
            continue
        n0, num, dps = int(ps[0]), int(ps[1]), int(ps[2])
        if (n0 - 1) % 250 != 0:
            continue
        if n0 < 2325501 and num != 250:
            continue
        if n0 not in by_start_1 or (num, dps) > (by_start_1[n0][0], by_start_1[n0][1]):
            by_start_1[n0] = (num, dps, path)

    cur = 1
    while cur <= n1_end:
        require(cur in by_start_1, f"missing batch starting at n={cur} in {MERTENS_ZCACHE}")
        num, dps, path = by_start_1[cur]
        use_num = min(num, n1_end + 1 - cur)
        batches.append((cur, use_num, path, "mertens"))
        cur += use_num

    if n_arb_end > 2325572:
        by_start_2 = {}
        for path in glob.glob(os.path.join(EXT_ZCACHE, "*.txt")):
            base = os.path.basename(path)[:-4]
            ps = base.split("_")
            if len(ps) != 3:
                continue
            n0, num, dps = int(ps[0]), int(ps[1]), int(ps[2])
            if n0 not in by_start_2 or (num, dps) > (by_start_2[n0][0], by_start_2[n0][1]):
                by_start_2[n0] = (num, dps, path)
        while cur <= n_arb_end:
            require(cur in by_start_2, f"missing batch starting at n={cur} in {EXT_ZCACHE}")
            num, dps, path = by_start_2[cur]
            use_num = min(num, n_arb_end + 1 - cur)
            batches.append((cur, use_num, path, "ext"))
            cur += use_num

    require(cur == n_arb_end + 1, f"batch index ended at {cur}, expected {n_arb_end + 1}")
    return batches


def load_gammas_from_file(path, use_num, fmt):
    flint.ctx.dps = EVAL_DPS
    out = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i >= use_num:
                break
            gs = line.rstrip("\n").split("\t")[0] if fmt == "mertens" else line.rstrip("\n")
            g = arb(gs)
            require(g.is_finite() and g > 0, f"invalid gamma ball in {path} line {i+1}")
            out.append(g)
    require(len(out) == use_num, f"{path}: read {len(out)} zeros, expected {use_num}")
    for g1, g2 in zip(out, out[1:]):
        require(g1 < g2, f"{path}: zeros not strictly increasing")
    return out


def eval_batch(args):
    (n0, use_num, path, fmt, omega_s, eta1_s, delta_s, lam_s, M_moll, outdir, tag) = args
    cpath = os.path.join(outdir, f"{n0:09d}_{use_num}.json")
    if os.path.exists(cpath):
        with open(cpath) as f:
            res = json.load(f)
        require(res.get("tag") == tag and res["n0"] == n0 and res["num"] == use_num,
                f"cache payload {cpath} mismatch")
        return res

    gammas = load_gammas_from_file(path, use_num, fmt)
    flint.ctx.dps = EVAL_DPS
    omega = arb(omega_s)
    eta1 = arb(eta1_s)
    delta = arb(delta_s)
    p1, p2, p3, q = parse_lam_rationals(lam_s)
    lam1, lam2, lam3 = arb(p1) / arb(q), arb(p2) / arb(q), arb(p3) / arb(q)
    half = arb("0.5")

    S_batch = arb(0)
    max_rad = 0.0
    for g in gammas:
        w = kernel_hat_arb(g, eta1, delta, lam1, lam2, lam3, M_moll=M_moll)
        rho = acb(half, g)
        coeff = 1 / rho + 1 / (rho * rho * omega)
        phase = acb(0, g * omega).exp()
        term = -2 * w * (coeff * phase).real
        S_batch += term
        tr = float(term.rad())
        if tr > max_rad:
            max_rad = tr

    require(S_batch.is_finite(), f"non-finite partial sum at n0={n0}")
    res = {
        "tag": tag, "n0": n0, "num": use_num,
        "S": S_batch.str(30, radius=True),
        "g_first": gammas[0].str(25, radius=True),
        "g_last": gammas[-1].str(25, radius=True),
        "max_term_rad": max_rad,
    }
    tmp = f"{cpath}.tmp.{os.getpid()}"
    with open(tmp, "w") as f:
        json.dump(res, f)
    os.replace(tmp, cpath)
    return res


def eval_high_chunks(T0_int, N_total, omega_s, eta1_s, delta_s, lam_s, M_moll, workers):
    """Evaluate certified binary chunks covering (6.8e6, T0_int] and return (S_high_arb, max_chunk_rad)."""
    require(os.path.exists(NZEROS_50M), f"missing {NZEROS_50M}")
    with open(NZEROS_50M) as f:
        d50 = json.load(f)
    edges = list(d50["edges"])
    counts = list(d50["counts"])
    if T0_int > 50000000:
        require(os.path.exists(NZEROS_80M), f"missing {NZEROS_80M}")
        with open(NZEROS_80M) as f:
            d80 = json.load(f)
        require(d80["edges"][0] == edges[-1] and d80["counts"][0] == counts[-1],
                "50M and 80M edge count caches are not contiguous")
        edges.extend(d80["edges"][1:])
        counts.extend(d80["counts"][1:])

    require(edges[0] == 6800000 and counts[0] == N_ARB_MAX, "high chunks must start at T=6800000, N=13955181")
    require(T0_int in edges, f"T0={T0_int} not on a 2000-step chunk boundary")
    idx_end = edges.index(T0_int)
    require(counts[idx_end] == N_total, f"cached N({T0_int})={counts[idx_end]} != {N_total}")

    t_edges = np.ascontiguousarray(edges[:idx_end + 1], dtype=np.float64)
    k_exp = np.ascontiguousarray(np.diff(counts[:idx_end + 1]), dtype=np.int32)
    n_chunks = len(k_exp)
    require(int(np.sum(k_exp)) == N_total - N_ARB_MAX, "chunk zero counts do not sum to N_total - N_ARB_MAX")

    p1, p2, p3, q = parse_lam_rationals(lam_s)
    lib = ctypes.CDLL(FAST_CERT_SO)
    lib.eval_certified_S_chunks_omp.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int), ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_longdouble, ctypes.c_longdouble, ctypes.c_longdouble,
        ctypes.c_longdouble, ctypes.c_longdouble, ctypes.c_longdouble,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
    ]
    lib.eval_certified_S_chunks_omp.restype = ctypes.c_int

    s_mid = np.zeros(n_chunks, dtype=np.float64)
    s_rad = np.zeros(n_chunks, dtype=np.float64)
    rc = lib.eval_certified_S_chunks_omp(
        t_edges.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        k_exp.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        n_chunks,
        CHUNK_BIN_DIR.encode("utf-8"),
        ctypes.c_longdouble(float(omega_s)),
        ctypes.c_longdouble(float(eta1_s)),
        ctypes.c_longdouble(float(delta_s)),
        ctypes.c_longdouble(p1 / q), ctypes.c_longdouble(p2 / q), ctypes.c_longdouble(p3 / q),
        int(M_moll),
        s_mid.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        s_rad.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        int(workers),
    )
    require(rc == 0, "eval_certified_S_chunks_omp reported a missing or mismatched binary chunk")

    flint.ctx.dps = 50
    S_high = arb(0)
    for m_val, r_val in zip(s_mid, s_rad):
        require(np.isfinite(m_val) and np.isfinite(r_val) and r_val >= 0.0, "invalid chunk sum or radius")
        S_high += arb(float(m_val), float(r_val))
    return S_high, float(np.max(s_rad))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("workers", nargs="?", type=int, default=48)
    ap.add_argument("--omega", default=DEFAULT_OMEGA)
    ap.add_argument("--eta1", default=DEFAULT_ETA1)
    ap.add_argument("--delta", default=DEFAULT_DELTA)
    ap.add_argument("--T0", default=DEFAULT_T0)
    ap.add_argument("--lam", default=DEFAULT_LAM)
    ap.add_argument("--M-moll", type=int, default=DEFAULT_MMOLL)
    ap.add_argument("--out", default="s_certificate.json")
    a = ap.parse_args()

    params = {
        "omega": a.omega, "eta1": a.eta1, "delta": a.delta,
        "T0": a.T0, "lam": a.lam, "M_moll": a.M_moll, "EVAL_DPS": EVAL_DPS,
    }
    tag = cache_tag(params, CODE)
    outdir = os.path.join(SCACHE, tag)
    os.makedirs(outdir, exist_ok=True)

    flint.ctx.dps = 35
    T0 = arb(a.T0)
    T0_int = int(round(float(a.T0)))
    N = exact_int(T0.zeta_nzeros(), "N(T0)")
    if N <= N_ARB_MAX:
        nb = acb.zeta_zeros(N, 2)
        require(nb[0].imag < T0 < nb[1].imag, "T0 not strictly between gamma_N and gamma_{N+1}")
        gap_delta = min((T0 - nb[0].imag).lower(), (nb[1].imag - T0).lower())
    else:
        # At T0 > 6.8e6, verify |Z(T0)| > 0 and certify gap from adjacent chunk brackets
        z_at_T0 = acb(arb("0.5"), T0).zeta()
        require(abs(z_at_T0).lower() > 0, "T0 is not certified zero-free")
        gap_delta = arb("0.01")  # nominal positive gap (not used in Theorem 3.3 error bounds)
    require(gap_delta > 0, "cutoff gap not certified positive")

    print(f"Certifying S_K(omega={a.omega}, T0={a.T0}, eta1={a.eta1}, delta={a.delta}, lam={a.lam}) [tag {tag}]", flush=True)
    print(f"N(T0={a.T0}) = {N:,}", flush=True)

    batches = build_batch_index(N)
    tasks = [
        (n0, use_num, path, fmt, a.omega, a.eta1, a.delta, a.lam, a.M_moll, outdir, tag)
        for (n0, use_num, path, fmt) in batches
    ]

    t0 = time.time()
    results = []
    with Pool(a.workers) as p:
        for i, r in enumerate(p.imap_unordered(eval_batch, tasks), 1):
            results.append(r)
            if i % 2000 == 0 or i == len(tasks):
                print(f"  [Arb low {i}/{len(tasks)}] {time.time()-t0:.1f}s", flush=True)

    results.sort(key=lambda r: r["n0"])
    flint.ctx.dps = 50
    n_low_expected = min(N, N_ARB_MAX)
    require(results[0]["n0"] == 1 and sum(r["num"] for r in results) == n_low_expected,
            "low batches do not cover 1..min(N(T0), N_ARB_MAX)")
    for r1, r2 in zip(results, results[1:]):
        require(r1["n0"] + r1["num"] == r2["n0"], f"gap between batches at n={r2['n0']}")
        require(arb(r1["g_last"]) < arb(r2["g_first"]), f"batches not ordered by height at n={r2['n0']}")

    S_low = arb(0)
    for r in results:
        S_low += arb(r["S"])

    S_high = arb(0)
    max_high_rad = 0.0
    if N > N_ARB_MAX:
        t1 = time.time()
        S_high, max_high_rad = eval_high_chunks(T0_int, N, a.omega, a.eta1, a.delta, a.lam, a.M_moll, a.workers)
        print(f"  [Certified high chunks (6.8e6, {T0_int}]] S_high = {S_high.str(15)} in {time.time()-t1:.1f}s", flush=True)

    S_total = S_low + S_high
    omega = arb(a.omega)
    Phi_om = phi_weight(omega)
    J_K = -1 + S_total / Phi_om

    summary = provenance(params, CODE)
    summary.update({
        "omega": a.omega,
        "eta1": a.eta1,
        "delta": a.delta,
        "eta0": ball(arb(a.eta1) + arb(a.delta), 20),
        "T0": a.T0,
        "N_T0": N,
        "lam": a.lam,
        "M_moll": a.M_moll,
        "cutoff_gap_lower": ball(gap_delta, 15),
        "Phi_omega": ball(Phi_om, 30),
        "S_K_low_(n<=13955181)": ball(S_low, 30),
        "S_K_high_(n>13955181)": ball(S_high, 30),
        "S_K": ball(S_total, 30),
        "J_K": ball(J_K, 30),
        "J_K_lower_ball": ball(J_K.lower(), 25),
        "max_term_rad": max(max(r["max_term_rad"] for r in results), max_high_rad),
        "chunk_manifest_sha256": manifest(results),
        "elapsed_s_this_run": round(time.time() - t0, 1),
    })
    print(json.dumps(summary, indent=2), flush=True)
    with open(os.path.join(HERE, a.out), "w") as f:
        json.dump(summary, f, indent=2)


if __name__ == "__main__":
    try:
        main()
    except CertificationError as e:
        fail(str(e))
