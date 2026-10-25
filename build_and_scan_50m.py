#!/usr/bin/env python3
"""
Certify all zeta zeros from T = 6.8e6 up to T = 5.0e7 (N(5.0e7) = 118,044,213 zeros!)
in chunks of dT = 2000.0 using:
  1) FLINT/Arb's certified `arb.zeta_nzeros(T_k)` for the exact zero count on every chunk [T_k, T_{k+1}],
  2) Fast C++/OpenMP Taylor-Riemann-Siegel + Gabcke Satz 1 certified sign-change bracketing (`libfast_cert.so`),
  3) Automatic FLINT/Arb `acb.zeta_zeros` fallback on any chunk where a Lehmer pair is closer than 1.5e-5.
Saves certified binary chunks in `<repo>/chunks_bin/` and the exact Taylor-bin
summary in `<repo>/bins_50m.npz`.
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
import ctypes
import json
import math
import struct
import sys
import time
from multiprocessing import Pool
import numpy as np
import flint
from flint import acb, arb
ROOT = os.path.dirname(os.path.abspath(__file__))

sys.path.insert(0, os.path.join(ROOT, "certification"))
from certio import exact_int, require

CHUNK_DIR = os.path.join(ROOT, "chunks_bin")
NZEROS_CACHE = os.path.join(ROOT, "nzeros_edges_6p8e6_50e6.json")
BINS_14M = os.path.join(ROOT, "bins_13955k.npz")
BINS_50M = os.path.join(ROOT, "bins_50m.npz")
U_REF = 727.95133585
DT_CHUNK = 2000.0
T_START = 6800000.0
T_END = 50000000.0


def arb_nzeros_worker(t_val):
    flint.ctx.dps = 25
    t_arb = arb(int(round(t_val)))
    nz = exact_int(t_arb.zeta_nzeros(), f"N({t_val})")
    return (int(round(t_val)), nz)


def get_edge_counts(t_edges):
    if os.path.exists(NZEROS_CACHE):
        with open(NZEROS_CACHE) as f:
            d = json.load(f)
        if len(d["edges"]) == len(t_edges) and d["edges"][0] == int(t_edges[0]) and d["edges"][-1] == int(t_edges[-1]):
            return np.array(d["counts"], dtype=np.int64)
    print(f"Computing exact certified N(T) at {len(t_edges)} chunk boundaries via Arb...", flush=True)
    t0 = time.time()
    with Pool(48) as p:
        res = p.map(arb_nzeros_worker, [float(x) for x in t_edges])
    res.sort(key=lambda x: x[0])
    counts = [x[1] for x in res]
    with open(NZEROS_CACHE, "w") as f:
        json.dump({"edges": [x[0] for x in res], "counts": counts}, f)
    print(f"Computed {len(counts)} exact N(T) endpoints (N_max = {counts[-1]:,}) in {time.time()-t0:.1f}s", flush=True)
    return np.array(counts, dtype=np.int64)


def arb_fallback_chunk(args):
    Ta, Tb, n0, K_exp = args
    flint.ctx.dps = 18
    zs = acb.zeta_zeros(n0, K_exp)
    omids = []
    rads = []
    prev_g = arb(Ta)
    for rho in zs:
        require(rho.real == arb("0.5"), f"Zero near {n0} not on critical line")
        g = rho.imag
        require(prev_g < g, f"Zeros near {n0} not increasing")
        prev_g = g
        off = g - arb(Ta)
        omid = float(off.mid())
        rad = float(off.rad()) + 2e-12
        omids.append(omid)
        rads.append(rad)
    require(prev_g < arb(Tb), f"Last zero in [{Ta}, {Tb}] >= {Tb}")
    fpath = os.path.join(CHUNK_DIR, f"chunk_{int(round(Ta)):08d}_{int(round(Tb)):08d}.bin")
    with open(fpath, "wb") as f:
        f.write(struct.pack("<id", K_exp, float(Ta)))
        for om, rd in zip(omids, rads):
            f.write(struct.pack("<df", om, rd))
    # Compute Taylor moments around U_REF and gb = 0.5*(Ta + Tb)
    gb = 0.5 * (Ta + Tb)
    g_arr = Ta + np.array(omids, dtype=np.float64)
    dg = g_arr - gb
    rho_c = 0.5 + 1j * g_arr
    # Use high-precision phase mod 2pi via flint or mpmath? Wait: in Python we can use arb for the phase!
    z0 = 0j
    z1 = 0j
    z2 = 0j
    u_arb = arb(str(U_REF))
    for i, om in enumerate(omids):
        g_a = arb(Ta) + arb(om)
        ph = float((g_a * u_arb).mid())
        g_f = float(Ta + om)
        dg_f = g_f - gb
        r_c = 0.5 + 1j * g_f
        z_val = -2.0 * (1.0 / r_c + 1.0 / (r_c**2 * U_REF)) * np.exp(1j * ph)
        z0 += z_val
        z1 += 1j * dg_f * z_val
        z2 += -0.5 * (dg_f**2) * z_val
    return (gb, z0, z1, z2, max(rads))


def main():
    os.makedirs(CHUNK_DIR, exist_ok=True)
    t_edges = np.arange(T_START, T_END + 1.0, DT_CHUNK, dtype=np.float64)
    counts = get_edge_counts(t_edges)
    k_exp_all = np.diff(counts).astype(np.int32)
    n_chunks_total = len(k_exp_all)
    print(f"Total chunks in [{T_START:.0f}, {T_END:.0f}]: {n_chunks_total} ({np.sum(k_exp_all):,} zeros, up to N={counts[-1]:,})", flush=True)

    lib = ctypes.CDLL(os.path.join(ROOT, "libfast_cert.so"))
    lib.certify_chunks_omp.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int), ctypes.c_int,
        ctypes.c_double, ctypes.c_char_p,
        ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.c_int,
    ]

    # Load existing 14M bins (T <= 6.8e6)
    d14 = np.load(BINS_14M)
    bg_list = [d14["bin_gammas"]]
    Z0_list = [d14["Z0"]]
    Z1_list = [d14["Z1"]]
    Z2_list = [d14["Z2"]]

    u_grid = np.linspace(727.951335200, 727.951336600, 1401)
    dy = u_grid - U_REF
    thresh = 1.0 + 2.0 / U_REF + 8.0 / (U_REF**2) + 58.56 / (U_REF**3)

    # Process in macro-blocks of 2000 chunks (dT = 4,000,000 per macro-block)
    BLK_CHUNKS = 2000
    t_run = time.time()
    for c0 in range(0, n_chunks_total, BLK_CHUNKS):
        c1 = min(n_chunks_total, c0 + BLK_CHUNKS)
        nb = c1 - c0
        sub_edges = np.ascontiguousarray(t_edges[c0:c1 + 1], dtype=np.float64)
        sub_k = np.ascontiguousarray(k_exp_all[c0:c1], dtype=np.int32)
        status = np.zeros(nb, dtype=np.int32)
        max_rad = np.zeros(nb, dtype=np.float64)
        bg = np.zeros(nb, dtype=np.float64)
        Z0r, Z0i = np.zeros(nb, dtype=np.float64), np.zeros(nb, dtype=np.float64)
        Z1r, Z1i = np.zeros(nb, dtype=np.float64), np.zeros(nb, dtype=np.float64)
        Z2r, Z2i = np.zeros(nb, dtype=np.float64), np.zeros(nb, dtype=np.float64)

        tb0 = time.time()
        lib.certify_chunks_omp(
            sub_edges.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            sub_k.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            nb, U_REF, CHUNK_DIR.encode("utf-8"),
            status.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            max_rad.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            bg.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            Z0r.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), Z0i.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            Z1r.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), Z1i.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            Z2r.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), Z2i.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            64,
        )
        failed_idx = np.where(status == 0)[0]
        if len(failed_idx) > 0:
            print(f"  Running Arb fallback on {len(failed_idx)} chunks with ultra-close Lehmer pairs...", flush=True)
            fb_tasks = [
                (float(sub_edges[idx]), float(sub_edges[idx + 1]), int(counts[c0 + idx]) + 1, int(sub_k[idx]))
                for idx in failed_idx
            ]
            with Pool(min(32, len(fb_tasks))) as p:
                fb_res = p.map(arb_fallback_chunk, fb_tasks)
            for idx, (gb_f, z0_f, z1_f, z2_f, mr_f) in zip(failed_idx, fb_res):
                status[idx] = 1
                max_rad[idx] = mr_f
                bg[idx] = gb_f
                Z0r[idx], Z0i[idx] = z0_f.real, z0_f.imag
                Z1r[idx], Z1i[idx] = z1_f.real, z1_f.imag
                Z2r[idx], Z2i[idx] = z2_f.real, z2_f.imag

        assert np.all(status == 1), "Some chunks failed certification!"
        bg_list.append(bg)
        Z0_list.append(Z0r + 1j * Z0i)
        Z1_list.append(Z1r + 1j * Z1i)
        Z2_list.append(Z2r + 1j * Z2i)

        all_bg = np.concatenate(bg_list)
        all_Z0 = np.concatenate(Z0_list)
        all_Z1 = np.concatenate(Z1_list)
        all_Z2 = np.concatenate(Z2_list)
        np.savez(BINS_50M, bin_gammas=all_bg, Z0=all_Z0, Z1=all_Z1, Z2=all_Z2)

        T_cur = sub_edges[-1]
        N_cur = counts[c1]
        dt_blk = time.time() - tb0
        rate = np.sum(sub_k) / max(dt_blk, 1e-3)
        # Compute current unweighted crossing and values at 727.95133580, 727.95133585, 727.95133590
        B_cur = np.real((all_Z0[:, None] + all_Z1[:, None] * dy[None, :] + all_Z2[:, None] * (dy[None, :]**2))
                        * np.exp(1j * np.outer(all_bg, dy)))
        J_unw = np.sum(B_cur, axis=0) / thresh - 1.0
        pos_u = np.where(J_unw > 0)[0]
        u_cross = u_grid[pos_u[0]] if len(pos_u) > 0 else 999.0
        j_580 = J_unw[int(np.argmin(np.abs(u_grid - 727.951335800)))]
        j_585 = J_unw[int(np.argmin(np.abs(u_grid - 727.951335850)))]
        j_590 = J_unw[int(np.argmin(np.abs(u_grid - 727.951335900)))]
        print(f"[T <= {T_cur/1e6:4.1f}e6 | N={N_cur:11,d}] {dt_blk:5.1f}s ({rate:,.0f} z/s, max_rad={np.max(max_rad):.2e}) "
              f"| u_cross={u_cross:.9f} | J(.580)={j_580:+.6f}, J(.585)={j_585:+.6f}, J(.590)={j_590:+.6f}", flush=True)

    print(f"\nALL {counts[-1]:,} ZEROS CERTIFIED UP TO T = {T_END:.0f} in {time.time()-t_run:.1f}s!", flush=True)


if __name__ == "__main__":
    main()
