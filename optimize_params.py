#!/usr/bin/env python3
"""
Load all certified zeros up to T0 (e.g. T0 = 4.3e6 or T0 = 6.8e6), build the exact
Taylor-binned evaluation matrix on u in [727.9513350, 727.9513368], and find the
exact parameters (omega, eta1, lam1, lam2, lam3) that minimize the certified upper bound
omega + eta0 with a comfortable positive margin I(omega, eta0) > 1e-5.
"""
import glob
import math
import os
import sys
import time
from multiprocessing import Pool
import numpy as np
import flint
from flint import arb
ROOT = os.path.dirname(os.path.abspath(__file__))

sys.path.insert(0, os.path.join(ROOT, "certification"))
from certio import exact_int
from certify_S import build_batch_index
from error_bounds import bound_R4, total_error, phi_weight

def read_batch_mid(args):
    n0, use_num, path, fmt = args
    out = np.empty(use_num, dtype=np.float64)
    with open(path) as f:
        for i in range(use_num):
            line = f.readline()
            gs = line.split("\t")[0] if fmt == "mertens" else line.strip()
            if gs.startswith("["):
                mid = gs[1:gs.index(" +/-")]
            else:
                mid = gs
            out[i] = float(mid)
    return n0, out

def load_all_gammas(N_target):
    batches = build_batch_index(N_target)
    with Pool(24) as p:
        res = p.map(read_batch_mid, batches)
    res.sort(key=lambda x: x[0])
    gammas = np.concatenate([r[1] for r in res])
    assert len(gammas) == N_target
    assert np.all(np.diff(gammas) > 0)
    return gammas

def kappa_hat_2m_np(xi, m):
    res = np.sinc(xi / np.pi)
    for k in range(1, m + 1):
        kp2 = (k * np.pi) ** 2
        res = res * (kp2 / (kp2 - xi**2))
    return res

def main():
    T0_str = sys.argv[1] if len(sys.argv) > 1 else "4300000"
    flint.ctx.dps = 30
    T0_arb = arb(T0_str)
    T0_val = float(T0_str)
    N_target = exact_int(T0_arb.zeta_nzeros(), "N(T0)")
    print(f"Loading N(T0={T0_str}) = {N_target:,} certified zeros...", flush=True)
    t0 = time.time()
    gammas = load_all_gammas(N_target)
    print(f"Loaded {len(gammas):,} zeros in {time.time()-t0:.2f}s (gamma_max = {gammas[-1]:.4f})", flush=True)

    um = 727.9513358
    rho = 0.5 + 1j * gammas
    z = -2.0 * (1.0 / rho + 1.0 / (rho**2 * um)) * np.exp(1j * gammas * um)

    N_BINS = 2000
    edges = np.linspace(0, N_target, N_BINS + 1, dtype=np.int64)
    bin_gammas = np.zeros(N_BINS, dtype=np.float64)
    Z0 = np.zeros(N_BINS, dtype=np.complex128)
    Z1 = np.zeros(N_BINS, dtype=np.complex128)
    Z2 = np.zeros(N_BINS, dtype=np.complex128)

    for b in range(N_BINS):
        i0, i1 = int(edges[b]), int(edges[b+1])
        gb = 0.5 * (gammas[i0] + gammas[i1 - 1])
        bin_gammas[b] = gb
        dg = gammas[i0:i1] - gb
        zb = z[i0:i1]
        Z0[b] = np.sum(zb)
        Z1[b] = np.sum(1j * dg * zb)
        Z2[b] = np.sum(-0.5 * (dg**2) * zb)

    u_grid = np.linspace(727.95133520, 727.95133660, 1401)  # step 1e-9
    dy = u_grid - um
    phase_mat = np.exp(1j * np.outer(bin_gammas, dy))
    B_mat = np.real((Z0[:, None] + Z1[:, None] * dy[None, :] + Z2[:, None] * (dy[None, :]**2)) * phase_mat)

    thresh = float(phi_weight(arb(um)).mid())
    s_unw = np.sum(B_mat, axis=0)
    idx_unw = np.where(s_unw > thresh)[0][0]
    i_peak = int(np.argmax(s_unw))
    print(f"Full {N_target:,}-zero unweighted profile:")
    print(f"  First u > Phi(u) : u = {u_grid[idx_unw]:.9f} (S = {s_unw[idx_unw]:.8f})")
    print(f"  Peak S_unw       : u = {u_grid[i_peak]:.9f} (S = {s_unw[i_peak]:.8f}, J = {s_unw[i_peak]/thresh - 1:+.8f})", flush=True)

    # Search over rational (p1, p2, p3) with p1 + p2 + p3 = 100 and eta1
    best_ub = 999.0
    best_cand = None
    for p1 in range(8, 20, 2):
        for p2 in range(64, 78, 2):
            p3 = 100 - p1 - p2
            if p3 < 0:
                continue
            l1, l2, l3 = p1 / 100.0, p2 / 100.0, p3 / 100.0
            for xi_T in np.linspace(10.1, 11.2, 23):
                eta1 = round(xi_T / T0_val, 9)
                if eta1 * T0_val <= 3 * math.pi * 1.01:
                    continue
                wb = (l1 * kappa_hat_2m_np(eta1 * bin_gammas, 1)
                      + l2 * kappa_hat_2m_np(eta1 * bin_gammas, 2)
                      + l3 * kappa_hat_2m_np(eta1 * bin_gammas, 3))
                out = wb @ B_mat
                J_vals = out / thresh - 1.0
                if np.max(J_vals) < 2.5e-4:
                    continue
                r_tot = bound_R4(
                    arb(um), arb(str(eta1)), T0_arb,
                    arb(p1)/100, arb(p2)/100, arb(p3)/100, n_sub=512
                )
                r_val = float(r_tot.mid()) + 1e-8
                if r_val > 0.0018:
                    continue
                net = J_vals - r_val
                pos = np.where(net > 1.5e-5)[0]
                if len(pos) > 0:
                    u_first = round(float(u_grid[pos[0]]), 9)
                    ub = u_first + eta1 + 1e-8
                    if ub < best_ub:
                        best_ub = ub
                        best_cand = (p1, p2, p3, eta1, u_first, float(J_vals[pos[0]]), r_val, float(net[pos[0]]))
                        print(f"  NEW BEST: u + eta0 = {best_ub:.9f} | omega={u_first:.9f}, eta1={eta1*1e6:.3f}e-6, lam={p1},{p2},{p3}, J={J_vals[pos[0]]:.6e}, R={r_val:.6e}, margin={net[pos[0]]:.6e}", flush=True)

    print("\n======================================================================")
    if best_cand:
        p1, p2, p3, eta1, u_first, J_v, r_v, marg = best_cand
        print(f"OPTIMAL PARAMETERS FOR T0 = {T0_str} (N = {N_target:,}):")
        print(f"  --T0 {T0_str} --omega {u_first:.9f} --eta1 {eta1:.9f} --delta 1e-8 --lam {p1},{p2},{p3}")
        print(f"  Predicted J_K = {J_v:.7e}, R_tot = {r_v:.7e}, Margin = {marg:.7e}")
        print(f"  Upper bound: omega + eta0 = {best_ub:.9f}")
    print("======================================================================")

if __name__ == "__main__":
    main()
