#!/usr/bin/env python3
"""
Complete, unaliased scan of u in [10, 728] at resolution du = 9.54e-7
using all 2,325,572 certified zeta zeros via Type-1 NUFFT (6-pt Lagrange + FFT),
parallelized across 18 worker processes, followed by exact polishing on top candidates.
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import time
import ctypes
from multiprocessing import Pool
import numpy as np
ROOT = os.path.dirname(os.path.abspath(__file__))

gammas = np.ascontiguousarray(np.load(os.path.join(ROOT, "gammas_2325572.npy")), dtype=np.float64)
N_total = len(gammas)
denom = 0.25 + gammas**2

L = 16.0
N_fft = 1 << 24  # 16,777,216 points -> du = 9.536743e-7
du = L / N_fft
N_freq = N_fft // 2 + 1
M_half = int(round(4.0 / du))  # 4,194,304 points -> active window [-4.0, 4.0] of width 8.0

z = gammas * (L / (2.0 * np.pi))
m0 = np.floor(z).astype(np.int64)
delta = z - m0

nodes = np.array([-2.0, -1.0, 0.0, 1.0, 2.0, 3.0])
weights = np.zeros((6, N_total), dtype=np.float64)
for idx_s, s in enumerate(nodes):
    w = np.ones(N_total, dtype=np.float64)
    for idx_t, t in enumerate(nodes):
        if idx_s != idx_t:
            w *= (delta - t) / (s - t)
    weights[idx_s] = w

bin_idx = np.zeros((6, N_total), dtype=np.int64)
for i, s in enumerate(nodes):
    bin_idx[i] = m0 + int(s)

def scan_window(uc):
    # Note: np.fft.irfft(X, n=N)[k] = (1/N) * (X[0] + 2 * Re(sum_{m=1}^{N/2-1} X[m] e^{i 2pi m k / N}))
    # So if X[m] carries (-1/rho - 1/(rho^2 * uc)) * e^{i gamma uc}, then N * irfft(X) gives
    # 2 * Re(sum (-1/rho - 1/(rho^2 * uc)) e^{i gamma u}) = S(u)!
    c_uc = (-1.0 / (0.5 + 1j * gammas) - 1.0 / ((0.5 + 1j * gammas)**2 * uc)) * np.exp(1j * gammas * uc)
    spectrum = np.zeros(N_freq, dtype=np.complex128)
    for i in range(6):
        np.add.at(spectrum, bin_idx[i], c_uc * weights[i])
    vals = np.fft.irfft(spectrum, n=N_fft) * N_fft
    act = np.concatenate([vals[N_fft - M_half:], vals[:M_half + 1]])
    u_min_w = uc - 4.0
    n_blk = len(act) // 4096
    blk_max = act[:n_blk * 4096].reshape(n_blk, 4096)
    max_per_blk = np.max(blk_max, axis=1)
    # Record any block with max > 0.85, plus the overall window max
    hot_blks = np.where(max_per_blk > 0.85)[0]
    cands = []
    for b in hot_blks:
        sub_idx = int(np.argmax(blk_max[b]))
        k_act = int(b * 4096 + sub_idx)
        u_est = u_min_w + k_act * du
        if 10.0 <= u_est <= 730.0:
            cands.append((u_est, float(act[k_act])))
    k_best = int(np.argmax(act))
    return uc, float(act[k_best]), u_min_w + k_best * du, cands

if __name__ == "__main__":
    t0 = time.time()
    centers = list(np.arange(14.0, 730.0, 8.0))
    print(f"Scanning {len(centers)} windows covering u in [10, 730] with 18 workers...", flush=True)
    all_cands = []
    with Pool(18) as p:
        for i, (uc, wmax, ubest, cands) in enumerate(p.imap_unordered(scan_window, centers), 1):
            all_cands.extend(cands)
            if wmax > 0.85 or i % 15 == 0 or i == len(centers):
                print(f"  [{i:2d}/{len(centers)}] uc={uc:5.1f} [{uc-4:.0f},{uc+4:.0f}] max={wmax:.6f} at u={ubest:.6f} (cands={len(cands)}) [{time.time()-t0:.1f}s]", flush=True)

    print(f"\nFFT scan finished in {time.time()-t0:.1f}s. Total candidates > 0.85: {len(all_cands)}", flush=True)
    # Polish all candidates using exact OpenMP evaluator
    lib = ctypes.CDLL(os.path.join(ROOT, "libfast_skewes.so"))
    lib.eval_skewes_uniform.argtypes = [
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
        ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.POINTER(ctypes.c_double),
    ]
    c1_cos = np.ascontiguousarray(-1.0 / denom, dtype=np.float64)
    c1_sin = np.ascontiguousarray(-2.0 * gammas / denom, dtype=np.float64)
    c2_cos = np.ascontiguousarray(-(0.5 - 2.0 * gammas**2) / (denom**2), dtype=np.float64)
    c2_sin = np.ascontiguousarray(-2.0 * gammas / (denom**2), dtype=np.float64)

    polished = []
    for u_est, s_est in all_cands:
        u0 = u_est - 5e-6
        u1 = u_est + 5e-6
        M = 1001
        du_p = (u1 - u0) / (M - 1)
        cc = np.ascontiguousarray(c1_cos + c2_cos / u_est, dtype=np.float64)
        cs = np.ascontiguousarray(c1_sin + c2_sin / u_est, dtype=np.float64)
        out = np.zeros(M, dtype=np.float64)
        lib.eval_skewes_uniform(
            gammas.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            cc.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            cs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_int(N_total),
            ctypes.c_double(u0), ctypes.c_double(du_p), ctypes.c_int(M),
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        )
        idx = int(np.argmax(out))
        u_pol = u0 + idx * du_p
        s_pol = float(out[idx])
        thresh = 1.0 + 2.0 / u_pol + 8.0 / (u_pol**2) + 48.0 / (u_pol**3)
        polished.append((u_pol, s_pol, thresh, s_pol - thresh))

    polished.sort(key=lambda x: -x[3])
    print("\nTop 30 peaks on u in [10, 730] sorted by (S - threshold):", flush=True)
    for u_pol, s_pol, thresh, diff in polished[:30]:
        log10_x = u_pol / np.log(10.0)
        print(f"  u = {u_pol:13.8f} (x ~ 10^{log10_x:8.4f}) | S = {s_pol:.8f} | thresh = {thresh:.8f} | S - thresh = {diff:+.8f}", flush=True)
