#!/usr/bin/env python3
"""Random-phase (LI) prediction of the first sign change of the team race modulo q.

Model (as in mertens_bound/search/ld_model.py): -E(u) - (r-1) = h(u) = sum_j w_j cos(gamma_j u + phi_j) with
w_j = 2/|rho_j| over the zeros (gamma_j <= T) of the real characters mod q; residues lead when h(u) >= H := (r-1).
Under LI the phases are independent and uniform, the saddle point gives P(h >= H), and Rice's formula the
rate of upcrossings per unit u. The predicted first crossing is u0 + 1/rate (rate is independent of u here).
usage: ld_team.py q [T]
"""
import json, math, os, sys
import numpy as np
from scipy.optimize import brentq
from scipy.special import i0e, i1e
import glob
HERE = os.path.dirname(os.path.abspath(__file__))
_RP = json.load(open(os.path.join(HERE, "real_prims.json")))   # written by team_scan.real_prims (needs flint)


def real_prims(q):
    return [tuple(x) for x in _RP[str(q)]]


def load_zeros(d, n, T):
    for p in sorted(glob.glob(os.path.join(HERE, "zeros", f"z_{d}_{n}_T*.npy")), key=lambda p: -float(p.split("_T")[1][:-4])):
        z = np.load(p)
        if z[-1] >= 0.98*T:
            return z[z <= T]
    raise FileNotFoundError(f"{d}.{n} to {T}")


def tilt(w, H):
    L1 = lambda t: np.sum(w*i1e(t*w)/i0e(t*w))
    return brentq(lambda t: L1(t) - H, 1e-8, 1e6)


def rate_and_tail(w, g, H):
    t = tilt(w, H)
    x = t*w
    r = i1e(x)/i0e(x)
    L = np.sum(np.log(i0e(x)) + x)
    L2 = np.sum(w*w*(1 - r/np.where(x > 0, x, 1) - r*r))
    I = t*H - L
    logdens = -I - 0.5*math.log(2*math.pi*L2)
    sd_deriv = math.sqrt(np.sum((w*g)**2*i1e(x)/(np.where(x > 0, x, 1)*i0e(x))))
    return logdens + math.log(sd_deriv/math.sqrt(2*math.pi)), logdens - math.log(t), t


def predict(q, T=1e4):
    g = np.concatenate([load_zeros(d, n, T) for d, n in real_prims(q)])
    w = 2/np.abs(0.5 + 1j*g)
    H = float(len(real_prims(q)))
    lrate, ltail, t = rate_and_tail(w, g, H)
    return dict(q=q, r=len(real_prims(q)) + 1, H=H, T=T, nzeros=len(g), sd=float(np.sqrt(np.sum(w*w)/2)),
                log_rate=lrate, u_first_pred=float(math.exp(-lrate)), log10_x_pred=float(math.exp(-lrate)/math.log(10)))


if __name__ == "__main__":
    q = int(sys.argv[1]); T = float(sys.argv[2]) if len(sys.argv) > 2 else 1e4
    print(json.dumps(predict(q, T)))
