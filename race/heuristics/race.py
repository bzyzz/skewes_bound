"""Race function E(u) = u e^{-u/2} phi(q) (pi(e^u;q,a) - pi(e^u;q,b)) from zeros of L(s,chi) (heuristic).

E(u) ~ -(N(a)-N(b)) Phi(u) - sum_chi (chi(a)-chi(b)) 2 Re sum_{gamma>0} e^{i gamma u} (1/rho + 1/(rho^2 u)),
N(a) = #{y mod q : y^2 = a}, Phi(u) = 1 + 2/u + 8/u^2 + 48/u^3.  Delta < 0  <=>  E < 0.
"""
import ctypes, os, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PRIM = {  # primitive real characters: conductor -> (q, n) label of saved zeros
    '-3': (3, 2), '-4': (4, 3), '-8': (8, 3), '8': (8, 5), '12': (12, 11), '-24': (24, 5), '24': (24, 11)}
MODCHARS = {3: ['-3'], 4: ['-4'], 8: ['-4', '-8', '8'], 12: ['-3', '-4', '12'],
            24: ['-3', '-4', '-8', '8', '12', '-24', '24']}
_cache = {}


def chival(lab, a):
    from flint import dirichlet_char
    d, n = PRIM[lab]
    z = dirichlet_char(d, n)(a % d)
    return 0 if z == 0 else int(round(float(z.real)))


def zeros(lab, suffix=''):
    key = (lab, suffix)
    if key not in _cache:
        d, n = PRIM[lab]
        _cache[key] = np.load(os.path.join(HERE, 'zeros', f'z_{d}_{n}{suffix}.npy'))
    return _cache[key]


def nsqrt(q, a):
    return sum(1 for y in range(q) if (y*y - a) % q == 0)


def weight(g, kind, T):
    x = g/T
    if kind == 'full':
        return (x <= 1).astype(float)
    if kind == 'jp':
        x = np.minimum(x, 1.0)
        return (1-x)*np.cos(np.pi*x) + np.sin(np.pi*x)/np.pi
    if kind == 'fejer':
        return np.maximum(1-x, 0)
    raise ValueError(kind)


def coeffs(q, a, b, kind='full', T=2e4, suffix=''):
    gs, A, B = [], [], []
    for lab in MODCHARS[q]:
        c = -(chival(lab, a) - chival(lab, b))
        if c == 0:
            continue
        g = zeros(lab, suffix)
        g = g[g <= T]
        w = weight(g, kind, T)
        rho = 0.5 + 1j*g
        gs.append(g); A.append(c*2*w/rho); B.append(c*2*w/rho**2)
    g = np.concatenate(gs); A = np.concatenate(A); B = np.concatenate(B)
    c0 = -(nsqrt(q, a) - nsqrt(q, b))
    return c0, g, A, B


_lib = ctypes.CDLL(os.path.join(HERE, 'libscan.so'))
_dp = ctypes.POINTER(ctypes.c_double)
_lib.scan.argtypes = [_dp]*5 + [ctypes.c_long, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                                ctypes.c_long, ctypes.c_long, _dp, _dp, _dp]


def scan(c0, g, A, B, u0, u1, h, blk=2048, full=False):
    npts = int(math.ceil((u1-u0)/h)) + 1
    nblk = (npts + blk - 1)//blk
    arrs = [np.ascontiguousarray(x, dtype=np.float64) for x in (g, A.real, A.imag, B.real, B.imag)]
    bmin = np.empty(nblk); bu = np.empty(nblk)
    out = np.empty(npts) if full else None
    p = lambda x: x.ctypes.data_as(_dp)
    _lib.scan(*[p(x) for x in arrs], len(g), float(c0), float(u0), float(h), npts, blk,
              p(bmin), p(bu), p(out) if full else None)
    return (bmin, bu, out) if full else (bmin, bu)
