# First integer x with Delta_K(x) < 0, last such x below X, and log-density of negativity, for Q(sqrt 5), Q(sqrt 13).
import math, sys
import numpy as np
from flint import acb, dirichlet_char
X = int(float(sys.argv[1])); EG = 0.5772156649015329
sieve = np.ones(X + 1, dtype=bool); sieve[:2] = False
for i in range(2, int(X**0.5) + 1):
    if sieve[i]: sieve[i*i::i] = False
P = np.nonzero(sieve)[0].astype(np.int64)
for name, q, n in [("Q(sqrt 5)", 5, 4), ("Q(sqrt 13)", 13, 12)]:
    c = dirichlet_char(q, n); kappa = float(acb.dirichlet_l(acb(1), c).real.mid())
    chi = np.array([int(round(float(c(a).real))) for a in range(q)])[P % q]
    split, inert, ram = P[chi == 1], P[chi == -1], P[chi == 0]
    inert = inert[inert.astype(np.float64)**2 <= X]
    norms = np.concatenate([split, inert**2, ram]).astype(np.int64)
    mult = np.concatenate([2*np.ones(len(split)), np.ones(len(inert)), np.ones(len(ram))])
    o = np.argsort(norms, kind="stable"); norms, mult = norms[o], mult[o]
    S = np.cumsum(-mult*np.log1p(-1/norms.astype(np.float64)))
    # integer x in [norms[k], norms[k+1]-1]: E3 decreasing, so the worst integer is norms[k+1]-1
    nxt = np.append(norms[1:], X + 1) - 1
    E3 = S - EG - math.log(kappa) - np.log(np.log(nxt.astype(np.float64)))
    neg = np.nonzero(E3 < 0)[0]
    # first integer x with E3 < 0: within [norms[k], nxt[k]], solve log log x > S_k - EG - log kappa
    k = neg[0]; xs = math.ceil(math.exp(math.exp(S[k] - EG - math.log(kappa))))
    xs = max(xs, int(norms[k]))
    # log-density of {x <= X: Delta_K(x) < 0}: measure in log x of the negative parts
    thr = np.exp(np.exp(S - EG - math.log(kappa)))             # E3 < 0 for x > thr within [norms[k], next norm)
    lo = np.maximum(thr, norms.astype(np.float64)); hi = np.append(norms[1:], X).astype(np.float64)
    meas = np.sum(np.where(hi > lo, np.log(hi) - np.log(np.maximum(lo, 1)), 0.0))
    dens = meas/ (math.log(X) - math.log(2))
    print(f"{name}: first integer x with Delta_K(x) < 0: x = {xs}; last negative interval below {X:.0e} ends at {nxt[neg[-1]]}; "
          f"log-density of negativity on [2, {X:.0e}] = {dens:.4f}", flush=True)
