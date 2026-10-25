#!/usr/bin/env python3
"""Combine the Arb batch sums of arb_recert.py into Arb-only certificates <name>_arb.json.

Checks:
- the batches cover n = 1 .. N(8e7) contiguously, with increasing, disjoint first/last zeros;
- arb_zeta_nzeros certifies N(T0), and gamma_N < T0 < gamma_{N+1}, for T0 = 5e7 and 8e7.

Each output has the same schema as <name>.json, with S_K, J_K and J_K_lower_ball replaced by the Arb-only
values. Check it with `python check_theorem.py <name>_arb.json`.
"""
import glob, hashlib, json, os
import flint
from flint import arb, acb
from certio import require

HERE = os.path.dirname(os.path.abspath(__file__))
BATCHES = os.path.join(HERE, "arb_batches")
NAMES = ["s_cert_T5e7", "s_certificate", "s_cert_T8e7_st"]


def main():
    flint.ctx.dps = 40
    paths = sorted(glob.glob(os.path.join(BATCHES, "*.json")))
    rs = sorted((json.load(open(p)) for p in paths), key=lambda r: r["n0"])
    require(rs[0]["n0"] == 1, "batches must start at n = 1")
    for a, b in zip(rs, rs[1:]):
        require(a["n0"] + a["k"] == b["n0"], f"gap after n = {a['n0']}")
        require(arb(a["last"]) < arb(b["first"]), f"zeros not increasing at n = {b['n0']}")
    N = rs[-1]["n0"] + rs[-1]["k"] - 1
    print(f"zeros covered: 1 .. {N}; max ball radius {max(r['maxrad'] for r in rs):.3g}")
    flint.ctx.dps = 25
    cut = {}
    for T0, n in ((50000000, 118488123), (80000000, 195565267)):
        require(arb(T0).zeta_nzeros().contains(n), f"N({T0}) != {n}")
        z = acb.zeta_zeros(n, 2)
        require(z[0].imag < arb(T0) < z[1].imag, f"cutoff check at {T0}")
        cut[T0] = n
        print(f"N({T0}) = {n} certified; gamma_N < T0 < gamma_N+1")
    manifest = hashlib.sha256("".join(open(p).read() for p in paths).encode()).hexdigest()
    flint.ctx.dps = 40
    for name in NAMES:
        d = json.load(open(os.path.join(HERE, f"{name}.json")))
        require(d["N_T0"] in cut.values() and N >= d["N_T0"], f"{name}: zeros do not reach N(T0)")
        S = arb(0)
        for r in rs:
            S += arb(r["S"][name])          # each batch sums only zeros with n <= N(T0) of this certificate
        om = arb(d["omega"])
        Phi = 1 + 2/om + 8/om**2 + arb("58.56")/om**3
        J = -1 + S/Phi
        out = dict(d)
        for k in [k for k in out if k.startswith("S_K_")]:
            del out[k]
        out["code_sha256"] = {f: hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest()
                              for f in ("arb_recert.py", "arb_combine.py", "certio.py", "error_bounds.py")}
        out.update(S_K=S.str(15), J_K=J.str(15), J_K_lower_ball=J.lower().str(25, radius=False),
                   zero_source="Arb only: zero caches (n <= 13,955,181) and acb.zeta_zeros re-isolation (arb_recert.py)",
                   batch_manifest_sha256=manifest, n_batches=len(rs))
        json.dump(out, open(os.path.join(HERE, f"{name}_arb.json"), "w"), indent=1)
        print(f"{name}_arb.json: S_K = {S.str(12)}, J_K >= {J.lower().str(10, radius=False)} (float-zero certificate: {d['J_K_lower_ball']})")


if __name__ == "__main__":
    main()
