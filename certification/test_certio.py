#!/usr/bin/env python3
"""
Tests for certio.py and error_bounds.py failure modes and outward-rounding guarantees.
"""
import os
import random
import subprocess
import sys
from fractions import Fraction
import flint
from flint import arb

from certio import CertificationError, ball, directed_decimal, exact_int, require
from error_bounds import bound_R4_envelope, bound_R4_stieltjes, kappa_2m_hat_arb, total_error

HERE = os.path.dirname(os.path.abspath(__file__))

flint.ctx.dps = 40


def must_fail(fn):
    try:
        fn()
    except CertificationError:
        return
    raise AssertionError(f"{fn} did not raise CertificationError")


def main():
    # 1. Outward ball serialisation encloses original ball
    x = arb.pi() / 7
    b = arb(ball(x, 20))
    assert b.contains(x), "serialised ball does not contain original"

    # 2. exact_int accepts unique integer and rejects wide/non-integer balls
    assert exact_int(arb(42), "exact") == 42
    assert exact_int(arb(42, arb("0.1")), "small rad") == 42
    must_fail(lambda: exact_int(arb("42.5", arb("0.6")), "straddling"))
    must_fail(lambda: exact_int(arb("42.5"), "non-int"))

    # 3. Removable singularity handling at xi = k*pi in kappa_2m_hat_arb
    pi = arb.pi()
    for m in [1, 2, 3]:
        # At xi = 0, hat{kappa}_{2m}(0) = 1
        v0 = kappa_2m_hat_arb(arb(0), m)
        assert v0.contains(1), f"hat{{kappa}}_{2*m}(0) does not contain 1"
        # At xi = pi, hat{kappa}_2(pi) = 1/2
        vp = kappa_2m_hat_arb(arb(pi, arb("1e-20")), m)
        assert vp.is_finite(), f"non-finite at removable singularity pi for m={m}"
    v_pi_m1 = kappa_2m_hat_arb(pi, 1)
    assert v_pi_m1.contains(arb("0.5")), f"hat{{kappa}}_2(pi) = {v_pi_m1}, expected 0.5"
    # Evenness, at and near the removable singularities +-k*pi and at generic points (review, 2026-10-30:
    # the negative-frequency branch used to return -1/2 at -pi)
    for m in [1, 2, 3, 40]:
        for k in range(0, min(m, 3) + 2):
            for x in (k * pi, k * pi + arb("1e-12"), k * pi + arb("0.3"), arb(pi * k, arb("1e-20"))):
                vp, vm = kappa_2m_hat_arb(x, m), kappa_2m_hat_arb(-x, m)
                assert vp.overlaps(vm), f"hat{{kappa}}_{2*m} not even at {x}: {vp} vs {vm}"
    for m in [1, 2, 3]:
        for k in range(1, m + 1):     # closed form: hat{kappa}_{2m}(k pi) = (-1)^(k-1) k^2 prod_{j != k} j^2/(j^2-k^2) / (2 k^2)
            exact = arb((-1) ** (k - 1)) / 2
            for j in range(1, m + 1):
                if j != k:
                    exact *= arb(j * j) / arb(j * j - k * k)
            for s in (1, -1):
                v = kappa_2m_hat_arb(s * k * pi, m)
                assert v.overlaps(exact), f"hat{{kappa}}_{2*m}({s}*{k}pi) = {v}, expected {exact}"

    # 4. Hypothesis checks in total_error
    omega = arb("727.9513358")
    eta1 = arb("1.55e-6")
    delta = arb("1e-8")
    T0 = arb("6800000")
    # Negative weight rejected
    must_fail(lambda: total_error(omega, eta1, delta, T0, arb("-0.1"), arb("0.6"), arb("0.5"), verbose=False))
    # Sum != 1 rejected
    must_fail(lambda: total_error(omega, eta1, delta, T0, arb("0.1"), arb("0.7"), arb("0.1"), verbose=False))
    # xi_T <= 3*pi rejected
    must_fail(lambda: total_error(omega, arb("1e-6"), delta, arb("3000000"), arb("0.12"), arb("0.70"), arb("0.18"), verbose=False))

    # 5. require() also triggers under python -O
    r = subprocess.run([sys.executable, "-O", "-c",
                        f"import sys; sys.path.insert(0, {HERE!r}); "
                        "from certio import require, fail, CertificationError\n"
                        "try:\n require(False, 'boom')\nexcept CertificationError as e:\n fail(str(e))"],
                       capture_output=True, text=True)
    assert r.returncode != 0 and "FAIL: boom" in r.stderr

    # 6. Directed decimal rounding: exact in both directions, for both signs (a nearest-rounded exponent such as
    #    47.9376 for 47.93760276 is below the bound it stands for)
    u = arb("110.3804095") / arb(10).log()
    assert directed_decimal(u, 4, upper=True) == "47.9377" and directed_decimal(u, 4) == "47.9376"
    assert directed_decimal(arb("-0.00019631"), 6, upper=True) == "-0.000196"
    assert directed_decimal(arb("-0.00019631"), 6) == "-0.000197"
    assert directed_decimal(arb("0.5"), 3, upper=True) == "0.500" and directed_decimal(arb(7), 0) == "7"
    rng = random.Random(1)
    for _ in range(2000):
        v = rng.uniform(-50, 50) * 10 ** rng.randint(-12, 3)
        x = arb(v, abs(v) * rng.choice([0, 1e-9, 1e-3]))
        places = rng.randint(0, 12)
        for up in (False, True):
            got = Fraction(directed_decimal(x, places, upper=up))
            m, e = (int(t) for t in (x.upper() if up else x.lower()).man_exp())
            end = Fraction(m) * Fraction(2) ** e
            assert (0 <= got - end < Fraction(1, 10 ** places)) if up else (0 <= end - got < Fraction(1, 10 ** places))
    must_fail(lambda: directed_decimal(arb("nan"), 3))
    must_fail(lambda: ball(arb("inf")))

    # 7. Finer R4 cells only tighten the bound (every cell bound is an upper bound)
    l = (arb("0.12"), arb("0.70"), arb("0.18"))
    for f in (bound_R4_envelope, bound_R4_stieltjes):
        coarse, fine = f(omega, eta1, T0, *l, n_sub=512), f(omega, eta1, T0, *l, n_sub=4096)
        assert fine < coarse, f"{f.__name__}: finer cells did not tighten the bound"

    print("test_certio.py: ALL TESTS PASSED.")


if __name__ == "__main__":
    main()
