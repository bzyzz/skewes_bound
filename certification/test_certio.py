#!/usr/bin/env python3
"""
Tests for certio.py and error_bounds.py failure modes and outward-rounding guarantees.
"""
import os
import subprocess
import sys
import flint
from flint import arb

from certio import CertificationError, ball, exact_int, require
from error_bounds import kappa_2m_hat_arb, total_error

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

    print("test_certio.py: ALL TESTS PASSED.")


if __name__ == "__main__":
    main()
