"""
certio.py -- shared helpers for the Skewes bound certification scripts.

* Outward-safe serialisation: every certified quantity is written as an Arb ball string
  (``x.str(n, radius=True, more=True)``), checked to enclose the ball, and read back with
  ``arb(s)``; directed endpoints are taken after reading (``.lower()``, ``.upper()``).
* ``directed_decimal`` rounds an endpoint down or up to a decimal grid in exact integer arithmetic; every
  decimal bound that is stated as a bound goes through it.
* ``exact_int`` turns an Arb ball into an integer only if the ball contains exactly one.
* ``require`` / ``fail`` stop with a nonzero exit status; they are not ``assert`` statements,
  so they also run under ``python -O``.
* ``provenance`` / ``cache_tag`` bind caches and certificates to the parameters, the code and
  the library versions that produced them.
"""
import hashlib
import json
import os
import sys

import flint
from flint import arb

SCHEMA = 2


class CertificationError(RuntimeError):
    pass


def require(cond, msg):
    if not cond:
        raise CertificationError(msg)


def fail(msg):
    print(f"FAIL: {msg}", file=sys.stderr)
    sys.exit(1)


def versions():
    return {"python": sys.version.split()[0], "python_flint": flint.__version__,
            "flint": getattr(flint, "__FLINT_VERSION__", "unknown")}


def file_hash(*paths):
    h = hashlib.sha256()
    for p in paths:
        with open(p, "rb") as f:
            h.update(f.read())
    return h.hexdigest()


def cache_tag(params, code_paths):
    """Short tag identifying a cache namespace: parameters, schema, library versions, code."""
    blob = json.dumps({"params": params, "schema": SCHEMA, "versions": versions(),
                       "code": file_hash(*code_paths)}, sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest()[:16]


def provenance(params, code_paths):
    return {"schema": SCHEMA, "params": params, "versions": versions(),
            "code_sha256": {os.path.basename(p): file_hash(p) for p in code_paths}}


def ball(x, digits=30):
    """Ball string that encloses x, with all `digits` midpoint digits kept (more=True). Without more=True Arb drops
    the midpoint digits below the radius and widens the printed radius, so a sum read back from the string can be
    several times wider than the computed one (review item 3). The enclosure is checked, not assumed."""
    require(x.is_finite(), "cannot serialise a non-finite ball")
    s = x.str(digits, radius=True, more=True)
    require(arb(s).contains(x), "serialised ball does not enclose its input")
    return s


def directed_decimal(x, places, upper=False):
    """Decimal string with `places` digits after the point that is <= x.lower() (upper=False) or >= x.upper()
    (upper=True). The endpoint of a ball is an exact binary number m * 2^e; it is rounded to the decimal grid in
    integer arithmetic, so the direction is exact. Arb's str(n, radius=False) rounds to nearest and must not be
    used for a stated bound."""
    require(x.is_finite(), "cannot round a non-finite bound")
    require(isinstance(places, int) and places >= 0, "places must be a nonnegative integer")
    m, e = (int(v) for v in (x.upper() if upper else x.lower()).man_exp())
    num = m * 10 ** places
    if e >= 0:
        r = num << e
    else:
        den = 1 << -e
        r = -((-num) // den) if upper else num // den        # ceiling or floor, also for negative numbers
    if places == 0:
        return str(r)
    d = str(abs(r)).rjust(places + 1, "0")
    return f"{'-' if r < 0 else ''}{d[:-places]}.{d[-places:]}"


def exact_int(x, what):
    """The unique integer in the ball x; fails if the ball contains none or several."""
    lo, hi = x.lower().ceil(), x.upper().floor()
    require(lo == hi, f"{what}: ball {x.str(20)} does not determine a unique integer")
    n = lo.unique_fmpz()
    require(n is not None, f"{what}: no unique integer")
    return int(n)


def manifest(payloads):
    """sha256 over the sorted JSON payloads of all chunks that a summary is built from."""
    h = hashlib.sha256()
    for p in payloads:
        h.update(json.dumps(p, sort_keys=True).encode())
    return h.hexdigest()


def cache_root(default_dir):
    """Cache root: $CERT_CACHE if set, else the script directory."""
    return os.environ.get("CERT_CACHE", default_dir)
