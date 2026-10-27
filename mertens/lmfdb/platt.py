"""Python wrapper for libplatt.so: exact decoding of Platt / LMFDB zeta zeros."""
import ctypes, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
_lib = ctypes.CDLL(os.path.join(HERE, "libplatt.so"))
_lib.platt_count.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_long)]
_lib.platt_count.restype = ctypes.c_long
P = lambda dt: np.ctypeslib.ndpointer(dtype=dt, flags="C_CONTIGUOUS")
_lib.platt_decode.argtypes = [ctypes.c_char_p, P(np.float64), P(np.uint64), P(np.uint64), P(np.uint64), ctypes.c_long,
                              P(np.float64), P(np.float64), P(np.uint64), P(np.uint64), ctypes.c_long]
_lib.platt_decode.restype = ctypes.c_long
EPS_EXP = -101      # gamma = t0 + Z * 2^-101
RAD_EXP = -101      # stored to +- 2^-102 (Platt 2017); we attach 2^-101 for slack

def decode(path):
    nb = ctypes.c_long(0)
    n = _lib.platt_count(path.encode(), ctypes.byref(nb))
    if n < 0: raise IOError(f"cannot read {path}")
    t0z = np.empty(n); hi = np.empty(n, np.uint64); lo = np.empty(n, np.uint64); idx = np.empty(n, np.uint64)
    bt0 = np.empty(nb.value); bt1 = np.empty(nb.value); bN0 = np.empty(nb.value, np.uint64); bN1 = np.empty(nb.value, np.uint64)
    m = _lib.platt_decode(path.encode(), t0z, hi, lo, idx, n, bt0, bt1, bN0, bN1, nb.value)
    if m != n: raise IOError(f"decode error in {path}")
    return dict(t0=t0z, hi=hi, lo=lo, idx=idx, bt0=bt0, bt1=bt1, bN0=bN0, bN1=bN1)

def gamma_ld(d):
    """gamma as numpy longdouble (for floating-point passes)."""
    ld = np.longdouble
    return d["t0"].astype(ld) + d["hi"].astype(ld)*ld(2.0)**(64 + EPS_EXP) + d["lo"].astype(ld)*ld(2.0)**EPS_EXP
