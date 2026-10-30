#!/usr/bin/env python3
"""List the primitive Dirichlet characters of conductor 3 <= d <= DMAX as 'd n parity is_real' lines."""
import math, sys
from flint import dirichlet_char
DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30
for d in range(3, DMAX + 1):
    for n in range(1, d + 1):
        if math.gcd(n, d) == 1:
            c = dirichlet_char(d, n)
            if c.is_primitive():
                print(d, n, c.parity(), int(c.is_real()))
