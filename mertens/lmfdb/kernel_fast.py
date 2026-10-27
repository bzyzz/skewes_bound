"""Fast rigorous hat{kappa}_{2M}(xi) for small |xi| (Arb).

hat kappa_{2M}(xi) = sinc(xi) * prod_{k<=M} (1 - xi^2/(k pi)^2)^{-1}
                   = sinc(xi) * exp( sum_{j>=1} (xi^2/pi^2)^j H_M(2j)/j ),  H_M(s) = sum_{k<=M} k^{-s},
valid for |xi| < pi. Truncating after J terms leaves a tail of at most sum_{j>J} q^j zeta(2j)/j <= 2 q^{J+1}/((J+1)(1-q)),
with q = xi^2/pi^2 <= 1/4 (since zeta(2j) <= zeta(2) < 2). The tail is added as a radius to the exponent.
"""
from flint import arb

class Mollifier:
    def __init__(self, M, J=10):
        self.M, self.J = M, J
        pi2 = arb.pi()**2
        self.c = []
        for j in range(1, J + 1):
            H = sum((arb(1)/arb(k)**(2*j) for k in range(1, M + 1)), arb(0))
            self.c.append(H/(j*pi2**j))           # coefficient of xi^{2j}
    def __call__(self, xi):
        x2 = xi*xi
        q = (x2/arb.pi()**2)
        if not (q < arb("0.25")):
            raise ValueError("mollifier argument too large for the series")
        s = arb(0); p = arb(1)
        for c in self.c:
            p = p*x2; s += c*p
        J = self.J
        tail = (2*q**(J + 1)/((J + 1)*(1 - q))).upper()
        s = s + arb(0, tail)
        return xi.sinc()*s.exp()
