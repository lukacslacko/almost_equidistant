"""Minimal rigorous interval arithmetic with outward rounding, used only to
construct the exact seed-simplex enclosures handed to the C kernel (which
contains the full certified interval engine).
IEEE-754 basic ops (+,-,*,/,sqrt) are correctly rounded, so nextafter-outward
gives true enclosures."""
import math

INF = math.inf

def _lo(x): return math.nextafter(x, -INF)
def _hi(x): return math.nextafter(x, INF)

class IV:
    __slots__ = ('a', 'b')
    def __init__(self, a, b=None):
        if b is None: b = a
        self.a = a; self.b = b
    def __repr__(self): return f"[{self.a!r},{self.b!r}]"
    def __truediv__(s, o):
        o = _c(o)
        if o.a <= 0 <= o.b: raise ZeroDivisionError("interval contains 0")
        p = (s.a / o.a, s.a / o.b, s.b / o.a, s.b / o.b)
        return IV(_lo(min(p)), _hi(max(p)))

def _c(x):
    return x if isinstance(x, IV) else IV(float(x))

def isqrt(x):
    x = _c(x)
    if x.b < 0: raise ValueError("sqrt of negative interval")
    a = max(x.a, 0.0)
    return IV(_lo(math.sqrt(a)), _hi(math.sqrt(x.b)))
