"""Rigorous float interval arithmetic with outward rounding.
IEEE-754 basic ops (+,-,*,/,sqrt) are correctly rounded, so nextafter-outward
gives true enclosures. cos/sin get 8-ulp padding (libm is well within that)."""
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
    def __add__(s, o):
        o = _c(o); return IV(_lo(s.a + o.a), _hi(s.b + o.b))
    __radd__ = __add__
    def __neg__(s): return IV(-s.b, -s.a)
    def __sub__(s, o):
        o = _c(o); return IV(_lo(s.a - o.b), _hi(s.b - o.a))
    def __rsub__(s, o): return _c(o).__sub__(s)
    def __mul__(s, o):
        o = _c(o)
        p = (s.a * o.a, s.a * o.b, s.b * o.a, s.b * o.b)
        return IV(_lo(min(p)), _hi(max(p)))
    __rmul__ = __mul__
    def __truediv__(s, o):
        o = _c(o)
        if o.a <= 0 <= o.b: raise ZeroDivisionError("interval contains 0")
        p = (s.a / o.a, s.a / o.b, s.b / o.a, s.b / o.b)
        return IV(_lo(min(p)), _hi(max(p)))
    def __rtruediv__(s, o): return _c(o).__truediv__(s)
    def sq(s):
        if s.a >= 0: return IV(_lo(s.a * s.a), _hi(s.b * s.b))
        if s.b <= 0: return IV(_lo(s.b * s.b), _hi(s.a * s.a))
        m = max(-s.a, s.b); return IV(0.0, _hi(m * m))
    def contains0(s): return s.a <= 0 <= s.b
    def gt0(s): return s.a > 0
    def lt0(s): return s.b < 0
    def width(s): return s.b - s.a
    def mid(s): return 0.5 * (s.a + s.b)

def _c(x):
    return x if isinstance(x, IV) else IV(float(x))

def isqrt(x):
    x = _c(x)
    if x.b < 0: raise ValueError("sqrt of negative interval")
    a = max(x.a, 0.0)
    return IV(_lo(math.sqrt(a)), _hi(math.sqrt(x.b)))

_ULP8 = 8
def _pad(x):
    y = x
    for _ in range(_ULP8): y = math.nextafter(y, -INF)
    z = x
    for _ in range(_ULP8): z = math.nextafter(z, INF)
    return y, z

def icos(x):
    x = _c(x)
    if x.width() >= 2 * math.pi: return IV(-1.0, 1.0)
    ca, cb = math.cos(x.a), math.cos(x.b)
    lo = min(_pad(ca)[0], _pad(cb)[0]); hi = max(_pad(ca)[1], _pad(cb)[1])
    # check for interior extrema: k*pi inside [a,b]
    import math as m
    k0 = m.floor(x.a / m.pi); k1 = m.floor(x.b / m.pi)
    for k in range(k0, k1 + 1):
        t = k * m.pi
        if x.a <= t <= x.b:
            if k % 2 == 0: hi = 1.0
            else: lo = -1.0
    return IV(max(lo, -1.0), min(hi, 1.0))

def isin(x):
    x = _c(x)
    return icos(IV(math.pi / 2, math.pi / 2) - x)

# ---- vectors ----
def vadd(u, v): return [a + b for a, b in zip(u, v)]
def vsub(u, v): return [a - b for a, b in zip(u, v)]
def vscale(t, u): return [t * a for a in u]
def vdot(u, v):
    s = IV(0.0)
    for a, b in zip(u, v): s = s + a * b
    return s
def vnorm2(u):
    s = IV(0.0)
    for a in u: s = s + a.sq()
    return s
