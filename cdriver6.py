"""Driver for the certified d=6 C kernel (ckernel6.c). Mirror of
cdriver5.py with d=6 constants: seeds K_7 (unit 6-simplex) or K_6,
sphere placement needs >= 6 placed neighbours, circle placement
exactly 5. Seed enclosures verified: pairwise squared distances
certifiably contain 1, tighter than 1e-12."""
import ctypes, math, os, shutil, subprocess, sys
from ival import IV, isqrt

TWO_PI = 2 * math.pi
D = 6
HERE = os.path.dirname(os.path.abspath(__file__))

LIBNAME = {"win32": "ckernel6.dll", "darwin": "ckernel6.dylib"}.get(
    sys.platform, "ckernel6.so")

def _compile_kernel(src, out):
    if sys.platform != "win32":
        subprocess.check_call(["cc", "-O2", "-shared", "-o", out, src, "-lm"],
                              cwd=HERE)
        return
    llvm = r"C:\Program Files\LLVM\bin\clang.exe"
    for gcc in ("cc", "gcc", "clang",
                llvm if os.path.exists(llvm) else None):
        if gcc and shutil.which(gcc):
            cmd = [gcc, "-O2", "-shared", "-o", out, src]
            if "clang" in os.path.basename(gcc):
                cmd.insert(1, "-fuse-ld=lld")
            subprocess.check_call(cmd, cwd=HERE)
            return
    raise RuntimeError("no C compiler found")

def ensure_kernel():
    lib = os.path.join(HERE, LIBNAME)
    src = os.path.join(HERE, "ckernel6.c")
    if (not os.path.exists(lib)
            or os.path.getmtime(lib) < os.path.getmtime(src)):
        print("compiling d=6 C kernel ...", flush=True)
        _compile_kernel("ckernel6.c", LIBNAME)

_lib = None
def _kernel():
    global _lib
    if _lib is None:
        ensure_kernel()
        _lib = ctypes.CDLL(os.path.join(HERE, LIBNAME))
        _lib.decide6_c.restype = ctypes.c_int64
        _lib.decide6_c.argtypes = [
            ctypes.c_int, ctypes.POINTER(ctypes.c_uint32),
            ctypes.c_int, ctypes.POINTER(ctypes.c_int),
            ctypes.c_int, ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_double, ctypes.c_double, ctypes.c_double,
            ctypes.c_int64,
            ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),
        ]
    return _lib

# unit-edge regular simplex in R^6, exact interval enclosures
def seed_coords(k):
    r3 = isqrt(IV(3.0)); r6 = isqrt(IV(6.0)); r10 = isqrt(IV(10.0))
    r15 = isqrt(IV(15.0)); r21 = isqrt(IV(21.0))
    Z = IV(0.0)
    pts = [
        [Z, Z, Z, Z, Z, Z],
        [IV(1.0), Z, Z, Z, Z, Z],
        [IV(0.5), r3 / 2, Z, Z, Z, Z],
        [IV(0.5), r3 / 6, r6 / 3, Z, Z, Z],
        [IV(0.5), r3 / 6, r6 / 12, r10 / 4, Z, Z],
        [IV(0.5), r3 / 6, r6 / 12, r10 / 20, r15 / 5, Z],
        [IV(0.5), r3 / 6, r6 / 12, r10 / 20, r15 / 30, r21 / 6],
    ]
    pts = pts[:k]
    for i in range(k):
        for j in range(i + 1, k):
            lo = hi = 0.0
            for t in range(D):
                dlo = pts[i][t].a - pts[j][t].b
                dhi = pts[i][t].b - pts[j][t].a
                m = max(abs(dlo), abs(dhi))
                mn = 0.0 if dlo <= 0 <= dhi else min(abs(dlo), abs(dhi))
                lo += mn * mn; hi += m * m
            assert lo <= 1.0 <= hi and hi - lo < 1e-12, (i, j, lo, hi)
    return pts

def _cliques(adj, n, k):
    res = []
    def ext(clq, cand, start):
        if len(clq) == k: res.append(tuple(clq)); return
        for v in range(start, n):
            if (cand >> v) & 1: ext(clq + [v], cand & adj[v], v + 1)
    ext([], (1 << n) - 1, 0)
    return res

def _layered_order(adj, n, seed, pick, force_circle_first=False):
    placed = set(seed); order = []
    forced_done = not force_circle_first
    while len(placed) < n:
        if not forced_done:
            c5 = [v for v in range(n) if v not in placed
                  and sum(1 for u in placed if (adj[v] >> u) & 1) == D - 1]
            if c5:
                v = pick(c5, placed); order.append(v); placed.add(v)
                forced_done = True
                continue
        progressed = False
        while True:
            cands = [(sum(1 for u in placed if (adj[v] >> u) & 1), v)
                     for v in range(n) if v not in placed]
            cands = [(k, v) for k, v in cands if k >= D]
            if not cands: break
            k, v = max(cands); order.append(v); placed.add(v)
            progressed = True
            if not forced_done:
                break
        if len(placed) == n: break
        if not forced_done and progressed:
            continue
        c5 = [v for v in range(n) if v not in placed
              and sum(1 for u in placed if (adj[v] >> u) & 1) == D - 1]
        if not c5: return None
        v = pick(c5, placed); order.append(v); placed.add(v)
        forced_done = True
    return order

def gen_orders(adj, n, kmax=12):
    seeds = _cliques(adj, n, 7) or _cliques(adj, n, 6)
    def closure_size(w, placed):
        P = set(placed); P.add(w); grew = True
        while grew:
            grew = False
            for x in range(n):
                if x not in P and sum(1 for u in P if (adj[x] >> u) & 1) >= D:
                    P.add(x); grew = True
        return len(P)
    picks = [lambda c, p: max(c, key=lambda w: closure_size(w, p)),
             lambda c, p: min(c, key=lambda w: closure_size(w, p)),
             lambda c, p: max(c), lambda c, p: min(c)]
    outs, seen = [], set()
    for seed in seeds:
        for pick in picks:
            for fc in (False, True):
                o = _layered_order(adj, n, seed, pick, force_circle_first=fc)
                if o is not None and (seed, tuple(o)) not in seen:
                    seen.add((seed, tuple(o))); outs.append((seed, o))
    def score(so):
        seed, order = so
        placed = set(seed); circ = []
        for i, v in enumerate(order):
            if sum(1 for u in placed if (adj[v] >> u) & 1) == D - 1:
                circ.append(i)
            placed.add(v)
        return (len(circ), [-c for c in circ])
    outs.sort(key=score)
    if len(outs) <= kmax:
        return outs
    zero = [so for so in outs if score(so)[0] == 0]
    circ = [so for so in outs if score(so)[0] > 0]
    def first_circle_pos(so):
        seed, order = so
        placed = set(seed)
        for i, v in enumerate(order):
            if sum(1 for u in placed if (adj[v] >> u) & 1) == D - 1:
                return i
            placed.add(v)
        return 1 << 30
    circ.sort(key=first_circle_pos)
    keep = zero[:kmax - min(len(circ), kmax // 2)] + circ[:kmax // 2]
    return keep[:kmax] if keep else outs[:kmax]

def decide6(adj, n, seed=None, order=None, theta_min=TWO_PI / (1 << 18),
            th0=(0.0, TWO_PI), max_nodes=500_000_000):
    if seed is None:
        seed, order = gen_orders(adj, n, kmax=1)[0]
    sc = seed_coords(len(seed))
    scarr = (ctypes.c_double * (len(seed) * D * 2))()
    for i in range(len(seed)):
        for j in range(D):
            scarr[(i * D + j) * 2] = sc[i][j].a
            scarr[(i * D + j) * 2 + 1] = sc[i][j].b
    adjarr = (ctypes.c_uint32 * n)(*adj)
    seedarr = (ctypes.c_int * len(seed))(*seed)
    orderarr = (ctypes.c_int * len(order))(*order)
    nodes = ctypes.c_int64(); unres = ctypes.c_int64()
    st = _kernel().decide6_c(n, adjarr, len(seed), seedarr,
                             len(order), orderarr,
                             scarr, theta_min, th0[0], th0[1], max_nodes,
                             ctypes.byref(nodes), ctypes.byref(unres))
    status = {0: "KILLED", 1: "SURVIVORS", 2: "ABORT"}[st]
    return status, nodes.value, unres.value
