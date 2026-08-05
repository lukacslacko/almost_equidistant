"""Driver for the certified d=5 C kernel (ckernel5.c).

Mirror of cdriver.py (the f(4)=12 driver) with the d=5 constants:
seeds are K_6 (unit 5-simplex) or K_5, a sphere placement needs >= 5
placed neighbours, a circle placement exactly 4. The seed enclosures are
constructed with the same minimal rigorous interval arithmetic (ival.py)
and verified: every pairwise squared distance interval must contain 1
and be tighter than 1e-12.
"""
import ctypes, math, os, shutil, subprocess, sys
from ival import IV, isqrt

TWO_PI = 2 * math.pi
D = 5
HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# kernel compilation and loading
# ---------------------------------------------------------------------------
LIBNAME = {"win32": "ckernel5.dll", "darwin": "ckernel5.dylib"}.get(
    sys.platform, "ckernel5.so")

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
                cmd.insert(1, "-fuse-ld=lld")   # link.exe is not on PATH
            subprocess.check_call(cmd, cwd=HERE)
            return
    # MSVC: locate via vswhere, compile from a vcvars64 environment
    vswhere = os.path.join(os.environ.get("ProgramFiles(x86)",
                                          r"C:\Program Files (x86)"),
                           "Microsoft Visual Studio", "Installer",
                           "vswhere.exe")
    vsroot = subprocess.check_output(
        [vswhere, "-latest", "-products", "*", "-requires",
         "Microsoft.VisualStudio.Component.VC.Tools.x86.x64",
         "-property", "installationPath"], text=True).strip()
    if not vsroot:
        raise RuntimeError("no C compiler found (need cc/gcc/clang or MSVC)")
    vcvars = os.path.join(vsroot, "VC", "Auxiliary", "Build", "vcvars64.bat")
    env = dict(os.environ)
    env["PATH"] += os.pathsep + os.path.dirname(vswhere)  # vcvars uses vswhere
    subprocess.check_call(
        f'"{vcvars}" >NUL && cl /nologo /O2 /LD {src} /Fe:{out}',
        shell=True, cwd=HERE, env=env)
    for junk in ("ckernel5.obj", "ckernel5.lib", "ckernel5.exp"):
        try: os.remove(os.path.join(HERE, junk))
        except OSError: pass

def ensure_kernel():
    lib = os.path.join(HERE, LIBNAME)
    src = os.path.join(HERE, "ckernel5.c")
    if (not os.path.exists(lib)
            or os.path.getmtime(lib) < os.path.getmtime(src)):
        print("compiling d=5 C kernel ...", flush=True)
        _compile_kernel("ckernel5.c", LIBNAME)

_lib = None
def _kernel():
    global _lib
    if _lib is None:
        ensure_kernel()
        _lib = ctypes.CDLL(os.path.join(HERE, LIBNAME))
        _lib.decide5_c.restype = ctypes.c_int64
        _lib.decide5_c.argtypes = [
            ctypes.c_int, ctypes.POINTER(ctypes.c_uint32),
            ctypes.c_int, ctypes.POINTER(ctypes.c_int),
            ctypes.c_int, ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_double, ctypes.c_double, ctypes.c_double,
            ctypes.c_int64,
            ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),
        ]
    return _lib

# ---------------------------------------------------------------------------
# seed simplex: unit-edge regular (k-1)-simplex in R^5, exact enclosures
# ---------------------------------------------------------------------------
def seed_coords(k):
    r3 = isqrt(IV(3.0)); r6 = isqrt(IV(6.0))
    r10 = isqrt(IV(10.0)); r15 = isqrt(IV(15.0))
    Z = IV(0.0)
    pts = [
        [Z, Z, Z, Z, Z],
        [IV(1.0), Z, Z, Z, Z],
        [IV(0.5), r3 / 2, Z, Z, Z],
        [IV(0.5), r3 / 6, r6 / 3, Z, Z],
        [IV(0.5), r3 / 6, r6 / 12, r10 / 4, Z],
        [IV(0.5), r3 / 6, r6 / 12, r10 / 20, r15 / 5],
    ]
    pts = pts[:k]
    # verify: all pairwise squared distances certifiably ~ 1
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

# ---------------------------------------------------------------------------
# elimination orders (layered, circles as late as possible; several tie-breaks)
# ---------------------------------------------------------------------------
def _cliques(adj, n, k):
    res = []
    def ext(clq, cand, start):
        if len(clq) == k: res.append(tuple(clq)); return
        for v in range(start, n):
            if (cand >> v) & 1: ext(clq + [v], cand & adj[v], v + 1)
    ext([], (1 << n) - 1, 0)
    return res

def _layered_order(adj, n, seed, pick):
    placed = set(seed); order = []
    while len(placed) < n:
        while True:
            cands = [(sum(1 for u in placed if (adj[v] >> u) & 1), v)
                     for v in range(n) if v not in placed]
            cands = [(k, v) for k, v in cands if k >= D]
            if not cands: break
            k, v = max(cands); order.append(v); placed.add(v)
        if len(placed) == n: break
        c4 = [v for v in range(n) if v not in placed
              and sum(1 for u in placed if (adj[v] >> u) & 1) == D - 1]
        if not c4: return None
        v = pick(c4, placed); order.append(v); placed.add(v)
    return order

def gen_orders(adj, n, kmax=8):
    """Candidate (seed, order) decompositions, best first: fewest circle
    steps, then circle steps as late as possible."""
    seeds = _cliques(adj, n, 6) or _cliques(adj, n, 5)
    def closure_size(w, placed):
        P = set(placed); P.add(w); grew = True
        while grew:
            grew = False
            for x in range(n):
                if x not in P and sum(1 for u in P if (adj[x] >> u) & 1) >= D:
                    P.add(x); grew = True
        return len(P)
    picks = [lambda c4, p: max(c4, key=lambda w: closure_size(w, p)),
             lambda c4, p: min(c4, key=lambda w: closure_size(w, p)),
             lambda c4, p: max(c4), lambda c4, p: min(c4)]
    outs, seen = [], set()
    for seed in seeds:
        for pick in picks:
            o = _layered_order(adj, n, seed, pick)
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
    return outs[:kmax]

# ---------------------------------------------------------------------------
# entry point
# ---------------------------------------------------------------------------
def decide5(adj, n, seed=None, order=None, theta_min=TWO_PI / (1 << 18),
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
    st = _kernel().decide5_c(n, adjarr, len(seed), seedarr,
                             len(order), orderarr,
                             scarr, theta_min, th0[0], th0[1], max_nodes,
                             ctypes.byref(nodes), ctypes.byref(unres))
    status = {0: "KILLED", 1: "SURVIVORS", 2: "ABORT"}[st]
    return status, nodes.value, unres.value
