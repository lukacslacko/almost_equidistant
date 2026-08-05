"""Driver for the certified C kernel (ckernel.c).

This module is the only Python<->C bridge: it compiles the kernel on first
use, constructs the exact interval enclosures of the seed simplex (the only
rigorous arithmetic done on the Python side; everything else is in C), and
generates the elimination orders that the search follows.
"""
import ctypes, math, os, shutil, subprocess, sys
from ival import IV, isqrt

TWO_PI = 2 * math.pi
HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# kernel compilation and loading
# ---------------------------------------------------------------------------
LIBNAME = {"win32": "ckernel.dll", "darwin": "ckernel.dylib"}.get(
    sys.platform, "ckernel.so")

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
    for junk in ("ckernel.obj", "ckernel.lib", "ckernel.exp"):
        try: os.remove(os.path.join(HERE, junk))
        except OSError: pass

def ensure_kernel():
    lib = os.path.join(HERE, LIBNAME)
    src = os.path.join(HERE, "ckernel.c")
    if (not os.path.exists(lib)
            or os.path.getmtime(lib) < os.path.getmtime(src)):
        print("compiling C kernel ...", flush=True)
        _compile_kernel("ckernel.c", LIBNAME)

_lib = None
def _kernel():
    global _lib
    if _lib is None:
        ensure_kernel()
        _lib = ctypes.CDLL(os.path.join(HERE, LIBNAME))
        _lib.decide_c.restype = ctypes.c_int64
        _lib.decide_c.argtypes = [
            ctypes.c_int, ctypes.POINTER(ctypes.c_ushort),
            ctypes.c_int, ctypes.POINTER(ctypes.c_int),
            ctypes.c_int, ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_double, ctypes.c_double, ctypes.c_double,
            ctypes.c_int64,
            ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),
        ]
    return _lib

# ---------------------------------------------------------------------------
# seed simplex: unit-edge regular (k-1)-simplex, exact interval enclosures
# ---------------------------------------------------------------------------
def seed_coords(k):
    r3 = isqrt(IV(3.0)); r6 = isqrt(IV(6.0)); r10 = isqrt(IV(10.0))
    pts = [
        [IV(0.0), IV(0.0), IV(0.0), IV(0.0)],
        [IV(1.0), IV(0.0), IV(0.0), IV(0.0)],
        [IV(0.5), r3 / 2, IV(0.0), IV(0.0)],
        [IV(0.5), r3 / 6, r6 / 3, IV(0.0)],
        [IV(0.5), r3 / 6, r6 / 12, r10 / 4],
    ]
    return pts[:k]

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
            cands = [(k, v) for k, v in cands if k >= 4]
            if not cands: break
            k, v = max(cands); order.append(v); placed.add(v)
        if len(placed) == n: break
        c3 = [v for v in range(n) if v not in placed
              and sum(1 for u in placed if (adj[v] >> u) & 1) == 3]
        if not c3: return None
        v = pick(c3, placed); order.append(v); placed.add(v)
    return order

def gen_orders(adj, n, kmax=8):
    """Candidate (seed, order) decompositions, best first: fewest circle
    steps, then circle steps as late as possible."""
    seeds = _cliques(adj, n, 5) or _cliques(adj, n, 4)
    def closure_size(w, placed):
        P = set(placed); P.add(w); grew = True
        while grew:
            grew = False
            for x in range(n):
                if x not in P and sum(1 for u in P if (adj[x] >> u) & 1) >= 4:
                    P.add(x); grew = True
        return len(P)
    picks = [lambda c3, p: max(c3, key=lambda w: closure_size(w, p)),
             lambda c3, p: min(c3, key=lambda w: closure_size(w, p)),
             lambda c3, p: max(c3), lambda c3, p: min(c3)]
    outs, seen = [], set()
    for seed in seeds:
        for pick in picks:
            o = _layered_order(adj, n, seed, pick)
            if o and (seed, tuple(o)) not in seen:
                seen.add((seed, tuple(o))); outs.append((seed, o))
    def score(so):
        seed, order = so
        placed = set(seed); circ = []
        for i, v in enumerate(order):
            if sum(1 for u in placed if (adj[v] >> u) & 1) == 3: circ.append(i)
            placed.add(v)
        return (len(circ), [-c for c in circ])
    outs.sort(key=score)
    return outs[:kmax]

# ---------------------------------------------------------------------------
# entry point
# ---------------------------------------------------------------------------
def decide_c(adj, n, seed=None, order=None, theta_min=TWO_PI / (1 << 18),
             th0=(0.0, TWO_PI), max_nodes=500_000_000):
    if seed is None:
        seed, order = gen_orders(adj, n, kmax=1)[0]
    sc = seed_coords(len(seed))
    scarr = (ctypes.c_double * (len(seed) * 4 * 2))()
    for i in range(len(seed)):
        for j in range(4):
            scarr[(i * 4 + j) * 2] = sc[i][j].a
            scarr[(i * 4 + j) * 2 + 1] = sc[i][j].b
    adjarr = (ctypes.c_ushort * n)(*adj)
    seedarr = (ctypes.c_int * len(seed))(*seed)
    orderarr = (ctypes.c_int * len(order))(*order)
    nodes = ctypes.c_int64(); unres = ctypes.c_int64()
    st = _kernel().decide_c(n, adjarr, len(seed), seedarr, len(order), orderarr,
                            scarr, theta_min, th0[0], th0[1], max_nodes,
                            ctypes.byref(nodes), ctypes.byref(unres))
    status = {0: "KILLED", 1: "SURVIVORS", 2: "ABORT"}[st]
    return status, nodes.value, unres.value
