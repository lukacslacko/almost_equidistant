"""Probe sub-arcs of a hot slice with the debug kernel, streaming results.
Usage: python probe_arc.py <gi> <lo_deg> <hi_deg> <parts> [cap]"""
import os, sys, ctypes, math
os.environ.setdefault("CK5_DUMP", "3")
import cdriver5
lib = ctypes.CDLL(os.path.join(cdriver5.HERE, "ckernel5dbg.dll"))
lib.decide5_c.restype = ctypes.c_int64
lib.decide5_c.argtypes = [
    ctypes.c_int, ctypes.POINTER(ctypes.c_uint32),
    ctypes.c_int, ctypes.POINTER(ctypes.c_int),
    ctypes.c_int, ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_int64,
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64)]
cdriver5._lib = lib
from cdriver5 import gen_orders, decide5, TWO_PI
from reproduce5 import load_level

gi = int(sys.argv[1]); lo0 = float(sys.argv[2]); hi0 = float(sys.argv[3])
parts = int(sys.argv[4]); cap = int(sys.argv[5]) if len(sys.argv) > 5 else 1_000_000
graphs = load_level(17)
adj = graphs[gi]
seed, order = gen_orders(adj, 17, kmax=12)[0]
deg = math.pi / 180
for k in range(parts):
    a = lo0 + k * (hi0 - lo0) / parts
    b = lo0 + (k + 1) * (hi0 - lo0) / parts
    st, nodes, unres = decide5(adj, 17, seed=seed, order=order,
                               th0=(a * deg, b * deg),
                               theta_min=TWO_PI / (1 << 8), max_nodes=cap)
    print(f"[{a:.4f},{b:.4f}] deg: {st} nodes={nodes} unres={unres}",
          flush=True)
