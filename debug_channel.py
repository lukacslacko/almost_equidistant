"""Dump surviving leaves from a hot slice using the debug kernel
(ckernel5dbg.dll, CK5_DUMP env). Coarse theta floor so survivors park
fast. Usage: python debug_channel.py <gi> <k96>"""
import os, sys, ctypes, math
os.environ["CK5_DUMP"] = "3"
import cdriver5
from cdriver5 import gen_orders, seed_coords, TWO_PI
from reproduce5 import load_level

# point cdriver5 at the debug dll (bypass ensure_kernel/rebuild)
lib = ctypes.CDLL(os.path.join(cdriver5.HERE, "ckernel5dbg.dll"))
lib.decide5_c.restype = ctypes.c_int64
lib.decide5_c.argtypes = [
    ctypes.c_int, ctypes.POINTER(ctypes.c_uint32),
    ctypes.c_int, ctypes.POINTER(ctypes.c_int),
    ctypes.c_int, ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_int64,
    ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int64),
]
cdriver5._lib = lib   # decide5 now uses the debug kernel

gi = int(sys.argv[1]) if len(sys.argv) > 1 else 11402
k = int(sys.argv[2]) if len(sys.argv) > 2 else 20
graphs = load_level(17)
adj = graphs[gi]
decs = gen_orders(adj, 17, kmax=12)
seed, order = decs[0]
print(f"graph {gi} dec0 seed={seed} order={order}", flush=True)
lo, hi = k * TWO_PI / 96, (k + 1) * TWO_PI / 96
st, nodes, unres = cdriver5.decide5(adj, 17, seed=seed, order=order,
                                    th0=(lo, hi),
                                    theta_min=TWO_PI / (1 << 8),
                                    max_nodes=30_000_000)
print(f"slice: {st} nodes={nodes} unres={unres}", flush=True)
