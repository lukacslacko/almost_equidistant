"""ctypes driver for the C certified kernel."""
import ctypes, json, math, os, time
from decide import seed_coords, make_order2

TWO_PI = 2 * math.pi
_here = os.path.dirname(os.path.abspath(__file__))
lib = ctypes.CDLL(os.path.join(_here, "ckernel.dylib"))
lib.decide_c.restype = ctypes.c_long
lib.decide_c.argtypes = [
    ctypes.c_int, ctypes.POINTER(ctypes.c_ushort),
    ctypes.c_int, ctypes.POINTER(ctypes.c_int),
    ctypes.c_int, ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_long,
    ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long),
]

def decide_c(adj, n, seed=None, order=None, theta_min=TWO_PI / (1 << 18),
             th0=(0.0, TWO_PI), max_nodes=500_000_000):
    if seed is None:
        seed, order, _ = make_order2(adj, n)
    sc = seed_coords(len(seed))
    scarr = (ctypes.c_double * (len(seed) * 4 * 2))()
    for i in range(len(seed)):
        for j in range(4):
            scarr[(i * 4 + j) * 2] = sc[i][j].a
            scarr[(i * 4 + j) * 2 + 1] = sc[i][j].b
    adjarr = (ctypes.c_ushort * n)(*adj)
    seedarr = (ctypes.c_int * len(seed))(*seed)
    orderarr = (ctypes.c_int * len(order))(*order)
    nodes = ctypes.c_long(); unres = ctypes.c_long()
    st = lib.decide_c(n, adjarr, len(seed), seedarr, len(order), orderarr,
                      scarr, theta_min, th0[0], th0[1], max_nodes,
                      ctypes.byref(nodes), ctypes.byref(unres))
    status = {0: "KILLED", 1: "SURVIVORS", 2: "ABORT"}[st]
    return status, nodes.value, unres.value

if __name__ == "__main__":
    import sys
    # validation: G10 and cross-polytope
    n = 10
    E = set()
    for i in range(8):
        for j in range(i + 1, 8):
            if j != i + 4: E.add((i, j))
    E.discard((2, 4))
    for u in (0, 1, 2, 3): E.add(tuple(sorted((8, u))))
    for u in (2, 3, 4, 5): E.add(tuple(sorted((9, u))))
    E.add((8, 9))
    adj = [0] * n
    for i, j in E: adj[i] |= 1 << j; adj[j] |= 1 << i
    t = time.time()
    print("G10 (expect KILLED):", decide_c(adj, n), f"{time.time()-t:.1f}s")

    n2 = 8
    adj3 = [0] * n2
    for i in range(8):
        for j in range(8):
            if i != j and j != (i + 4) % 8 and i != (j + 4) % 8:
                adj3[i] |= 1 << j
    t = time.time()
    print("crosspoly (expect SURVIVORS):",
          decide_c(adj3, n2, theta_min=TWO_PI / (1 << 12)), f"{time.time()-t:.1f}s")

    # full 59 comparison with Python verdicts
    data = json.load(open("aeq_d4_n13.json"))
    graphs = data["graphs"]
    t0 = time.time()
    summary = {}
    for gi, a in enumerate(graphs):
        t = time.time()
        st, nodes, unres = decide_c(a, 13, max_nodes=200_000_000)
        summary[gi] = st
        print(f"[{gi:2d}] {st} nodes={nodes} unres={unres} t={time.time()-t:.1f}s",
              flush=True)
    print("total", time.time() - t0)
    killed = [k for k, v in summary.items() if v == "KILLED"]
    print(f"KILLED {len(killed)}/59; not killed:",
          {k: v for k, v in summary.items() if v != "KILLED"})
    json.dump(summary, open("c_results.json", "w"))
