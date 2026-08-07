"""Escalation for frontier survivors: for each surviving level-(n+1) graph C
and every vertex u, enumerate ALL maximal triangle-free completions of
comp(C-u); each complement must be a residue member (certified graphs kill
C).  Branching rule: for the first addable pair p=(a,b), a maximal
completion either contains p, or contains a blocking cherry a-z-b; both
branches are explored, duplicates deduped at the leaves.  Canonical forms
are computed by ext6 canon (the C tool), used as an external oracle.

Usage: python escalate_frontier.py <frontier_out> <residue_corpus>
"""
import subprocess, sys, tempfile, os

EXT6 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ext6.exe")

def canon_batch(graphs):
    """graphs: list of (n, masks). Returns list of canonical line strings."""
    if not graphs:
        return []
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False,
                                     dir=os.path.dirname(EXT6)) as tf:
        for n, masks in graphs:
            tf.write(str(n) + " " + " ".join(str(m) for m in masks) + "\n")
        name = tf.name
    try:
        out = subprocess.run([EXT6, "canon", name], capture_output=True,
                             text=True, check=True).stdout
    finally:
        os.unlink(name)
    return out.strip().splitlines()

def addable(H, n, a, b):
    return not (H[a] >> b) & 1 and not (H[a] & H[b])

def first_free(H, n):
    for a in range(n):
        for b in range(a + 1, n):
            if addable(H, n, a, b):
                return (a, b)
    return None

def all_max_completions(H, n, cap=200000):
    """Yield edge-set-frozensets of all maximal TF supergraphs of H."""
    seen = set()
    nodes = 0
    def rec(H):
        nonlocal nodes
        nodes += 1
        if nodes > cap:
            raise RuntimeError("completion enumeration cap exceeded")
        p = first_free(H, n)
        if p is None:
            key = tuple(H)
            if key not in seen:
                seen.add(key)
                yield H[:]
            return
        a, b = p
        H2 = H[:]
        H2[a] |= 1 << b; H2[b] |= 1 << a
        yield from rec(H2)
        for z in range(n):
            if z == a or z == b:
                continue
            H3 = H[:]
            ok = True
            for (x, y) in ((a, z), (b, z)):
                if (H3[x] >> y) & 1:
                    continue
                if H3[x] & H3[y]:
                    ok = False
                    break
                H3[x] |= 1 << y; H3[y] |= 1 << x
            if ok and not (H3[a] >> b) & 1 and not addable(H3, n, a, b):
                yield from rec(H3)
    yield from rec(H)

def main():
    frontier_fn, residue_fn = sys.argv[1], sys.argv[2]
    res_lines = open(residue_fn).read().split("\n")
    res_graphs = []
    for line in res_lines:
        p = line.split()
        if p:
            res_graphs.append((int(p[0]), [int(x) for x in p[1:]]))
    residue_canon = set(canon_batch(res_graphs))
    rn = res_graphs[0][0]
    print(f"residue: {len(residue_canon)} canonical forms on {rn} vertices")

    survivors = []
    for line in open(frontier_fn):
        if line.startswith("SURVIVES"):
            p = line.split()
            n = int(p[2])
            survivors.append((int(p[1]), n, [int(x) for x in p[3:3 + n]]))
    print(f"escalating {len(survivors)} survivors")

    true_frontier = []
    for tag, n, C in survivors:
        assert n == rn + 1
        verdict = "SURVIVES"
        why = ""
        total_comps = 0
        for u in range(n):
            keep = [i for i in range(n) if i != u]
            Cd = []
            for i in keep:
                m = 0
                for jpos, j in enumerate(keep):
                    if (C[i] >> j) & 1:
                        m |= 1 << jpos
                Cd.append(m)
            nd = n - 1
            full = (1 << nd) - 1
            H = [full & ~Cd[i] & ~(1 << i) for i in range(nd)]
            comps = list(all_max_completions(H, nd))
            total_comps += len(comps)
            cgraphs = [(nd, [full & ~Hm[i] & ~(1 << i) for i in range(nd)])
                       for Hm in comps]
            for ci, cl in enumerate(canon_batch(cgraphs)):
                if cl not in residue_canon:
                    verdict = "KILLED"
                    why = f"del={u} completion#{ci} outside residue"
                    break
            if verdict == "KILLED":
                break
        print(f"graph {tag}: {verdict} {why} (completions tried {total_comps})")
        if verdict == "SURVIVES":
            true_frontier.append((tag, n, C))
    print(f"TRUE FRONTIER: {len(true_frontier)} graphs")
    with open("frontier20_true.txt", "w") as out:
        for tag, n, C in true_frontier:
            out.write(str(n) + " " + " ".join(str(m) for m in C) + "\n")

if __name__ == "__main__":
    main()
