#!/usr/bin/env python3
"""Independent small controls for the exact d=6 graph-only filters.

This intentionally does not import or call profile_d6.c.  It reconstructs the
18-point Larman--Rogers-plus-two-apices lower-bound graph and checks every link
bound in plain Python.  It also checks synthetic witnesses for the K5-link and
K7-reflection rules.
"""
from itertools import combinations


BOUNDS = {2: 16, 3: 12, 4: 10, 5: 4, 6: 2, 7: 0}


def cliques(adj, k):
    n = len(adj)
    for q in combinations(range(n), k):
        if all((adj[u] >> v) & 1 for u, v in combinations(q, 2)):
            yield q


def link_violations(adj):
    out = []
    full = (1 << len(adj)) - 1
    for k, bound in BOUNDS.items():
        for q in cliques(adj, k):
            common = full
            for v in q:
                common &= adj[v]
            if common.bit_count() > bound:
                out.append((k, q, common.bit_count()))
                break
    return out


def has_k7_reflection_violation(adj):
    n = len(adj)
    for q in cliques(adj, 7):
        seed = sum(1 << v for v in q)
        forced = []
        for v in range(n):
            if (seed >> v) & 1:
                continue
            defects = seed & ~adj[v]
            if defects.bit_count() == 1:
                forced.append((v, defects.bit_length() - 1))
        for (u, tu), (v, tv) in combinations(forced, 2):
            if tu == tv or ((adj[u] >> v) & 1):
                return True
    return False


def lower_bound_18_graph():
    base = [x for x in range(32) if x.bit_count() % 2 == 1]
    adj = [0] * 18
    for i, x in enumerate(base):
        for j, y in enumerate(base[:i]):
            if (x ^ y).bit_count() == 2:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    for apex in (16, 17):
        for v in range(16):
            adj[apex] |= 1 << v
            adj[v] |= 1 << apex
    return adj


def k5_bad_link_graph():
    adj = [0] * 10
    for u, v in combinations(range(5), 2):
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    for u in range(5):
        for v in range(5, 10):
            adj[u] |= 1 << v
            adj[v] |= 1 << u
    for i in range(5):
        u, v = 5 + i, 5 + (i + 1) % 5
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    return adj


def reflection_bad_graph(collision=False):
    adj = [0] * 9
    for u, v in combinations(range(7), 2):
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    omitted = (0, 0 if collision else 1)
    for x, miss in zip((7, 8), omitted):
        for v in range(7):
            if v != miss:
                adj[x] |= 1 << v
                adj[v] |= 1 << x
    if not collision:
        adj[7] |= 1 << 8
        adj[8] |= 1 << 7
    return adj


def main():
    positive = lower_bound_18_graph()
    assert not link_violations(positive)
    assert not has_k7_reflection_violation(positive)
    assert any(True for _ in cliques(positive, 6))
    assert not any(True for _ in cliques(positive, 7))

    bad_link = link_violations(k5_bad_link_graph())
    assert any(k == 5 and count == 5 for k, _, count in bad_link)
    assert has_k7_reflection_violation(reflection_bad_graph(False))
    assert has_k7_reflection_violation(reflection_bad_graph(True))
    print("profile_d6 controls: PASS")


if __name__ == "__main__":
    main()
