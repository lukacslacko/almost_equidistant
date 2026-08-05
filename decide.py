"""Certified branch-and-prune deciding unit-distance realizability in R^4.

For a graph G, a realization assigns points in R^4 with |x_i - x_j| = 1 for
every edge. We place a seed clique (K4/K5, canonical coordinates, exact
interval enclosures), then follow an elimination order:
  - vertex with >= 4 placed neighbours: intersection of >= 4 unit spheres ->
    at most 2 candidate points (branch), each an interval enclosure of the
    exact solutions; remaining placed-neighbour constraints are kill-checks.
  - vertex with exactly 3 placed neighbours: intersection = circle, exactly
    parametrized x(t) = c + r(cos t * u + sin t * v) with certified interval
    enclosures of the exact c, r, u, v; adaptive subdivision over t.
A branch dies when an interval constraint certifiably excludes 0 (or a
discriminant / radius-squared is certifiably negative).
If ALL branches die, G has NO realization in R^4 (certified).
Surviving sub-cells at max depth are reported as unresolved.
"""
import math, time, json, sys
from itertools import combinations
from ival import IV, isqrt, icos, isin, vadd, vsub, vscale, vdot, vnorm2

D = 4
TWO_PI = 2 * math.pi

def seed_coords(k):
    """unit-edge regular (k-1)-simplex, exact interval enclosures, k <= 5"""
    r3 = isqrt(IV(3.0)); r6 = isqrt(IV(6.0)); r10 = isqrt(IV(10.0))
    pts = [
        [IV(0.0), IV(0.0), IV(0.0), IV(0.0)],
        [IV(1.0), IV(0.0), IV(0.0), IV(0.0)],
        [IV(0.5), r3 / 2, IV(0.0), IV(0.0)],
        [IV(0.5), r3 / 6, r6 / 3, IV(0.0)],
        [IV(0.5), r3 / 6, r6 / 12, r10 / 4],
    ]
    return pts[:k]

def mig(x):
    if x.contains0(): return 0.0
    return min(abs(x.a), abs(x.b))

def gauss_solve(rows, nuns=D):
    """rows: list of (coeff IV list, rhs IV). Returns (xp, nulls) or None.
    xp: particular solution (free vars = 0); nulls: basis of null space."""
    m = len(rows)
    A = [list(r[0]) for r in rows]; b = [r[1] for r in rows]
    piv = []  # (row, col)
    usedr, usedc = set(), set()
    for _ in range(m):
        best, br, bc = 0.0, -1, -1
        for i in range(m):
            if i in usedr: continue
            for j in range(nuns):
                if j in usedc: continue
                g = mig(A[i][j])
                if g > best: best, br, bc = g, i, j
        if br < 0 or best == 0.0: return None
        usedr.add(br); usedc.add(bc); piv.append((br, bc))
        for i in range(m):
            if i == br or i in usedr - {br}: pass
            if i == br: continue
            f = A[i][bc] / A[br][bc]
            for j in range(nuns):
                A[i][j] = A[i][j] - f * A[br][j]
            b[i] = b[i] - f * b[br]
            A[i][bc] = IV(0.0)
    freec = [j for j in range(nuns) if j not in usedc]
    def backsolve(rhs_extra, free_assign):
        x = [IV(0.0)] * nuns
        for j, val in zip(freec, free_assign): x[j] = val
        for (i, j) in reversed(piv):
            s = rhs_extra[i]
            for jj in range(nuns):
                if jj != j and not (x[jj].a == 0.0 == x[jj].b and False):
                    s = s - A[i][jj] * x[jj]
            x[j] = s / A[i][j]
        return x
    xp = backsolve(b, [IV(0.0)] * len(freec))
    nulls = []
    for fi in range(len(freec)):
        assign = [IV(1.0) if t == fi else IV(0.0) for t in range(len(freec))]
        nv = backsolve([IV(0.0)] * m, assign)
        nulls.append(nv)
    return xp, nulls

def sphere_rows(nbr_pts):
    """|x-p0|^2 = |x-pi|^2  ->  2(pi-p0).x = |pi|^2 - |p0|^2"""
    p0 = nbr_pts[0]
    rows = []
    for p in nbr_pts[1:]:
        coef = [2 * (a - b) for a, b in zip(p, p0)]
        rhs = vnorm2(p) - vnorm2(p0)
        rows.append((coef, rhs))
    return rows, p0

def edge_check(x, p):
    e = vnorm2(vsub(x, p)) - 1
    return not (e.gt0() or e.lt0())   # True = consistent (contains 0)

import numpy as np

def _best_quad(x, nbr_pts):
    """pick 4 sphere centres with best midpoint Jacobian conditioning (float)"""
    xm = np.array([xi.mid() for xi in x])
    P = np.array([[p.mid() for p in pt] for pt in nbr_pts])
    k = len(nbr_pts)
    if k == 4: return (0, 1, 2, 3)
    best, bq = 0.0, (0, 1, 2, 3)
    for q in combinations(range(k), 4):
        J = 2 * (xm[None, :] - P[list(q)])
        dt = abs(np.linalg.det(J))
        if dt > best: best, bq = dt, q
        if best > 0.75: break
    return bq

def newton_refine(x, nbr_pts, iters=3):
    """Krawczyk contraction for F_i(y)=|y-p_i|^2-1 on 4 chosen spheres.
    Returns contracted enclosure or None if certifiably empty.
    K(X) = m - Y F(m) + (I - Y J(X))(X - m); solutions in X lie in K(X)."""
    q = _best_quad(x, nbr_pts)
    ps = [nbr_pts[i] for i in q]
    for _ in range(iters):
        mid = [xi.mid() for xi in x]
        m = [IV(v) for v in mid]
        # float approximate inverse of midpoint Jacobian
        Pm = np.array([[p.mid() for p in pt] for pt in ps])
        Jm = 2 * (np.array(mid)[None, :] - Pm)
        try:
            Y = np.linalg.inv(Jm)
        except np.linalg.LinAlgError:
            return x
        Fm = [vnorm2(vsub(m, p)) - 1 for p in ps]           # exact IVs
        # YF = Y * F(m) (interval)
        YF = []
        for i in range(4):
            s = IV(0.0)
            for j in range(4):
                s = s + IV(Y[i, j]) * Fm[j]
            YF.append(s)
        # M = I - Y*J(X) (interval), J(X) rows: 2(X - p_i)
        JX = [[2 * (x[c] - ps[r][c]) for c in range(4)] for r in range(4)]
        d = [xi - mi for xi, mi in zip(x, m)]               # X - m
        K = []
        for i in range(4):
            s = m[i] - YF[i]
            for j in range(4):
                mij = (IV(1.0) if i == j else IV(0.0))
                acc = mij
                for c in range(4):
                    acc = acc - IV(Y[i, c]) * JX[c][j]
                s = s + acc * d[j]
            K.append(s)
        newx, empty, shrunk = [], False, False
        for xi, ki in zip(x, K):
            lo = max(xi.a, ki.a); hi = min(xi.b, ki.b)
            if lo > hi: return None                        # certified empty
            if hi - lo < (xi.b - xi.a) - 1e-15: shrunk = True
            newx.append(IV(lo, hi))
        x = newx
        if not shrunk: break
    return x


def det4(rows):
    """interval determinant of 4x4 (list of 4 IV-vectors)"""
    def det2(a, b, c, d): return a * d - b * c
    def det3(m):
        return (m[0][0] * det2(m[1][1], m[1][2], m[2][1], m[2][2])
                - m[0][1] * det2(m[1][0], m[1][2], m[2][0], m[2][2])
                + m[0][2] * det2(m[1][0], m[1][1], m[2][0], m[2][1]))
    tot = None
    for j in range(4):
        cols = [k for k in range(4) if k != j]
        sub = [[rows[i][k] for k in cols] for i in (1, 2, 3)]
        term = rows[0][j] * det3(sub)
        if j % 2: term = -term
        tot = term if tot is None else tot + term
    return tot

def forced_coincident(adj, pts, u, v):
    """True if x_u = x_v is FORCED in every true solution with these boxes:
    if x_u != x_v, all common neighbours of u,v lie on the perpendicular
    bisector hyperplane; 5 common neighbours affinely spanning R^4
    (certified nonzero volume) contradict that."""
    common = [w for w in pts if w not in (u, v)
              and (adj[v] >> w) & 1 and (adj[u] >> w) & 1]
    if len(common) < 5: return False
    from itertools import combinations as _comb
    for W in _comb(common, 5):
        base = pts[W[0]]
        rows = [vsub(pts[w], base) for w in W[1:]]
        dt = det4(rows)
        if dt.gt0() or dt.lt0():
            return True
    return False

def boxes_overlap(x, y):
    return all(xi.a <= yi.b and yi.a <= xi.b for xi, yi in zip(x, y))

def place_sphere(nbr_pts):
    """>=4 neighbour points -> list of candidate x enclosures (0..2) or None
    if ill-conditioned. Uses first 4; caller passes others as checks."""
    rows, p0 = sphere_rows(nbr_pts[:4])
    sol = gauss_solve(rows)
    if sol is None: return None
    xp, nulls = sol
    if len(nulls) != 1: return None
    n = nulls[0]
    w = vsub(xp, p0)
    a = vnorm2(n); bq = 2 * vdot(n, w); c = vnorm2(w) - 1
    if not a.gt0(): return None
    disc = bq.sq() - 4 * a * c
    if disc.lt0(): return []            # certified: no intersection
    sq = isqrt(disc)
    if disc.contains0():
        # roots not separated: one merged enclosure instead of two fat copies
        t = ((-bq) + IV(-sq.b, sq.b)) / (2 * a)
        return [(vadd(xp, vscale(t, n)), (xp, n, t))]
    out = []
    for sgn in (1, -1):
        t = ((-bq) + (sq if sgn > 0 else -sq)) / (2 * a)
        out.append((vadd(xp, vscale(t, n)), (xp, n, t)))
    return out

def circle_setup(nbr_pts):
    """3 neighbour points -> (c, r, u, v) exact-enclosure circle params,
    or [] if certifiably empty, or None if ill-conditioned/degenerate."""
    p0, p1, p2 = nbr_pts
    d1, d2 = vsub(p1, p0), vsub(p2, p0)
    # circumcenter c = p0 + a d1 + b d2 : (c - (p0+pi)/2) . di = 0
    r11, r12, r22 = vdot(d1, d1), vdot(d1, d2), vdot(d2, d2)
    rows = [([r11, r12, IV(0.0), IV(0.0)], r11 / 2),
            ([r12, r22, IV(0.0), IV(0.0)], r22 / 2)]
    sol = gauss_solve(rows, nuns=2)
    if sol is None: return None
    (ab, _) = sol
    a, b = ab
    c = vadd(p0, vadd(vscale(a, d1), vscale(b, d2)))
    R2 = vnorm2(vsub(c, p0))
    r2 = 1 - R2
    if r2.lt0(): return []             # certified: spheres don't meet
    r = isqrt(r2)
    # orthonormal basis of complement of span(d1,d2)
    # orthonormalize the hull directions FIRST (classical Gram-Schmidt),
    # then complete to an orthonormal basis of the complement
    obasis = []
    for q in (d1, d2):
        e = list(q)
        for b in obasis:
            e = vsub(e, vscale(vdot(e, b), b))
        nn = vnorm2(e)
        if not (nn.gt0() and mig(nn) > 1e-8): return None
        obasis.append([x / isqrt(nn) for x in e])
    out = []
    for k in range(D):
        e = [IV(1.0) if i == k else IV(0.0) for i in range(D)]
        for b in obasis + out:
            e = vsub(e, vscale(vdot(e, b), b))
        nn = vnorm2(e)
        if nn.gt0() and mig(nn) > 1e-6:
            out.append([x / isqrt(nn) for x in e])
        if len(out) == 2: break
    if len(out) != 2: return None
    return (c, r, out[0], out[1])

class Prover:
    def __init__(self, adj, order, seed, max_nodes=6_000_000,
                 theta_min=TWO_PI / (1 << 16), max_width=1.0):
        self.adj = adj; self.order = order; self.seed = seed
        self.n = len(adj)
        self.max_nodes = max_nodes; self.theta_min = theta_min
        self.max_width = max_width
        self.nodes = 0
        self.unresolved = []
        self.aborted = False

    def run(self):
        pts = {v: p for v, p in zip(self.seed, seed_coords(len(self.seed)))}
        self.dfs(pts, 0)
        return len(self.unresolved) == 0 and not self.aborted

    def dfs(self, pts, oi, cellw=None):
        if self.aborted: return
        self.nodes += 1
        if self.nodes > self.max_nodes:
            self.aborted = True; return
        if oi == len(self.order):
            if self.final_sweep(pts):
                self.unresolved.append(dict(pts))   # full survivor
            return
        v = self.order[oi]
        nbrs = [u for u in pts if (self.adj[v] >> u) & 1]
        nbr_pts = [pts[u] for u in nbrs]
        if len(nbr_pts) >= 4:
            cands = None
            # try subsets of 4 for conditioning
            for quad in combinations(range(len(nbr_pts)), 4):
                cands = place_sphere([nbr_pts[i] for i in quad])
                if cands is not None: break
            if cands is None:
                # ill-conditioned 4-sphere system: if we are already inside a
                # fine circle cell, resolve v by a nested circle; else defer.
                if cellw is not None and cellw <= 64 * self.theta_min:
                    self.circle_stage(pts, oi, v, nbrs, sub=True, outerw=cellw)
                else:
                    self.unresolved.append(("illcond", v))
                return
            for x, seg in cands:
                x = newton_refine(x, nbr_pts)
                if x is None: continue          # certified empty
                if max(xi.width() for xi in x) > self.max_width:
                    if cellw is not None and cellw <= 8 * self.theta_min:
                        self.segment_stage(pts, oi, v, nbrs, seg, cellw)
                    else:
                        self.unresolved.append(("wide", v))
                    continue
                if self.noninjective(v, x, pts):
                    continue                     # branch forces x_v = x_u: no
                                                 # 13 distinct points here
                if all(edge_check(x, pts[u]) for u in nbrs):
                    pts2 = dict(pts); pts2[v] = x
                    if cellw is not None and not self.local_sweep(pts2, v):
                        continue                 # certified kill by propagation
                    self.dfs(pts2, oi + 1, cellw)
        elif len(nbr_pts) == 3:
            self.circle_stage(pts, oi, v, nbrs, outerw=cellw)
        else:
            raise RuntimeError(f"vertex {v} has only {len(nbr_pts)} placed nbrs")

    TIGHT = 1e-9

    def noninjective(self, v, x, pts):
        """True if this candidate branch certifiably forces x_v to coincide
        with an already-placed point u NOT adjacent to v (so the branch can
        contain no realization with all points distinct), or certifiably
        contains no true solution at all.

        Rule: pick 4 placed common neighbours Q of v and u. In any true
        solution both x_v and x_u solve the 4-sphere system over Q (all
        edges to Q are exactly unit). If the system's two root enclosures are
        disjoint and x-box and u-box overlap the SAME single root enclosure
        only, then x_v = x_u in every true solution of this branch."""
        for u in pts:
            if (self.adj[v] >> u) & 1: continue      # adjacent: edge_check rules
            if not boxes_overlap(x, pts[u]): continue
            pts_tmp = dict(pts); pts_tmp[v] = x
            if forced_coincident(self.adj, pts_tmp, u, v):
                return True                      # forced x_v = x_u
            common = [w for w in pts
                      if w != u and (self.adj[v] >> w) & 1 and (self.adj[u] >> w) & 1]
            if len(common) < 4: continue
            for Q in combinations(common, 4):
                rts = place_sphere([pts[w] for w in Q])
                if rts is None: continue
                if rts == []: return True            # v must solve Q-system: dead
                if len(rts) != 2: continue
                r1, r2 = rts[0][0], rts[1][0]
                if boxes_overlap(r1, r2): continue   # roots not separated
                d1, d2 = boxes_overlap(pts[u], r1), boxes_overlap(pts[u], r2)
                c1, c2 = boxes_overlap(x, r1), boxes_overlap(x, r2)
                if not d1 and not d2: return True    # u solves neither: dead
                if not c1 and not c2: return True    # v solves neither: dead
                if (d1 and not d2 and c1 and not c2) or \
                   (d2 and not d1 and c2 and not c1):
                    return True                      # forced coincidence
                break
        return False

    def local_sweep(self, pts, v0, rounds=2):
        """Contract v0's placed neighbourhood after placing v0 (Gauss-Seidel,
        local). False = certifiably empty."""
        seedset = set(self.seed)
        frontier = [v0]
        for _ in range(rounds):
            nxt = []
            for v in frontier:
                for u in list(pts):
                    if u in seedset or not (self.adj[v] >> u) & 1: continue
                    if max(c.width() for c in pts[u]) < self.TIGHT: continue
                    nbrs_u = [w for w in pts if (self.adj[u] >> w) & 1]
                    if len(nbrs_u) < 4: continue
                    x = newton_refine(pts[u], [pts[w] for w in nbrs_u], iters=1)
                    if x is None: return False
                    if any(a.width() < b.width() - 1e-14
                           for a, b in zip(x, pts[u])):
                        nxt.append(u)
                    pts[u] = x
                    for w in nbrs_u:
                        if not edge_check(pts[u], pts[w]): return False
            frontier = nxt
            if not frontier: break
        return True

    def final_sweep(self, pts, sweeps=6):
        """Gauss-Seidel interval-Newton contraction over all non-seed vertices
        using ALL placed neighbours. Returns False if certifiably empty."""
        pts = dict(pts)
        seedset = set(self.seed)
        for _ in range(sweeps):
            changed = False
            for v in pts:
                if v in seedset: continue
                if max(c.width() for c in pts[v]) < self.TIGHT: continue
                nbrs = [u for u in pts if (self.adj[v] >> u) & 1]
                if len(nbrs) < 4: continue
                x = newton_refine(pts[v], [pts[u] for u in nbrs], iters=2)
                if x is None: return False       # certified kill
                if any(a.width() < b.width() - 1e-15
                       for a, b in zip(x, pts[v])): changed = True
                pts[v] = x
                for u in nbrs:
                    if not edge_check(pts[v], pts[u]): return False
            if not changed: break
        return self.injective_or_unresolved(pts)

    def injective_or_unresolved(self, pts):
        """Leaf filter: False = cell certifiably contains no realization with
        all points distinct (forced coincidence, or emptiness). True = keep."""
        for u in pts:
            for v in pts:
                if u >= v or (self.adj[v] >> u) & 1: continue
                if not boxes_overlap(pts[u], pts[v]): continue
                if forced_coincident(self.adj, pts, u, v):
                    return False                 # forced x_u = x_v: degenerate
                common = [w for w in pts
                          if w not in (u, v)
                          and (self.adj[v] >> w) & 1 and (self.adj[u] >> w) & 1]
                if len(common) < 4: continue
                for Q in combinations(common, 4):
                    rts = place_sphere([pts[w] for w in Q])
                    if rts is None: continue
                    if rts == []: return False       # both must solve: empty
                    if len(rts) != 2: continue
                    r1, r2 = rts[0][0], rts[1][0]
                    if boxes_overlap(r1, r2): continue
                    du = (boxes_overlap(pts[u], r1), boxes_overlap(pts[u], r2))
                    dv = (boxes_overlap(pts[v], r1), boxes_overlap(pts[v], r2))
                    if not any(du) or not any(dv): return False   # empty
                    if du == dv and sum(du) == 1:
                        return False                 # forced x_u = x_v
                    break
        return True


    def segment_stage(self, pts, oi, v, nbrs, seg, cellw):
        """Wide sphere-root: subdivide the root segment x = xp + t n over the
        t-interval; each sub-cell is checked and recursed like a circle cell."""
        xp, n, tint = seg
        floor = max(1e-7, (cellw / 4.0) if cellw is not None else 1e-7)
        NIN = 16
        lo, hi = tint.a, tint.b
        if not (hi > lo): return
        stack = [(lo + i * (hi - lo) / NIN, lo + (i + 1) * (hi - lo) / NIN)
                 for i in range(NIN)]
        while stack:
            if self.aborted: return
            t0, t1 = stack.pop()
            tc = IV(t0, t1)
            x = vadd(xp, vscale(tc, n))
            if any(not edge_check(x, pts[u]) for u in nbrs):
                continue
            if self.noninjective(v, x, pts):
                continue
            before = len(self.unresolved)
            pts2 = dict(pts); pts2[v] = x
            if not self.local_sweep(pts2, v):
                continue
            self.dfs(pts2, oi + 1, cellw)
            if len(self.unresolved) > before:
                if t1 - t0 > floor:
                    del self.unresolved[before:]
                    tm = 0.5 * (t0 + t1)
                    stack.append((t0, tm)); stack.append((tm, t1))
                # else: keep unresolved records

    def circle_stage(self, pts, oi, v, nbrs, sub=False, outerw=None):
        nbr_pts = [pts[u] for u in nbrs]
        setup = None; tri = None
        for tri in combinations(range(len(nbr_pts)), 3):
            setup = circle_setup([nbr_pts[i] for i in tri])
            if setup is not None: break
        if setup is None:
            self.unresolved.append(("degenerate-circle", v)); return
        if setup == []: return          # certified kill
        tri_ids = {nbrs[i] for i in tri}
        c, r, u1, u2 = setup
        NINIT = 32
        floor = self.theta_min
        if outerw is not None:                    # nested: refine only somewhat
            floor = max(floor, outerw / 4.0)      # below the enclosing scale
        lo0, hi0 = (0.0, TWO_PI)
        if outerw is None and getattr(self, 'theta0_range', None) is not None \
                and not getattr(self, '_theta0_used', False):
            lo0, hi0 = self.theta0_range
            self._theta0_used = True
        stack = [(lo0 + i * (hi0 - lo0) / NINIT, lo0 + (i + 1) * (hi0 - lo0) / NINIT)
                 for i in range(NINIT)]
        cells = 0
        import time as _t
        tstart = _t.time()
        while stack:
            if self.aborted: return
            t0, t1 = stack.pop()
            cells += 1
            if cells % 500 == 0:
                print(f"      [circle v{v}] cells={cells} stack={len(stack)} "
                      f"width={t1-t0:.2e} nodes={self.nodes} "
                      f"unres={len(self.unresolved)} t={_t.time()-tstart:.0f}s",
                      flush=True)
            th = IV(t0, t1)
            x = vadd(c, vadd(vscale(r * icos(th), u1), vscale(r * isin(th), u2)))
            bad = False
            for u in nbrs:
                if not edge_check(x, pts[u]):
                    if u in tri_ids:
                        raise RuntimeError(
                            f"BUG: circle point violates defining sphere {u}")
                    bad = True; break
            if bad:
                continue                # cell dies on genuine extra constraints
            before = len(self.unresolved)
            pts2 = dict(pts); pts2[v] = x
            if not self.local_sweep(pts2, v):
                continue                    # cell certifiably dead
            self.dfs(pts2, oi + 1, t1 - t0)
            if len(self.unresolved) > before:
                # something survived; refine if allowed
                if t1 - t0 > floor:
                    del self.unresolved[before:]
                    tm = 0.5 * (t0 + t1)
                    stack.append((t0, tm)); stack.append((tm, t1))
                # else: keep unresolved records

def make_order(adj, n):
    """Dijkstra over subsets minimizing circle steps; returns (seed, order)."""
    import heapq
    # seeds: K5s then K4s
    def cliques(k):
        res = []
        def ext(clq, cand, start):
            if len(clq) == k: res.append(tuple(clq)); return
            for v in range(start, n):
                if (cand >> v) & 1: ext(clq + [v], cand & adj[v], v + 1)
        ext([], (1 << n) - 1, 0)
        return res
    seeds = cliques(5) or cliques(4)
    full = (1 << n) - 1
    best = {}; pq = []; pred = {}
    for sd in seeds:
        s0 = 0
        for v in sd: s0 |= 1 << v
        if best.get(s0, 99) > 0:
            best[s0] = 0; heapq.heappush(pq, (0, s0)); pred[s0] = (None, None, sd)
    target = None
    while pq:
        cst, P = heapq.heappop(pq)
        if cst > best.get(P, 99): continue
        if P == full: target = P; break
        for v in range(n):
            if (P >> v) & 1: continue
            k = bin(adj[v] & P).count('1')
            stepc = 0 if k >= 4 else (1 if k == 3 else None)
            if stepc is None: continue
            Q = P | (1 << v); nc = cst + stepc
            if nc < best.get(Q, 99):
                best[Q] = nc; pred[Q] = (P, v, None)
                heapq.heappush(pq, (nc, Q))
    if target is None: return None
    order = []
    P = full
    while pred[P][0] is not None:
        Pp, v, _ = pred[P]
        order.append(v); P = Pp
    seed = pred[P][2]
    order.reverse()
    return seed, order, best[full]

def make_order2(adj, n):
    """Layered circle-late order: maximize sphere placements before each circle.
    Returns (seed, order, ncircles) with order = vertex list."""
    def cliques(k):
        res = []
        def ext(clq, cand, start):
            if len(clq) == k: res.append(tuple(clq)); return
            for v in range(start, n):
                if (cand >> v) & 1: ext(clq + [v], cand & adj[v], v + 1)
        ext([], (1 << n) - 1, 0)
        return res
    best = None
    for seed in (cliques(5) or cliques(4)):
        placed = set(seed); order = []; ncirc = 0; ok = True
        while len(placed) < n:
            # close under >=4-neighbour placements, most-constrained first
            while True:
                cands = [(sum(1 for u in placed if (adj[v] >> u) & 1), v)
                         for v in range(n) if v not in placed]
                cands = [(k, v) for k, v in cands if k >= 4]
                if not cands: break
                k, v = max(cands)
                order.append(v); placed.add(v)
            if len(placed) == n: break
            # pick a 3-neighbour circle vertex maximizing subsequent closure
            c3 = [v for v in range(n) if v not in placed
                  and sum(1 for u in placed if (adj[v] >> u) & 1) == 3]
            if not c3: ok = False; break
            def closure_size(w):
                P = set(placed); P.add(w)
                grew = True
                while grew:
                    grew = False
                    for x in range(n):
                        if x not in P and sum(1 for u in P if (adj[x] >> u) & 1) >= 4:
                            P.add(x); grew = True
                return len(P)
            v = max(c3, key=closure_size)
            order.append(v); placed.add(v); ncirc += 1
        if not ok: continue
        # score: fewer circles, then longer prefix before first circle
        # prefix length = index of first circle vertex in order
        circ_positions = []
        placed2 = set(seed)
        for i, v in enumerate(order):
            k = sum(1 for u in placed2 if (adj[v] >> u) & 1)
            if k == 3: circ_positions.append(i)
            placed2.add(v)
        key = (ncirc, [-p for p in circ_positions])
        if best is None or key < best[0]:
            best = (key, seed, order, ncirc)
    if best is None: return None
    return best[1], best[2], best[3]

def decide_graph(adj, n, verbose=True, **kw):
    mo = make_order2(adj, n) or make_order(adj, n)
    if mo is None: return "no-order", None
    seed, order, ncirc = mo
    pr = Prover(adj, order, seed, **kw)
    t0 = time.time()
    killed = pr.run()
    dt = time.time() - t0
    status = "KILLED" if killed else ("ABORT" if pr.aborted else "SURVIVORS")
    if verbose:
        print(f"    seed={seed} circles={ncirc} nodes={pr.nodes} "
              f"unresolved={len(pr.unresolved)} t={dt:.1f}s -> {status}", flush=True)
    return status, pr

if __name__ == "__main__":
    fn = sys.argv[1]; idxs = [int(x) for x in sys.argv[2:]] if len(sys.argv) > 2 else None
    data = json.load(open(fn))
    graphs = data["graphs"]; n = data["n"]
    minimal = set(data.get("minimal", []))
    todo = idxs if idxs is not None else range(len(graphs))
    summary = {}
    for gi in todo:
        print(f"[{gi}] minimal={gi in minimal}", flush=True)
        st, pr = decide_graph(graphs[gi], n)
        summary[gi] = st
    print(json.dumps(summary))
