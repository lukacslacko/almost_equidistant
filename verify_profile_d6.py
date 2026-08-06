#!/usr/bin/env python3
"""Independent controls and report checks for exact d=6 graph filters.

This intentionally does not import or call profile_d6.c.  It reconstructs the
18-point Larman--Rogers-plus-two-apices lower-bound graph and checks every link
bound in plain Python.  It also checks synthetic witnesses for the K5-link and
K7-reflection rules.
"""
import hashlib
import json
from itertools import combinations
from pathlib import Path


BOUNDS = {2: 16, 3: 12, 4: 10, 5: 4, 6: 2, 7: 0}
ROOT = Path(__file__).resolve().parent


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


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_committed_report():
    """Check fixed totals, population partitions, and manifest hashes."""

    report_path = ROOT / "d6_profile.json"
    manifest_path = ROOT / "d6_profile_manifest.json"
    with report_path.open(encoding="utf-8") as stream:
        report = json.load(stream)
    with manifest_path.open(encoding="utf-8") as stream:
        manifest = json.load(stream)

    assert report["schema"] == 2
    assert report["n"] == 19 and report["dimension"] == 6
    assert report["candidate_count_expected"] == 3_971_787
    assert report["candidate_count_profiled"] == 3_971_787
    assert report["kill_log_unique_indices"] == 3_055_474

    populations = report["populations"]
    fixed = {
        "all": (3_971_787, 1_017_404, 3_680_379, 73, 3_680_381, 291_406),
        "certified": (3_055_474, 1_004_938, 3_053_023, 69, 3_053_023, 2_451),
        "deferred": (916_313, 12_466, 627_356, 4, 627_358, 288_955),
    }
    for name, (total, previous, cover, tight, exact, residue) in fixed.items():
        item = populations[name]
        assert item["total"] == total
        assert item["cumulative_previous_exact_rejections"] == previous
        assert item["K7_clique_Hall_rejections"] == 0
        assert item["K6_clique_Hall_rejections"] == 0
        assert item["K7_disjoint_edge_bounded_cover_rejections"] == cover
        assert item["K7_tight_cover_matching_rejections"] == tight
        assert item["cumulative_exact_rejections"] == exact
        assert sum(item["clique_number"].values()) == total
        assert sum(item["exact_rejections_by_clique_number"].values()) == exact
        assert sum(item["residue_by_clique_number"].values()) == residue
        assert exact + residue == total
        for witnesses in item["new_rule_witnesses"].values():
            assert len(witnesses) <= 5
            assert [w["index"] for w in witnesses] == sorted(
                w["index"] for w in witnesses
            )

    scalar_fields = (
        "total",
        "cumulative_link_rejections",
        "K7_reflection_rejections",
        "K7_defect_CSP_rejections",
        "K7_clique_Hall_rejections",
        "K7_disjoint_edge_bounded_cover_rejections",
        "K7_tight_cover_matching_rejections",
        "K6_clique_Hall_rejections",
        "cumulative_previous_exact_rejections",
        "cumulative_exact_rejections",
    )
    for field in scalar_fields:
        assert populations["all"][field] == (
            populations["certified"][field]
            + populations["deferred"][field]
        )
    for field in (
        "clique_number",
        "link_rejections",
        "exact_rejections_by_clique_number",
        "residue_by_clique_number",
    ):
        for key, value in populations["all"][field].items():
            assert value == (
                populations["certified"][field][key]
                + populations["deferred"][field][key]
            )

    profiler = manifest["profiler"]
    assert profiler["source_sha256"] == sha256(ROOT / "profile_d6.c")
    assert profiler["result_sha256"] == sha256(report_path)
    reference = manifest["reference_crosscheck"]
    assert reference["sample_sha256"] == sha256(
        ROOT / "d6_reference_sample.json"
    )


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
    verify_committed_report()
    print("profile_d6 controls and committed report: PASS")


if __name__ == "__main__":
    main()
