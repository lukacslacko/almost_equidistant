#!/usr/bin/env python3
"""Plain-Python reference checks for the new dimension-six graph filters.

This file is intentionally independent of ``profile_d6.c``.  It favors a
small, transparent implementation over throughput and is intended for
controls and a committed cross-check sample, not the 3,971,787-graph run.

The graph representation is a list of Python integer adjacency bit masks.
Candidate nonedges are used only as *allowed* defect coordinates; they are
never required to have non-unit distance.

One simplifying consequence is made explicit in every Hall witness.  For K7,
an outside clique of size ``|A|+1`` with all defects in ``A`` joins the
``7-|A|`` unaffected seed vertices to form a required K8.  For K6 the sizes
are ``|A|+2`` and ``6-|A|``, again giving a K8.  Thus both requested Hall
filters are independently useful algebra controls but add no coverage to this
already K8-free corpus.

The tight-cover equality rule is kept separate from the old bound.  It is
consulted only for a seed whose eligible cover number is exactly seven, then
enumerates every eligible seven-cover and independently searches for a perfect
matching from its seven allowed defect masks onto the seven seed coordinates.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable, Iterator, Sequence


@dataclass
class FilterDecision:
    """One graph-level filter decision, including a checkable witness."""

    name: str
    applicable: bool
    rejected: bool
    seeds_checked: int
    witness: dict | None = None
    note: str | None = None

    def jsonable(self) -> dict:
        return asdict(self)


def validate_graph(adj: Sequence[int]) -> None:
    """Raise ``ValueError`` unless ``adj`` is a simple undirected graph."""

    n = len(adj)
    full = (1 << n) - 1
    for u, row in enumerate(adj):
        if row & ~full:
            raise ValueError(f"adjacency row {u} has a bit outside the graph")
        if row & (1 << u):
            raise ValueError(f"loop at vertex {u}")
        for v in range(u):
            if bool(row & (1 << v)) != bool(adj[v] & (1 << u)):
                raise ValueError(f"asymmetric adjacency at {u},{v}")


def graph_from_edges(n: int, edges: Iterable[tuple[int, int]]) -> list[int]:
    """Build adjacency bit masks from an edge iterable."""

    adj = [0] * n
    for u, v in edges:
        if not (0 <= u < n and 0 <= v < n) or u == v:
            raise ValueError(f"bad edge {(u, v)} for n={n}")
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    validate_graph(adj)
    return adj


def add_edge(adj: list[int], u: int, v: int) -> None:
    """Add one undirected edge to an adjacency-mask graph."""

    if u == v:
        raise ValueError("a simple graph cannot have a loop")
    adj[u] |= 1 << v
    adj[v] |= 1 << u


def add_clique(adj: list[int], vertices: Iterable[int]) -> None:
    """Add every edge on ``vertices``."""

    for u, v in combinations(tuple(vertices), 2):
        add_edge(adj, u, v)


def vertices(mask: int) -> list[int]:
    """Return the increasing list of set-bit positions in ``mask``."""

    out = []
    while mask:
        bit = mask & -mask
        out.append(bit.bit_length() - 1)
        mask ^= bit
    return out


def is_clique(adj: Sequence[int], vs: Sequence[int]) -> bool:
    """Return whether ``vs`` is a required-edge clique."""

    return all(adj[u] & (1 << v) for u, v in combinations(vs, 2))


def clique_masks(
    adj: Sequence[int], size: int, candidate_mask: int | None = None
) -> Iterator[int]:
    """Enumerate each clique of the requested size exactly once."""

    if candidate_mask is None:
        candidate_mask = (1 << len(adj)) - 1

    def visit(candidates: int, need: int, chosen: int) -> Iterator[int]:
        if need == 0:
            yield chosen
            return
        while candidates.bit_count() >= need:
            bit = candidates & -candidates
            candidates ^= bit
            v = bit.bit_length() - 1
            yield from visit(candidates & adj[v], need - 1, chosen | bit)

    yield from visit(candidate_mask, size, 0)


def find_clique_mask(adj: Sequence[int], candidate_mask: int, size: int) -> int:
    """Return one clique mask of ``size``, or zero if none exists."""

    if size == 0:
        return 0

    def visit(candidates: int, need: int, chosen: int) -> int:
        if need == 0:
            return chosen
        while candidates.bit_count() >= need:
            bit = candidates & -candidates
            candidates ^= bit
            v = bit.bit_length() - 1
            found = visit(candidates & adj[v], need - 1, chosen | bit)
            if found:
                return found
        return 0

    return visit(candidate_mask, size, 0)


def _seed_defects(
    adj: Sequence[int], seed_mask: int
) -> tuple[list[int], list[int], list[int]]:
    """Return seed vertices, outside vertices, and relative defect masks."""

    seed = vertices(seed_mask)
    seed_position = {v: i for i, v in enumerate(seed)}
    outside = [v for v in range(len(adj)) if not seed_mask & (1 << v)]
    defects = []
    for x in outside:
        absolute = seed_mask & ~adj[x]
        relative = 0
        for q in vertices(absolute):
            relative |= 1 << seed_position[q]
        defects.append(relative)
    return seed, outside, defects


def _hall_seed(
    adj: Sequence[int], seed_mask: int, universal_coordinates: int
) -> dict | None:
    """Check one seed by the coordinate-mask form of Hall's theorem.

    ``universal_coordinates`` is zero for a K7 seed and one for the virtual
    normal coordinate of a K6 seed.  The implementation deliberately scans
    every coordinate mask, independently of the equivalent K8 shortcut.
    """

    seed, outside, defects = _seed_defects(adj, seed_mask)
    coordinate_count = len(seed)
    for coordinate_mask in range(1 << coordinate_count):
        induced = 0
        for x, allowed in zip(outside, defects):
            if allowed & ~coordinate_mask == 0:
                induced |= 1 << x
        limit = coordinate_mask.bit_count() + universal_coordinates
        bad_clique = find_clique_mask(adj, induced, limit + 1)
        if not bad_clique:
            continue

        clique = vertices(bad_clique)
        selected_coordinates = [
            seed[i]
            for i in range(coordinate_count)
            if coordinate_mask & (1 << i)
        ]
        unaffected_seed = [q for q in seed if q not in selected_coordinates]
        forced_k8 = unaffected_seed + clique
        # This is a useful independent sanity check: each Hall failure is
        # literally a required K8, not a new obstruction on a K8-free corpus.
        assert len(forced_k8) == 8
        assert is_clique(adj, forced_k8)
        return {
            "seed": seed,
            "coordinate_subset": selected_coordinates,
            "outside_clique": clique,
            "allowed_defects": {
                str(x): [
                    seed[i]
                    for i in range(coordinate_count)
                    if allowed & (1 << i)
                ]
                for x, allowed in zip(outside, defects)
                if x in clique
            },
            "hall_limit": limit,
            "forced_K8": forced_k8,
        }
    return None


def k7_clique_hall_seed(adj: Sequence[int], seed: Sequence[int]) -> dict | None:
    """Return a K7 Hall-failure witness for one required K7, if present."""

    if len(seed) != 7 or not is_clique(adj, seed):
        raise ValueError("seed is not a required K7")
    return _hall_seed(adj, sum(1 << v for v in seed), 0)


def k6_clique_hall_seed(adj: Sequence[int], seed: Sequence[int]) -> dict | None:
    """Return a K6 Hall-failure witness for one required K6, if present."""

    if len(seed) != 6 or not is_clique(adj, seed):
        raise ValueError("seed is not a required K6")
    return _hall_seed(adj, sum(1 << v for v in seed), 1)


def k7_clique_hall(adj: Sequence[int]) -> FilterDecision:
    """Check every required K7 seed for the simplex-clique Hall rule."""

    checked = 0
    for seed_mask in clique_masks(adj, 7):
        checked += 1
        witness = _hall_seed(adj, seed_mask, 0)
        if witness is not None:
            return FilterDecision(
                "K7_clique_Hall", True, True, checked, witness
            )
    return FilterDecision(
        "K7_clique_Hall",
        checked > 0,
        False,
        checked,
        note=None if checked else "graph has no required K7 seed",
    )


def _minimum_eligible_cover(
    edge_masks: Sequence[int], eligible_mask: int, maximum_size: int | None = None
) -> int | None:
    """Return the exact minimum eligible vertex-cover mask, if one exists."""

    eligible = vertices(eligible_mask)
    last_size = len(eligible) if maximum_size is None else min(
        len(eligible), maximum_size
    )
    for size in range(last_size + 1):
        for chosen in combinations(eligible, size):
            cover = sum(1 << i for i in chosen)
            if all(edge & cover for edge in edge_masks):
                return cover
    return None


def _k7_cover_instance(
    adj: Sequence[int], seed: Sequence[int]
) -> tuple[list[int], list[int], list[int], list[int], list[tuple[int, int]], int]:
    """Construct L and its eligibility mask for one required K7 seed."""

    if len(seed) != 7 or not is_clique(adj, seed):
        raise ValueError("seed is not a required K7")
    seed_mask = sum(1 << v for v in seed)
    seed_vertices, outside, defects = _seed_defects(adj, seed_mask)
    edge_masks = []
    edge_vertices = []
    for i, j in combinations(range(len(outside)), 2):
        x, y = outside[i], outside[j]
        if adj[x] & (1 << y) and defects[i] & defects[j] == 0:
            edge_masks.append((1 << i) | (1 << j))
            edge_vertices.append((x, y))
    eligible_mask = sum(
        1 << i for i, allowed in enumerate(defects) if allowed.bit_count() >= 3
    )
    return (
        seed_vertices,
        outside,
        defects,
        edge_masks,
        edge_vertices,
        eligible_mask,
    )


def k7_bounded_cover_seed(
    adj: Sequence[int], seed: Sequence[int], cover_bound: int = 7
) -> dict | None:
    """Check the disjoint-defect edge cover rule for one required K7."""

    (
        seed_vertices,
        outside,
        defects,
        edge_masks,
        edge_vertices,
        eligible_mask,
    ) = _k7_cover_instance(adj, seed)
    minimum = _minimum_eligible_cover(
        edge_masks, eligible_mask, maximum_size=cover_bound
    )
    if minimum is not None:
        return None

    # This unrestricted search is used only to make a small diagnostic
    # witness.  It does not participate in the decision above.
    unrestricted_minimum = _minimum_eligible_cover(edge_masks, eligible_mask)

    ineligible_edge = next(
        (
            edge_vertices[i]
            for i, edge in enumerate(edge_masks)
            if edge & eligible_mask == 0
        ),
        None,
    )
    failure_kind = (
        "ineligible_disjoint_edge" if ineligible_edge is not None
        else "cover_cardinality"
    )
    return {
        "seed": seed_vertices,
        "failure_kind": failure_kind,
        "cover_bound": cover_bound,
        "minimum_eligible_cover_size": (
            None
            if unrestricted_minimum is None
            else unrestricted_minimum.bit_count()
        ),
        "minimum_eligible_cover": (
            None
            if unrestricted_minimum is None
            else [outside[i] for i in vertices(unrestricted_minimum)]
        ),
        "eligible_vertices": [outside[i] for i in vertices(eligible_mask)],
        "ineligible_edge": (
            None if ineligible_edge is None else list(ineligible_edge)
        ),
        "L_edges": [list(edge) for edge in edge_vertices],
        "allowed_defects": {
            str(x): [
                seed_vertices[i]
                for i in range(7)
                if allowed & (1 << i)
            ]
            for x, allowed in zip(outside, defects)
        },
    }


def _defect_perfect_matching(
    cover_mask: int, defects: Sequence[int]
) -> dict[int, int] | None:
    """Match seven selected local vertices to all seven seed coordinates."""

    selected = vertices(cover_mask)
    if len(selected) != 7:
        raise ValueError("the tight-cover matching requires seven vertices")
    # Smallest domains first is only an ordering optimization.  The recursion
    # still enumerates every possible injection using exact bit operations.
    selected.sort(key=lambda i: (defects[i].bit_count(), i))
    assignment: dict[int, int] = {}

    def visit(at: int, used: int) -> bool:
        if at == len(selected):
            return used == 0x7F
        local_vertex = selected[at]
        choices = defects[local_vertex] & ~used & 0x7F
        while choices:
            coordinate = choices & -choices
            choices ^= coordinate
            assignment[local_vertex] = coordinate.bit_length() - 1
            if visit(at + 1, used | coordinate):
                return True
        assignment.pop(local_vertex, None)
        return False

    return assignment.copy() if visit(0, 0) else None


def k7_tight_cover_matching_seed(
    adj: Sequence[int], seed: Sequence[int]
) -> dict | None:
    """Reject a tight size-seven cover equality case with no support matching.

    The rule applies to this seed only when L has an eligible cover of size
    seven but none of size at most six.  Every eligible size-seven cover is
    then checked.  The seed passes as soon as one cover can be perfectly
    matched to all seven allowed simplex coordinates; rejection occurs only
    after exhaustive failure of every such cover.
    """

    (
        seed_vertices,
        outside,
        defects,
        edge_masks,
        edge_vertices,
        eligible_mask,
    ) = _k7_cover_instance(adj, seed)
    if (
        _minimum_eligible_cover(edge_masks, eligible_mask, maximum_size=6)
        is not None
    ):
        return None
    if (
        _minimum_eligible_cover(edge_masks, eligible_mask, maximum_size=7)
        is None
    ):
        # The old bounded-cover rule already rejects this seed.
        return None

    eligible = vertices(eligible_mask)
    covers_checked = 0
    failed_cover_examples = []
    for chosen in combinations(eligible, 7):
        cover = sum(1 << i for i in chosen)
        if not all(edge & cover for edge in edge_masks):
            continue
        covers_checked += 1
        matching = _defect_perfect_matching(cover, defects)
        if matching is not None:
            return None
        if len(failed_cover_examples) < 8:
            coordinate_union = 0
            for i in chosen:
                coordinate_union |= defects[i]
            failed_cover_examples.append(
                {
                    "cover": [outside[i] for i in chosen],
                    "coordinate_union": [
                        seed_vertices[q]
                        for q in vertices(coordinate_union)
                    ],
                    "allowed_defects": {
                        str(outside[i]): [
                            seed_vertices[q] for q in vertices(defects[i])
                        ]
                        for i in chosen
                    },
                }
            )
    assert covers_checked > 0
    return {
        "seed": seed_vertices,
        "failure_kind": "tight_cover_no_perfect_matching",
        "size7_covers_checked": covers_checked,
        "failed_cover_examples": failed_cover_examples,
        "eligible_vertices": [outside[i] for i in eligible],
        "L_edges": [list(edge) for edge in edge_vertices],
    }


def k7_disjoint_edge_bounded_cover(adj: Sequence[int]) -> FilterDecision:
    """Check every K7 seed for an eligible vertex cover of size at most 7."""

    checked = 0
    for seed_mask in clique_masks(adj, 7):
        checked += 1
        seed = vertices(seed_mask)
        witness = k7_bounded_cover_seed(adj, seed)
        if witness is not None:
            return FilterDecision(
                "K7_disjoint_edge_bounded_cover",
                True,
                True,
                checked,
                witness,
            )
    return FilterDecision(
        "K7_disjoint_edge_bounded_cover",
        checked > 0,
        False,
        checked,
        note=None if checked else "graph has no required K7 seed",
    )


def k7_tight_cover_matching(adj: Sequence[int]) -> FilterDecision:
    """Check every K7 seed's tight size-seven cover equality case."""

    checked = 0
    for seed_mask in clique_masks(adj, 7):
        checked += 1
        seed = vertices(seed_mask)
        witness = k7_tight_cover_matching_seed(adj, seed)
        if witness is not None:
            return FilterDecision(
                "K7_tight_cover_matching", True, True, checked, witness
            )
    return FilterDecision(
        "K7_tight_cover_matching",
        checked > 0,
        False,
        checked,
        note=None if checked else "graph has no required K7 seed",
    )


def k6_clique_hall(adj: Sequence[int]) -> FilterDecision:
    """Check K6 Hall on the K6-only population, as specified in steering."""

    k7 = find_clique_mask(adj, (1 << len(adj)) - 1, 7)
    if k7:
        return FilterDecision(
            "K6_clique_Hall",
            False,
            False,
            0,
            note="contains a required K7; K6 rule is restricted to K6-only",
        )

    checked = 0
    for seed_mask in clique_masks(adj, 6):
        checked += 1
        witness = _hall_seed(adj, seed_mask, 1)
        if witness is not None:
            return FilterDecision(
                "K6_clique_Hall", True, True, checked, witness
            )
    return FilterDecision(
        "K6_clique_Hall",
        checked > 0,
        False,
        checked,
        note=None if checked else "graph has no required K6 seed",
    )


def evaluate(adj: Sequence[int]) -> dict[str, FilterDecision]:
    """Run the four requested graph-level reference filters."""

    validate_graph(adj)
    return {
        "K7_clique_Hall": k7_clique_hall(adj),
        "K7_disjoint_edge_bounded_cover": k7_disjoint_edge_bounded_cover(adj),
        "K7_tight_cover_matching": k7_tight_cover_matching(adj),
        "K6_clique_Hall": k6_clique_hall(adj),
    }


def evaluate_sample(sample: dict, include_graph_results: bool = False) -> dict:
    """Evaluate a ``d6_reference_sample.json``-format object."""

    names = (
        "K7_clique_Hall",
        "K7_disjoint_edge_bounded_cover",
        "K7_tight_cover_matching",
        "K6_clique_Hall",
    )
    strata: dict[str, dict] = {}
    graph_results = []
    for graph in sample["graphs"]:
        decisions = evaluate(graph["adjacency"])
        stratum = graph["stratum"]
        summary = strata.setdefault(
            stratum,
            {
                "graphs": 0,
                "applicable": {name: 0 for name in names},
                "rejected": {name: 0 for name in names},
            },
        )
        summary["graphs"] += 1
        for name, decision in decisions.items():
            summary["applicable"][name] += int(decision.applicable)
            summary["rejected"][name] += int(decision.rejected)
        if include_graph_results:
            graph_results.append(
                {
                    "index": graph["index"],
                    "stratum": stratum,
                    "decisions": {
                        name: decision.jsonable()
                        for name, decision in decisions.items()
                    },
                }
            )
    result = {
        "schema": 1,
        "sample_schema": sample["schema"],
        "sample_graphs": len(sample["graphs"]),
        "strata": strata,
    }
    if include_graph_results:
        result["graph_results"] = graph_results
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample", type=Path)
    parser.add_argument(
        "--details", action="store_true", help="include per-graph witnesses"
    )
    args = parser.parse_args()
    with args.sample.open(encoding="utf-8") as stream:
        sample = json.load(stream)
    print(json.dumps(evaluate_sample(sample, args.details), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
