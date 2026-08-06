#!/usr/bin/env python3
"""Independent exact reference for the post-cover K7 support/rank filters.

This deliberately does not import the production profiler or the earlier
Python reference.  It works with one required K7 at a time, enumerates every
possible zero-factor set (including nonminimal vertex covers), and applies
only exact bit-mask arguments.  The zero-factor vectors are also orthogonal
to every nonzero-factor vector, sharpening the K/B rank upper bounds by the
cover size.

Original candidate nonedges remain optional unit distances.  The exact zero
patterns used below arise only after fixing ``N = {x : c_x != 0}``: on N the
required-edge graph is exactly the intersection graph of the allowed defect
masks.  Consequently the normalized Gram matrix K has off-diagonal graph G[N]
and B=K-J has off-diagonal graph complement(G[N]).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass, field
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Iterable, Iterator, Sequence


def bits(mask: int) -> Iterator[int]:
    """Yield set-bit positions in increasing order."""

    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def validate_graph(adj: Sequence[int], require_alpha_two: bool = True) -> None:
    """Validate the graph and, by default, the almost-equidistant premise."""

    n = len(adj)
    full = (1 << n) - 1
    for u, row in enumerate(adj):
        if row & ~full:
            raise ValueError(f"row {u} has bits outside the graph")
        if row & (1 << u):
            raise ValueError(f"loop at {u}")
        for v in range(u):
            if bool(row & (1 << v)) != bool(adj[v] & (1 << u)):
                raise ValueError(f"asymmetric pair {u},{v}")
    if require_alpha_two:
        for u in range(n):
            non = full & ~adj[u] & ~(1 << u)
            todo = non
            while todo:
                bit = todo & -todo
                todo ^= bit
                v = bit.bit_length() - 1
                if todo & ~adj[v]:
                    raise ValueError("graph has an independent triple")


def clique_masks(adj: Sequence[int], size: int) -> Iterator[int]:
    """Enumerate required cliques exactly once."""

    def visit(candidates: int, need: int, chosen: int) -> Iterator[int]:
        if need == 0:
            yield chosen
            return
        while candidates.bit_count() >= need:
            bit = candidates & -candidates
            candidates ^= bit
            v = bit.bit_length() - 1
            yield from visit(candidates & adj[v], need - 1, chosen | bit)

    yield from visit((1 << len(adj)) - 1, size, 0)


def seed_instance(
    adj: Sequence[int], seed_mask: int
) -> tuple[list[int], list[int], tuple[int, ...], tuple[int, ...], int]:
    """Return seed/outside labels, defect masks, L adjacency, and eligibility."""

    seed = list(bits(seed_mask))
    if len(seed) != 7:
        raise ValueError("seed must have seven vertices")
    if any(not (adj[u] & (1 << v)) for u, v in combinations(seed, 2)):
        raise ValueError("seed is not a required K7")
    position = {q: i for i, q in enumerate(seed)}
    outside = [v for v in range(len(adj)) if not (seed_mask & (1 << v))]
    defects: list[int] = []
    for x in outside:
        allowed = 0
        for q in seed:
            if not (adj[x] & (1 << q)):
                allowed |= 1 << position[q]
        defects.append(allowed)

    ladj = [0] * len(outside)
    for i, j in combinations(range(len(outside)), 2):
        x, y = outside[i], outside[j]
        if adj[x] & (1 << y) and not (defects[i] & defects[j]):
            ladj[i] |= 1 << j
            ladj[j] |= 1 << i
    eligible = sum(
        1 << i for i, allowed in enumerate(defects)
        if allowed.bit_count() >= 3
    )
    return seed, outside, tuple(defects), tuple(ladj), eligible


def is_vertex_cover(ladj: Sequence[int], zmask: int) -> bool:
    """Return whether zmask meets every edge of L."""

    remaining = ((1 << len(ladj)) - 1) & ~zmask
    return all(not (ladj[v] & remaining) for v in bits(remaining))


def eligible_covers(ladj: Sequence[int], eligible: int, cap: int = 7) -> list[int]:
    """Enumerate *all* eligible covers, not merely minimum/minimal ones."""

    return [
        zmask
        for zmask in range(1 << len(ladj))
        if not (zmask & ~eligible)
        and zmask.bit_count() <= cap
        and is_vertex_cover(ladj, zmask)
    ]


SUPPORT_DOMAINS: tuple[tuple[int, ...], ...] = tuple(
    tuple(
        support
        for support in range(128)
        if not (support & ~allowed) and support.bit_count() >= 3
    )
    for allowed in range(128)
)


def matching_size(supports: Sequence[int]) -> int:
    """Maximum matching from columns with the given coordinate masks."""

    reachable = {0}
    for support in supports:
        following = set(reachable)  # This column may remain unmatched.
        for used in reachable:
            available = support & ~used & 0x7F
            while available:
                bit = available & -available
                available ^= bit
                following.add(used | bit)
        reachable = following
    return max((mask.bit_count() for mask in reachable), default=0)


def support_family_valid(supports: Sequence[int]) -> bool:
    """Check all necessary orthogonal zero-pattern rules on one family."""

    if any(s.bit_count() < 3 for s in supports):
        return False
    if any((s & t).bit_count() == 1 for s, t in combinations(supports, 2)):
        return False
    if matching_size(supports) != len(supports):
        return False
    if len(supports) != 7:
        return True
    rows = [
        sum(1 << column for column, support in enumerate(supports)
            if support & (1 << coordinate))
        for coordinate in range(7)
    ]
    if any(row.bit_count() in (0, 2) for row in rows):
        return False
    return not any(
        (first & second).bit_count() == 1
        for first, second in combinations(rows, 2)
    )


class SupportSolver:
    """Memoized exact support-CSP solver with deterministic witnesses."""

    def __init__(self) -> None:
        self.cache: dict[tuple[int, ...], tuple[int, ...] | None] = {}
        self.cache_hits = 0
        self.nodes = 0

    @staticmethod
    def _row_partial_possible(
        selected: Sequence[int], remaining_domains: Sequence[Sequence[int]]
    ) -> bool:
        """Safe forward checks for the extra square-orthogonal row rules."""

        if len(selected) + len(remaining_domains) != 7:
            return True
        for coordinate in range(7):
            current = sum(bool(s & (1 << coordinate)) for s in selected)
            possible = sum(
                any(s & (1 << coordinate) for s in domain)
                for domain in remaining_domains
            )
            if not any(
                total not in (0, 2)
                for total in range(current, current + possible + 1)
            ):
                return False
        for first, second in combinations(range(7), 2):
            both = (1 << first) | (1 << second)
            current = sum((s & both) == both for s in selected)
            if current == 1 and not any(
                any((s & both) == both for s in domain)
                for domain in remaining_domains
            ):
                return False
        return True

    def solve(self, allowed_masks: Iterable[int]) -> tuple[int, ...] | None:
        """Return one feasible actual-support family, or None if exhaustive."""

        allowed = tuple(sorted(allowed_masks))
        cached = self.cache.get(allowed, ...)
        if cached is not ...:
            self.cache_hits += 1
            return cached
        count = len(allowed)
        if count > 7:
            self.cache[allowed] = None
            return None
        domains = [list(SUPPORT_DOMAINS[d]) for d in allowed]
        if any(not domain for domain in domains):
            self.cache[allowed] = None
            return None
        assignment = [0] * count

        def visit(unassigned: tuple[int, ...], selected: tuple[int, ...]) -> bool:
            self.nodes += 1
            if not unassigned:
                return support_family_valid(assignment)

            candidates_by_variable: list[tuple[int, list[int]]] = []
            for variable in unassigned:
                # Equal allowed masks are interchangeable.  Always assign the
                # first remaining member of an equal-domain block.
                if any(
                    other < variable and allowed[other] == allowed[variable]
                    for other in unassigned
                ):
                    continue
                lower = 0
                if variable and allowed[variable - 1] == allowed[variable]:
                    lower = assignment[variable - 1]
                candidates = [
                    support for support in domains[variable]
                    if support >= lower
                    and all((support & previous).bit_count() != 1
                            for previous in selected)
                ]
                if not candidates:
                    return False
                candidates_by_variable.append((variable, candidates))
            variable, candidates = min(
                candidates_by_variable, key=lambda item: (len(item[1]), item[0])
            )
            following = tuple(v for v in unassigned if v != variable)
            for support in candidates:
                assignment[variable] = support
                chosen = selected + (support,)
                if matching_size(chosen) != len(chosen):
                    continue
                remaining_domains = []
                forward_ok = True
                for other in following:
                    domain = [
                        candidate for candidate in domains[other]
                        if all((candidate & prior).bit_count() != 1
                               for prior in chosen)
                    ]
                    if not domain:
                        forward_ok = False
                        break
                    remaining_domains.append(domain)
                if forward_ok and self._row_partial_possible(
                    chosen, remaining_domains
                ) and visit(following, chosen):
                    return True
            assignment[variable] = 0
            return False

        answer = tuple(assignment) if visit(tuple(range(count)), ()) else None
        self.cache[allowed] = answer
        return answer


@dataclass(frozen=True)
class ZeroForcingResult:
    number: int
    initial: int
    forces: tuple[tuple[int, int], ...]


class ZeroForcingSolver:
    """Exact ordinary zero forcing by initial-subset enumeration."""

    def __init__(self) -> None:
        self.cache: dict[tuple[int, ...], ZeroForcingResult] = {}
        self.cache_hits = 0
        self.initial_sets_checked = 0

    @staticmethod
    def closure(adj: Sequence[int], initial: int) -> tuple[int, tuple[tuple[int, int], ...]]:
        black = initial
        forces: list[tuple[int, int]] = []
        while True:
            made_force = False
            for vertex in range(len(adj)):
                if not (black & (1 << vertex)):
                    continue
                white = adj[vertex] & ~black
                if white and not (white & (white - 1)):
                    target = white.bit_length() - 1
                    black |= white
                    forces.append((vertex, target))
                    made_force = True
                    break
            if not made_force:
                return black, tuple(forces)

    def solve(self, adj: Sequence[int]) -> ZeroForcingResult:
        key = tuple(adj)
        if key in self.cache:
            self.cache_hits += 1
            return self.cache[key]
        n = len(adj)
        full = (1 << n) - 1
        for size in range(n + 1):
            for chosen in combinations(range(n), size):
                initial = sum(1 << v for v in chosen)
                self.initial_sets_checked += 1
                closure, forces = self.closure(adj, initial)
                if closure == full:
                    answer = ZeroForcingResult(size, initial, forces)
                    self.cache[key] = answer
                    return answer
        raise AssertionError("the full vertex set must be zero forcing")


def induced_graph(adj: Sequence[int], selected: Sequence[int]) -> tuple[int, ...]:
    """Relabel an induced graph consecutively."""

    position = {vertex: i for i, vertex in enumerate(selected)}
    return tuple(
        sum(1 << position[other] for other in selected
            if adj[vertex] & (1 << other))
        for vertex in selected
    )


def complement_graph(adj: Sequence[int]) -> tuple[int, ...]:
    full = (1 << len(adj)) - 1
    return tuple(full & ~row & ~(1 << vertex) for vertex, row in enumerate(adj))


def clique_number(adj: Sequence[int]) -> int:
    """Return the exact clique number of a small bit-mask graph."""

    for size in range(len(adj), 0, -1):
        if next(clique_masks(adj, size), None) is not None:
            return size
    return 0


def rational_nullspace(
    rows: Sequence[Sequence[int]], ncolumns: int
) -> tuple[tuple[Fraction, ...], ...]:
    """Return a deterministic exact rational basis for a matrix kernel."""

    matrix = [list(map(Fraction, row)) for row in rows]
    if any(len(row) != ncolumns for row in matrix):
        raise ValueError("inconsistent rational-nullspace row length")
    pivots: list[int] = []
    pivot_row = 0
    for column in range(ncolumns):
        selected = next(
            (
                row
                for row in range(pivot_row, len(matrix))
                if matrix[row][column]
            ),
            None,
        )
        if selected is None:
            continue
        matrix[pivot_row], matrix[selected] = (
            matrix[selected],
            matrix[pivot_row],
        )
        pivot = matrix[pivot_row][column]
        matrix[pivot_row] = [value / pivot for value in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row or not matrix[row][column]:
                continue
            multiplier = matrix[row][column]
            matrix[row] = [
                value - multiplier * basis_value
                for value, basis_value in zip(
                    matrix[row], matrix[pivot_row], strict=True
                )
            ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break

    free = [column for column in range(ncolumns) if column not in pivots]
    basis = []
    for free_column in free:
        vector = [Fraction(0)] * ncolumns
        vector[free_column] = Fraction(1)
        for row, column in enumerate(pivots):
            vector[column] = -matrix[row][free_column]
        basis.append(tuple(vector))
    return tuple(basis)


def basis_kernel_compatibility(
    adj: Sequence[int], clique_mask: int
) -> tuple[bool, dict | None, int]:
    """Test ``(D+A) ker(P)=0`` exactly for one saturating clique.

    The returned witness is JSON-friendly and records the first row that
    cannot admit a common rational diagonal scalar greater than one.
    """

    clique = list(bits(clique_mask))
    if any(
        not (adj[first] & (1 << second))
        for first, second in combinations(clique, 2)
    ):
        raise ValueError("basis-kernel test requires a clique")
    remainder = [
        vertex for vertex in range(len(adj))
        if not (clique_mask & (1 << vertex))
    ]
    p_matrix = [
        [int(bool(adj[vertex] & (1 << other))) for other in remainder]
        for vertex in clique
    ]
    kernel = rational_nullspace(p_matrix, len(remainder))
    dimension = len(kernel)
    if not kernel:
        return True, None, 0

    for local_row, vertex in enumerate(remainder):
        kernel_row = tuple(vector[local_row] for vector in kernel)
        action_row = tuple(
            sum(
                vector[local_column]
                for local_column, other in enumerate(remainder)
                if adj[vertex] & (1 << other)
            )
            for vector in kernel
        )
        nonzero = next(
            (column for column, value in enumerate(kernel_row) if value),
            None,
        )
        common = {
            "clique": clique,
            "row": vertex,
            "kernel_dimension": dimension,
        }
        if nonzero is None:
            if any(action_row):
                return False, {
                    **common,
                    "failure_kind": "zero_kernel_row_nonzero_action",
                }, dimension
            continue
        diagonal = -action_row[nonzero] / kernel_row[nonzero]
        if any(
            -action != diagonal * value
            for value, action in zip(kernel_row, action_row, strict=True)
        ):
            return False, {
                **common,
                "failure_kind": "no_common_diagonal_scalar",
            }, dimension
        if diagonal <= 1:
            return False, {
                **common,
                "failure_kind": "diagonal_not_greater_than_one",
                "forced_diagonal": [
                    diagonal.numerator,
                    diagonal.denominator,
                ],
            }, dimension
    return True, None, dimension


def saturating_clique_mask_compatibility(
    adj: Sequence[int], clique_mask: int
) -> tuple[bool, dict | None]:
    """Apply the exact Sherman--Morrison mask rules for one rank basis."""

    clique = list(bits(clique_mask))
    if any(
        not (adj[first] & (1 << second))
        for first, second in combinations(clique, 2)
    ):
        raise ValueError("saturating-mask test requires a clique")
    remainder = [
        vertex for vertex in range(len(adj))
        if not (clique_mask & (1 << vertex))
    ]
    masks = {vertex: adj[vertex] & clique_mask for vertex in remainder}
    for vertex in remainder:
        if not masks[vertex]:
            return False, {
                "clique": clique,
                "vertices": [vertex],
                "failure_kind": "empty_basis_neighbour_mask",
            }
        if masks[vertex] == clique_mask:
            return False, {
                "clique": clique,
                "vertices": [vertex],
                "failure_kind": "full_basis_neighbour_mask",
            }
    for first, second in combinations(remainder, 2):
        first_mask, second_mask = masks[first], masks[second]
        common = {"clique": clique, "vertices": [first, second]}
        if first_mask == second_mask:
            return False, {
                **common,
                "failure_kind": "duplicate_basis_neighbour_masks",
            }
        if not (first_mask & second_mask):
            return False, {
                **common,
                "failure_kind": "disjoint_basis_neighbour_masks",
            }
        required_edge = bool(adj[first] & (1 << second))
        if not required_edge and (
            not (first_mask & ~second_mask)
            or not (second_mask & ~first_mask)
        ):
            return False, {
                **common,
                "failure_kind": "nonedge_comparable_basis_masks",
            }
        if required_edge and first_mask | second_mask == clique_mask:
            return False, {
                **common,
                "failure_kind": "edge_basis_masks_cover_clique",
            }
    return True, None


@dataclass(frozen=True)
class CliqueStructureResult:
    clique_number: int
    f_maximum_degree: int
    saturating_cliques_checked: int
    basis_kernel_failures: int
    basis_kernel_maximum_nullity: int
    basis_kernel_failure_reasons: dict[str, int]
    basis_kernel_first_failure: dict | None
    saturating_mask_failures: int
    saturating_mask_failure_reasons: dict[str, int]
    saturating_mask_first_failure: dict | None


class CliqueStructureSolver:
    """Memoize exact clique, kernel, and Schur-mask tests on ``G[N]``."""

    def __init__(self) -> None:
        self.cache: dict[tuple[tuple[int, ...], int], CliqueStructureResult] = {}
        self.cache_hits = 0
        self.saturating_cliques_checked = 0
        self.rational_kernel_systems_checked = 0

    def solve(
        self, adj: Sequence[int], rank_upper: int
    ) -> CliqueStructureResult:
        key = (tuple(adj), rank_upper)
        if key in self.cache:
            self.cache_hits += 1
            return self.cache[key]
        fgraph = complement_graph(adj)
        maximum_degree = max((row.bit_count() for row in fgraph), default=0)
        saturating = (
            list(clique_masks(adj, rank_upper))
            if 0 < rank_upper <= len(adj)
            else []
        )
        basis_failures = 0
        mask_failures = 0
        basis_reasons: Counter[str] = Counter()
        mask_reasons: Counter[str] = Counter()
        maximum_nullity = 0
        first_basis = None
        first_mask = None
        for clique_mask in saturating:
            self.saturating_cliques_checked += 1
            compatible, witness, nullity = basis_kernel_compatibility(
                adj, clique_mask
            )
            self.rational_kernel_systems_checked += 1
            maximum_nullity = max(maximum_nullity, nullity)
            if not compatible:
                basis_failures += 1
                assert witness is not None
                basis_reasons[witness["failure_kind"]] += 1
                if first_basis is None:
                    first_basis = witness
            mask_compatible, mask_witness = (
                saturating_clique_mask_compatibility(adj, clique_mask)
            )
            if not mask_compatible:
                mask_failures += 1
                assert mask_witness is not None
                mask_reasons[mask_witness["failure_kind"]] += 1
                if first_mask is None:
                    first_mask = mask_witness
        result = CliqueStructureResult(
            clique_number=clique_number(adj),
            f_maximum_degree=maximum_degree,
            saturating_cliques_checked=len(saturating),
            basis_kernel_failures=basis_failures,
            basis_kernel_maximum_nullity=maximum_nullity,
            basis_kernel_failure_reasons=dict(basis_reasons),
            basis_kernel_first_failure=first_basis,
            saturating_mask_failures=mask_failures,
            saturating_mask_failure_reasons=dict(mask_reasons),
            saturating_mask_first_failure=first_mask,
        )
        self.cache[key] = result
        return result


def components(adj: Sequence[int]) -> list[tuple[int, ...]]:
    remaining = (1 << len(adj)) - 1
    output = []
    while remaining:
        queue = remaining & -remaining
        component = 0
        while queue:
            bit = queue & -queue
            queue ^= bit
            if component & bit:
                continue
            component |= bit
            vertex = bit.bit_length() - 1
            queue |= adj[vertex] & remaining & ~component
        remaining &= ~component
        output.append(tuple(bits(component)))
    return output


def component_inertia_nullity_caps(
    fgraph: Sequence[int], zero_forcing: ZeroForcingSolver
) -> tuple[dict[str, int], dict[str, int]]:
    """Return exact and cheap-tier component/inertia nullity bounds.

    A connected m-vertex graph always has the all-but-one zero-forcing set of
    size m-1.  Tiers 1--3 replace that fallback by the exact zero-forcing
    number only when it is at most the tier.  The second result is a
    conservative upper bound on the initial subsets a capped search examines.
    """

    nontrivial = [part for part in components(fgraph) if len(part) > 1]
    if not nontrivial:
        names = ("trivial", "tier1", "tier2", "tier3", "exact")
        return ({name: 0 for name in names}, {name: 0 for name in names})
    count = len(nontrivial)
    data = []
    for part in nontrivial:
        order = len(part)
        exact = zero_forcing.solve(induced_graph(fgraph, part)).number
        data.append((order, exact))
    caps: dict[str, int] = {}
    checks: dict[str, int] = {
        "trivial": 0,
        "exact": sum(
            sum(math.comb(order, size) for size in range(exact + 1))
            for order, exact in data
        ),
    }
    bounds = {
        "trivial": [order - 1 for order, _ in data],
        "exact": [exact for _, exact in data],
    }
    for tier in range(1, 4):
        name = f"tier{tier}"
        bounds[name] = [
            exact if exact <= tier else order - 1
            for order, exact in data
        ]
        checks[name] = sum(
            sum(math.comb(order, size) for size in range(min(tier, exact) + 1))
            for order, exact in data
        )
    for name, values in bounds.items():
        caps[name] = max(
            count,
            max(count - 1 + value for value in values),
        )
    return caps, checks


def component_inertia_nullity_cap(
    fgraph: Sequence[int], zero_forcing: ZeroForcingSolver
) -> int:
    """Compatibility helper returning the exact tier only."""

    return component_inertia_nullity_caps(fgraph, zero_forcing)[0]["exact"]


@dataclass
class CoverAnalysis:
    size: int
    direct_cap_failure: str | None
    support_failed: bool
    subspace_k_failed: bool
    pd_clique_failed: bool
    perpendicular_degree_failed: bool
    basis_kernel_failed: bool
    saturating_mask_failed: bool
    ordinary_b_failed: bool
    component_b_failed: bool
    joint_failed: bool
    enhanced_joint_failed: bool
    n_size: int
    term_rank: int
    total_term_rank: int
    k_rank_upper: int
    k_clique_number: int
    f_maximum_degree: int
    k_zero_forcing: int
    b_zero_forcing: int
    saturating_cliques_checked: int
    basis_kernel_clique_failures: int
    basis_kernel_maximum_nullity: int
    basis_kernel_failure_reasons: dict[str, int]
    basis_kernel_first_failure: dict | None
    saturating_mask_clique_failures: int
    saturating_mask_failure_reasons: dict[str, int]
    saturating_mask_first_failure: dict | None
    component_nullity_cap: int
    component_tier_failures: dict[str, bool]
    component_tier_check_upper_bounds: dict[str, int]


def rank_upper_bounds(
    n_size: int,
    term_rank: int,
    total_term_rank: int,
    z_size: int,
) -> tuple[int, int]:
    """Return the support/subspace upper bounds for rank(K) and rank(B)."""

    if not 0 <= z_size <= 7:
        raise ValueError("a K7 zero-factor cover has size between zero and seven")
    k_upper = min(
        n_size,
        term_rank,
        max(0, total_term_rank - z_size),
        7 - z_size,
    )
    return k_upper, min(n_size, k_upper + 1, 8 - z_size)


def direct_cover_cap_failure(z_size: int) -> str | None:
    """Return the separately proved direct cover-cap transition, if any.

    The complement of the orthonormal Z-family has dimension ``7-z_size``.
    Size seven leaves no dimension for the nonempty N family.  Size six puts
    all nonzero normalized columns on one line, incompatible with K off-
    diagonals 0/1 and diagonals greater than one.  At size five, the
    orthogonality graph F is a matching in a plane, forcing rank(B)>=4>3.
    At size four, F has maximum degree two; the component/inertia equality
    cases reduce to 4K2 and C4+2K2, both excluded respectively by PSD rank
    monotonicity for K=B+J and by the opposite-pair geometry in R3.
    """

    return {
        7: "cap7_to6",
        6: "cap6_to5",
        5: "cap5_to4",
        4: "cap4_to3",
    }.get(z_size)


def analyze_cover(
    adj: Sequence[int],
    outside: Sequence[int],
    defects: Sequence[int],
    zmask: int,
    support_solver: SupportSolver,
    zero_forcing: ZeroForcingSolver,
    total_term_rank: int | None = None,
    clique_solver: CliqueStructureSolver | None = None,
) -> CoverAnalysis:
    """Apply all exact support/rank screens to one eligible cover."""

    zvertices = list(bits(zmask))
    support_failed = support_solver.solve(defects[i] for i in zvertices) is None
    nvertices = [i for i in range(len(outside)) if not (zmask & (1 << i))]
    ndefects = [defects[i] for i in nvertices]
    term_rank = matching_size(ndefects)
    if total_term_rank is None:
        total_term_rank = matching_size(defects)
    global_n = [outside[i] for i in nvertices]
    graph_n = induced_graph(adj, global_n)

    # Cover + alpha(G)<=2 makes adjacency exactly mask intersection on N.
    for i, j in combinations(range(len(nvertices)), 2):
        intersection = bool(ndefects[i] & ndefects[j])
        required_edge = bool(graph_n[i] & (1 << j))
        assert intersection == required_edge

    k_zf = zero_forcing.solve(graph_n)
    k_upper, b_upper = rank_upper_bounds(
        len(nvertices), term_rank, total_term_rank, len(zvertices)
    )
    if clique_solver is None:
        clique_solver = CliqueStructureSolver()
    structure = clique_solver.solve(graph_n, k_upper)
    # The term-rank part alone is a proved-redundant control: the 0/1 mask-
    # incidence Gram matrix has this graph and rank at most term_rank.  The
    # new 7-|Z| subspace bound is genuinely stronger.
    assert len(nvertices) - k_zf.number <= min(
        len(nvertices), term_rank
    )
    subspace_k_failed = (
        bool(nvertices) and len(zvertices) == 7
    ) or len(nvertices) - k_zf.number > k_upper
    pd_clique_failed = structure.clique_number > k_upper
    perpendicular_degree_failed = bool(nvertices) and (
        structure.f_maximum_degree > k_upper - 1
    )
    basis_kernel_failed = bool(structure.basis_kernel_failures)
    saturating_mask_failed = bool(structure.saturating_mask_failures)

    fgraph = complement_graph(graph_n)
    b_zf = zero_forcing.solve(fgraph)
    ordinary_failed = len(nvertices) - b_zf.number > b_upper

    nullity_caps, tier_checks = component_inertia_nullity_caps(
        fgraph, zero_forcing
    )
    nullity_cap = nullity_caps["exact"]
    component_failed = len(nvertices) - nullity_cap > b_upper
    tier_failures = {
        name: len(nvertices) - cap > b_upper
        for name, cap in nullity_caps.items()
    }
    assert tier_failures["exact"] == component_failed
    assert not ordinary_failed or component_failed
    cap_failure = direct_cover_cap_failure(len(zvertices))
    return CoverAnalysis(
        size=len(zvertices),
        direct_cap_failure=cap_failure,
        support_failed=support_failed,
        subspace_k_failed=subspace_k_failed,
        pd_clique_failed=pd_clique_failed,
        perpendicular_degree_failed=perpendicular_degree_failed,
        basis_kernel_failed=basis_kernel_failed,
        saturating_mask_failed=saturating_mask_failed,
        ordinary_b_failed=ordinary_failed,
        component_b_failed=component_failed,
        joint_failed=(
            cap_failure is not None
            or support_failed
            or subspace_k_failed
            or component_failed
        ),
        enhanced_joint_failed=(
            cap_failure is not None
            or support_failed
            or subspace_k_failed
            or component_failed
            or pd_clique_failed
            or perpendicular_degree_failed
            or basis_kernel_failed
            or saturating_mask_failed
        ),
        n_size=len(nvertices),
        term_rank=term_rank,
        total_term_rank=total_term_rank,
        k_rank_upper=k_upper,
        k_clique_number=structure.clique_number,
        f_maximum_degree=structure.f_maximum_degree,
        k_zero_forcing=k_zf.number,
        b_zero_forcing=b_zf.number,
        saturating_cliques_checked=structure.saturating_cliques_checked,
        basis_kernel_clique_failures=structure.basis_kernel_failures,
        basis_kernel_maximum_nullity=(
            structure.basis_kernel_maximum_nullity
        ),
        basis_kernel_failure_reasons=(
            structure.basis_kernel_failure_reasons
        ),
        basis_kernel_first_failure=structure.basis_kernel_first_failure,
        saturating_mask_clique_failures=structure.saturating_mask_failures,
        saturating_mask_failure_reasons=(
            structure.saturating_mask_failure_reasons
        ),
        saturating_mask_first_failure=(
            structure.saturating_mask_first_failure
        ),
        component_nullity_cap=nullity_cap,
        component_tier_failures=tier_failures,
        component_tier_check_upper_bounds=tier_checks,
    )


@dataclass
class SeedAnalysis:
    seed: list[int]
    covers: int
    cover_sizes: dict[str, int]
    direct_cap_failures: dict[str, int]
    support_failures: int
    subspace_k_failures: int
    pd_clique_failures: int
    perpendicular_degree_failures: int
    basis_kernel_failures: int
    saturating_mask_failures: int
    saturating_cliques_checked: int
    basis_kernel_clique_failures: int
    saturating_mask_clique_failures: int
    basis_kernel_failure_reasons: dict[str, int]
    saturating_mask_failure_reasons: dict[str, int]
    ordinary_b_failures: int
    component_b_failures: int
    component_tier_failures: dict[str, int]
    marginal_support_failures: int
    marginal_subspace_k_failures: int
    marginal_ordinary_b_failures: int
    marginal_component_b_failures: int
    marginal_pd_clique_failures: int
    marginal_perpendicular_degree_failures: int
    marginal_basis_kernel_failures: int
    marginal_saturating_mask_failures: int
    marginal_component_tier_failures: dict[str, int]
    component_tier_check_upper_bounds: dict[str, int]
    joint_passes: int
    enhanced_joint_passes: int
    old_cover_rejected: bool
    cap7_to6_rejected: bool
    cap6_to5_rejected: bool
    cap5_to4_rejected: bool
    cap4_to3_rejected: bool
    support_only_rejected: bool
    subspace_k_only_rejected: bool
    pd_clique_only_rejected: bool
    perpendicular_degree_only_rejected: bool
    basis_kernel_only_rejected: bool
    saturating_mask_only_rejected: bool
    ordinary_b_only_rejected: bool
    component_b_only_rejected: bool
    component_tier_only_rejected: dict[str, bool]
    joint_rejected: bool
    joint_pd_clique_rejected: bool
    joint_pd_clique_degree_rejected: bool
    joint_pd_clique_degree_basis_rejected: bool
    enhanced_joint_rejected: bool


def analyze_seed(
    adj: Sequence[int],
    seed_mask: int,
    support_solver: SupportSolver,
    zero_forcing: ZeroForcingSolver,
    clique_solver: CliqueStructureSolver,
) -> SeedAnalysis:
    seed, outside, defects, ladj, eligible = seed_instance(adj, seed_mask)
    covers = eligible_covers(ladj, eligible)
    total_term_rank = matching_size(defects)
    analyses = [
        analyze_cover(
            adj,
            outside,
            defects,
            cover,
            support_solver,
            zero_forcing,
            total_term_rank,
            clique_solver,
        )
        for cover in covers
    ]
    support_pass = [entry for entry in analyses if not entry.support_failed]
    subspace_pass = [entry for entry in analyses if not entry.subspace_k_failed]
    pd_clique_pass = [entry for entry in analyses if not entry.pd_clique_failed]
    perpendicular_degree_pass = [
        entry for entry in analyses
        if not entry.perpendicular_degree_failed
    ]
    basis_kernel_pass = [
        entry for entry in analyses if not entry.basis_kernel_failed
    ]
    saturating_mask_pass = [
        entry for entry in analyses if not entry.saturating_mask_failed
    ]
    ordinary_pass = [entry for entry in analyses if not entry.ordinary_b_failed]
    component_pass = [entry for entry in analyses if not entry.component_b_failed]
    direct_pass = [
        entry for entry in analyses if entry.direct_cap_failure is None
    ]
    after_direct_support = [
        entry for entry in direct_pass if not entry.support_failed
    ]
    after_support_subspace = [
        entry for entry in after_direct_support
        if not entry.subspace_k_failed
    ]
    after_support_subspace_ordinary = [
        entry for entry in after_support_subspace
        if not entry.ordinary_b_failed
    ]
    joint_pass = [
        entry for entry in after_support_subspace
        if not entry.component_b_failed
    ]
    after_current_pd_clique = [
        entry for entry in joint_pass if not entry.pd_clique_failed
    ]
    after_perpendicular_degree = [
        entry for entry in after_current_pd_clique
        if not entry.perpendicular_degree_failed
    ]
    after_basis_kernel = [
        entry for entry in after_perpendicular_degree
        if not entry.basis_kernel_failed
    ]
    enhanced_joint_pass = [
        entry for entry in after_basis_kernel
        if not entry.saturating_mask_failed
    ]
    basis_failure_reasons: Counter[str] = Counter()
    mask_failure_reasons: Counter[str] = Counter()
    for entry in analyses:
        basis_failure_reasons.update(entry.basis_kernel_failure_reasons)
        mask_failure_reasons.update(entry.saturating_mask_failure_reasons)
    return SeedAnalysis(
        seed=seed,
        covers=len(covers),
        cover_sizes={
            str(size): sum(entry.size == size for entry in analyses)
            for size in range(8)
            if any(entry.size == size for entry in analyses)
        },
        direct_cap_failures={
            name: sum(entry.direct_cap_failure == name for entry in analyses)
            for name in ("cap7_to6", "cap6_to5", "cap5_to4", "cap4_to3")
        },
        support_failures=sum(entry.support_failed for entry in analyses),
        subspace_k_failures=sum(
            entry.subspace_k_failed for entry in analyses
        ),
        pd_clique_failures=sum(entry.pd_clique_failed for entry in analyses),
        perpendicular_degree_failures=sum(
            entry.perpendicular_degree_failed for entry in analyses
        ),
        basis_kernel_failures=sum(
            entry.basis_kernel_failed for entry in analyses
        ),
        saturating_mask_failures=sum(
            entry.saturating_mask_failed for entry in analyses
        ),
        saturating_cliques_checked=sum(
            entry.saturating_cliques_checked for entry in analyses
        ),
        basis_kernel_clique_failures=sum(
            entry.basis_kernel_clique_failures for entry in analyses
        ),
        saturating_mask_clique_failures=sum(
            entry.saturating_mask_clique_failures for entry in analyses
        ),
        basis_kernel_failure_reasons=dict(basis_failure_reasons),
        saturating_mask_failure_reasons=dict(mask_failure_reasons),
        ordinary_b_failures=sum(entry.ordinary_b_failed for entry in analyses),
        component_b_failures=sum(entry.component_b_failed for entry in analyses),
        component_tier_failures={
            name: sum(entry.component_tier_failures[name] for entry in analyses)
            for name in ("trivial", "tier1", "tier2", "tier3", "exact")
        },
        marginal_support_failures=sum(
            entry.support_failed for entry in direct_pass
        ),
        marginal_subspace_k_failures=sum(
            entry.subspace_k_failed for entry in after_direct_support
        ),
        marginal_ordinary_b_failures=sum(
            entry.ordinary_b_failed for entry in after_support_subspace
        ),
        marginal_component_b_failures=sum(
            entry.component_b_failed
            for entry in after_support_subspace_ordinary
        ),
        marginal_pd_clique_failures=sum(
            entry.pd_clique_failed for entry in joint_pass
        ),
        marginal_perpendicular_degree_failures=sum(
            entry.perpendicular_degree_failed
            for entry in after_current_pd_clique
        ),
        marginal_basis_kernel_failures=sum(
            entry.basis_kernel_failed for entry in after_perpendicular_degree
        ),
        marginal_saturating_mask_failures=sum(
            entry.saturating_mask_failed for entry in after_basis_kernel
        ),
        marginal_component_tier_failures={
            name: sum(
                entry.component_tier_failures[name]
                for entry in after_support_subspace_ordinary
            )
            for name in ("trivial", "tier1", "tier2", "tier3", "exact")
        },
        component_tier_check_upper_bounds={
            name: sum(
                entry.component_tier_check_upper_bounds[name]
                for entry in analyses
            )
            for name in ("trivial", "tier1", "tier2", "tier3", "exact")
        },
        joint_passes=len(joint_pass),
        enhanced_joint_passes=len(enhanced_joint_pass),
        old_cover_rejected=not analyses,
        cap7_to6_rejected=bool(analyses) and not any(
            entry.size <= 6 for entry in analyses
        ),
        cap6_to5_rejected=bool(analyses) and not any(
            entry.size <= 5 for entry in analyses
        ),
        cap5_to4_rejected=bool(analyses) and not any(
            entry.size <= 4 for entry in analyses
        ),
        cap4_to3_rejected=bool(analyses) and not direct_pass,
        support_only_rejected=bool(analyses) and not support_pass,
        subspace_k_only_rejected=bool(analyses) and not subspace_pass,
        pd_clique_only_rejected=bool(analyses) and not pd_clique_pass,
        perpendicular_degree_only_rejected=(
            bool(analyses) and not perpendicular_degree_pass
        ),
        basis_kernel_only_rejected=bool(analyses) and not basis_kernel_pass,
        saturating_mask_only_rejected=(
            bool(analyses) and not saturating_mask_pass
        ),
        ordinary_b_only_rejected=bool(analyses) and not ordinary_pass,
        component_b_only_rejected=bool(analyses) and not component_pass,
        component_tier_only_rejected={
            name: bool(analyses) and not any(
                not entry.component_tier_failures[name] for entry in analyses
            )
            for name in ("trivial", "tier1", "tier2", "tier3", "exact")
        },
        joint_rejected=bool(analyses) and not joint_pass,
        joint_pd_clique_rejected=(
            bool(analyses) and not after_current_pd_clique
        ),
        joint_pd_clique_degree_rejected=(
            bool(analyses) and not after_perpendicular_degree
        ),
        joint_pd_clique_degree_basis_rejected=(
            bool(analyses) and not after_basis_kernel
        ),
        enhanced_joint_rejected=bool(analyses) and not enhanced_joint_pass,
    )


@dataclass
class GraphAnalysis:
    index: int | None
    stratum: str | None
    applicable: bool
    seeds: int
    covers: int
    elapsed_seconds: float
    old_cover_rejected: bool
    cap7_to6_rejected: bool
    cap6_to5_rejected: bool
    cap5_to4_rejected: bool
    cap4_to3_rejected: bool
    support_only_rejected: bool
    subspace_k_only_rejected: bool
    pd_clique_only_rejected: bool
    perpendicular_degree_only_rejected: bool
    basis_kernel_only_rejected: bool
    saturating_mask_only_rejected: bool
    ordinary_b_only_rejected: bool
    component_b_only_rejected: bool
    component_tier_only_rejected: dict[str, bool]
    component_tier_check_upper_bounds: dict[str, int]
    joint_rejected: bool
    incremental_joint_rejected: bool
    joint_pd_clique_rejected: bool
    joint_pd_clique_degree_rejected: bool
    joint_pd_clique_degree_basis_rejected: bool
    enhanced_joint_rejected: bool
    incremental_enhanced_over_current_joint: bool
    saturating_cliques_checked: int
    basis_kernel_clique_failures: int
    saturating_mask_clique_failures: int
    basis_kernel_failure_reasons: dict[str, int]
    saturating_mask_failure_reasons: dict[str, int]
    cover_sizes: dict[str, int] = field(default_factory=dict)
    cover_failures: dict[str, int] = field(default_factory=dict)
    marginal_cover_failures: dict[str, int] = field(default_factory=dict)
    first_rejecting_seeds: dict[str, list[int]] = field(default_factory=dict)


def analyze_graph(
    adj: Sequence[int],
    support_solver: SupportSolver,
    zero_forcing: ZeroForcingSolver,
    clique_solver: CliqueStructureSolver,
    *,
    index: int | None = None,
    stratum: str | None = None,
) -> GraphAnalysis:
    validate_graph(adj)
    started = time.perf_counter()
    seeds = [
        analyze_seed(
            adj,
            seed_mask,
            support_solver,
            zero_forcing,
            clique_solver,
        )
        for seed_mask in clique_masks(adj, 7)
    ]
    cover_sizes: Counter[str] = Counter()
    failures: Counter[str] = Counter()
    marginal: Counter[str] = Counter()
    for seed in seeds:
        cover_sizes.update(seed.cover_sizes)
        failures.update(seed.direct_cap_failures)
        failures.update(
            {f"component_B_{name}": value
             for name, value in seed.component_tier_failures.items()}
        )
        failures.update(
            support=seed.support_failures,
            subspace_K=seed.subspace_k_failures,
            PD_clique=seed.pd_clique_failures,
            perpendicular_degree=seed.perpendicular_degree_failures,
            basis_kernel=seed.basis_kernel_failures,
            saturating_masks=seed.saturating_mask_failures,
            basis_kernel_cliques=seed.basis_kernel_clique_failures,
            saturating_mask_cliques=seed.saturating_mask_clique_failures,
            ordinary_B=seed.ordinary_b_failures,
            component_B=seed.component_b_failures,
        )
        failures.update(
            {
                f"basis_kernel_reason_{name}": value
                for name, value in seed.basis_kernel_failure_reasons.items()
            }
        )
        failures.update(
            {
                f"saturating_mask_reason_{name}": value
                for name, value in seed.saturating_mask_failure_reasons.items()
            }
        )
        marginal.update(
            cap7_to6=seed.direct_cap_failures["cap7_to6"],
            cap6_to5=seed.direct_cap_failures["cap6_to5"],
            cap5_to4=seed.direct_cap_failures["cap5_to4"],
            cap4_to3=seed.direct_cap_failures["cap4_to3"],
            support_after_cap3=seed.marginal_support_failures,
            subspace_K_after_cap3_and_support=(
                seed.marginal_subspace_k_failures
            ),
            ordinary_B_after_cap3_support_and_subspace=(
                seed.marginal_ordinary_b_failures
            ),
            component_B_after_cap3_support_subspace_and_ordinary=(
                seed.marginal_component_b_failures
            ),
            PD_clique_after_current_joint=seed.marginal_pd_clique_failures,
            perpendicular_degree_after_current_joint_and_PD_clique=(
                seed.marginal_perpendicular_degree_failures
            ),
            basis_kernel_after_current_joint_PD_clique_and_degree=(
                seed.marginal_basis_kernel_failures
            ),
            saturating_masks_after_current_joint_and_prior_new_rules=(
                seed.marginal_saturating_mask_failures
            ),
        )
        marginal.update(
            {
                f"component_B_{name}_after_prior": value
                for name, value in seed.marginal_component_tier_failures.items()
            }
        )
    predicates = {
        "old_cover": lambda seed: seed.old_cover_rejected,
        "cap7_to6": lambda seed: seed.cap7_to6_rejected,
        "cap6_to5": lambda seed: seed.cap6_to5_rejected,
        "cap5_to4": lambda seed: seed.cap5_to4_rejected,
        "cap4_to3": lambda seed: seed.cap4_to3_rejected,
        "support": lambda seed: seed.support_only_rejected,
        "subspace_K": lambda seed: seed.subspace_k_only_rejected,
        "PD_clique": lambda seed: seed.pd_clique_only_rejected,
        "perpendicular_degree": (
            lambda seed: seed.perpendicular_degree_only_rejected
        ),
        "basis_kernel": lambda seed: seed.basis_kernel_only_rejected,
        "saturating_masks": lambda seed: seed.saturating_mask_only_rejected,
        "ordinary_B": lambda seed: seed.ordinary_b_only_rejected,
        "component_B": lambda seed: seed.component_b_only_rejected,
        "joint": lambda seed: seed.joint_rejected,
        "joint_PD_clique": lambda seed: seed.joint_pd_clique_rejected,
        "joint_PD_clique_degree": (
            lambda seed: seed.joint_pd_clique_degree_rejected
        ),
        "joint_PD_clique_degree_basis": (
            lambda seed: seed.joint_pd_clique_degree_basis_rejected
        ),
        "enhanced_joint": lambda seed: seed.enhanced_joint_rejected,
    }
    first = {
        name: next((seed.seed for seed in seeds if predicate(seed)), [])
        for name, predicate in predicates.items()
    }
    old_rejected = bool(first["old_cover"])
    joint_rejected = bool(first["joint"])
    joint_pd_clique_rejected = bool(first["joint_PD_clique"])
    joint_pd_clique_degree_rejected = bool(first["joint_PD_clique_degree"])
    joint_pd_clique_degree_basis_rejected = bool(
        first["joint_PD_clique_degree_basis"]
    )
    enhanced_joint_rejected = bool(first["enhanced_joint"])
    return GraphAnalysis(
        index=index,
        stratum=stratum,
        applicable=bool(seeds),
        seeds=len(seeds),
        covers=sum(seed.covers for seed in seeds),
        elapsed_seconds=time.perf_counter() - started,
        old_cover_rejected=old_rejected,
        cap7_to6_rejected=bool(first["cap7_to6"]),
        cap6_to5_rejected=bool(first["cap6_to5"]),
        cap5_to4_rejected=bool(first["cap5_to4"]),
        cap4_to3_rejected=bool(first["cap4_to3"]),
        support_only_rejected=bool(first["support"]),
        subspace_k_only_rejected=bool(first["subspace_K"]),
        pd_clique_only_rejected=bool(first["PD_clique"]),
        perpendicular_degree_only_rejected=bool(
            first["perpendicular_degree"]
        ),
        basis_kernel_only_rejected=bool(first["basis_kernel"]),
        saturating_mask_only_rejected=bool(first["saturating_masks"]),
        ordinary_b_only_rejected=bool(first["ordinary_B"]),
        component_b_only_rejected=bool(first["component_B"]),
        component_tier_only_rejected={
            name: any(
                seed.component_tier_only_rejected[name] for seed in seeds
            )
            for name in ("trivial", "tier1", "tier2", "tier3", "exact")
        },
        component_tier_check_upper_bounds={
            name: sum(
                seed.component_tier_check_upper_bounds[name] for seed in seeds
            )
            for name in ("trivial", "tier1", "tier2", "tier3", "exact")
        },
        joint_rejected=joint_rejected,
        incremental_joint_rejected=joint_rejected and not old_rejected,
        joint_pd_clique_rejected=joint_pd_clique_rejected,
        joint_pd_clique_degree_rejected=(
            joint_pd_clique_degree_rejected
        ),
        joint_pd_clique_degree_basis_rejected=(
            joint_pd_clique_degree_basis_rejected
        ),
        enhanced_joint_rejected=enhanced_joint_rejected,
        incremental_enhanced_over_current_joint=(
            enhanced_joint_rejected and not joint_rejected
        ),
        saturating_cliques_checked=sum(
            seed.saturating_cliques_checked for seed in seeds
        ),
        basis_kernel_clique_failures=sum(
            seed.basis_kernel_clique_failures for seed in seeds
        ),
        saturating_mask_clique_failures=sum(
            seed.saturating_mask_clique_failures for seed in seeds
        ),
        basis_kernel_failure_reasons=dict(
            sum(
                (Counter(seed.basis_kernel_failure_reasons) for seed in seeds),
                Counter(),
            )
        ),
        saturating_mask_failure_reasons=dict(
            sum(
                (
                    Counter(seed.saturating_mask_failure_reasons)
                    for seed in seeds
                ),
                Counter(),
            )
        ),
        cover_sizes=dict(sorted(cover_sizes.items(), key=lambda item: int(item[0]))),
        cover_failures=dict(failures),
        marginal_cover_failures=dict(marginal),
        first_rejecting_seeds={name: value for name, value in first.items() if value},
    )


def solver_counters(
    support_solver: SupportSolver,
    zero_forcing: ZeroForcingSolver,
    clique_solver: CliqueStructureSolver,
) -> dict[str, int]:
    return {
        "support_CSP_nodes": support_solver.nodes,
        "support_cache_entries": len(support_solver.cache),
        "support_cache_hits": support_solver.cache_hits,
        "zero_forcing_initial_sets_checked": zero_forcing.initial_sets_checked,
        "zero_forcing_cache_entries": len(zero_forcing.cache),
        "zero_forcing_cache_hits": zero_forcing.cache_hits,
        "clique_structure_cache_entries": len(clique_solver.cache),
        "clique_structure_cache_hits": clique_solver.cache_hits,
        "saturating_cliques_actually_checked": (
            clique_solver.saturating_cliques_checked
        ),
        "rational_kernel_systems_actually_checked": (
            clique_solver.rational_kernel_systems_checked
        ),
    }


def isolated_graph_analysis(graph: dict) -> tuple[GraphAnalysis, dict[str, int]]:
    """Process-worker entry point with graph-local transparent caches."""

    support_solver = SupportSolver()
    zero_forcing = ZeroForcingSolver()
    clique_solver = CliqueStructureSolver()
    result = analyze_graph(
        graph["adjacency"],
        support_solver,
        zero_forcing,
        clique_solver,
        index=graph.get("index"),
        stratum=graph.get("stratum"),
    )
    return result, solver_counters(
        support_solver, zero_forcing, clique_solver
    )


def sample_report(
    sample: dict, *, workers: int = 1, include_per_graph: bool = True
) -> dict:
    if workers < 1:
        raise ValueError("workers must be positive")
    started = time.perf_counter()
    counter_totals: Counter[str] = Counter()
    if workers == 1:
        support_solver = SupportSolver()
        zero_forcing = ZeroForcingSolver()
        clique_solver = CliqueStructureSolver()
        graphs = [
            analyze_graph(
                graph["adjacency"],
                support_solver,
                zero_forcing,
                clique_solver,
                index=graph.get("index"),
                stratum=graph.get("stratum"),
            )
            for graph in sample["graphs"]
        ]
        counter_totals.update(
            solver_counters(support_solver, zero_forcing, clique_solver)
        )
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = list(
                executor.map(
                    isolated_graph_analysis,
                    sample["graphs"],
                    chunksize=1,
                )
            )
        graphs = [result for result, _ in results]
        for _, counters in results:
            counter_totals.update(counters)
    elapsed = time.perf_counter() - started
    aggregate_cover_sizes: Counter[str] = Counter()
    aggregate_failures: Counter[str] = Counter()
    aggregate_marginal: Counter[str] = Counter()
    tier_check_totals: Counter[str] = Counter()
    by_stratum: dict[str, Counter[str]] = {}
    for graph in graphs:
        aggregate_cover_sizes.update(graph.cover_sizes)
        aggregate_failures.update(graph.cover_failures)
        aggregate_marginal.update(graph.marginal_cover_failures)
        tier_check_totals.update(graph.component_tier_check_upper_bounds)
        counts = by_stratum.setdefault(graph.stratum or "unspecified", Counter())
        counts.update(
            graphs=1,
            applicable=int(graph.applicable),
            old_cover_rejected=int(graph.old_cover_rejected),
            cap7_to6_rejected=int(graph.cap7_to6_rejected),
            cap6_to5_rejected=int(graph.cap6_to5_rejected),
            cap5_to4_rejected=int(graph.cap5_to4_rejected),
            cap4_to3_rejected=int(graph.cap4_to3_rejected),
            support_only_rejected=int(graph.support_only_rejected),
            subspace_K_only_rejected=int(graph.subspace_k_only_rejected),
            PD_clique_only_rejected=int(graph.pd_clique_only_rejected),
            perpendicular_degree_only_rejected=int(
                graph.perpendicular_degree_only_rejected
            ),
            basis_kernel_only_rejected=int(
                graph.basis_kernel_only_rejected
            ),
            saturating_mask_only_rejected=int(
                graph.saturating_mask_only_rejected
            ),
            ordinary_B_only_rejected=int(graph.ordinary_b_only_rejected),
            component_B_only_rejected=int(graph.component_b_only_rejected),
            joint_rejected=int(graph.joint_rejected),
            incremental_joint_rejected=int(graph.incremental_joint_rejected),
            joint_PD_clique_rejected=int(
                graph.joint_pd_clique_rejected
            ),
            joint_PD_clique_degree_rejected=int(
                graph.joint_pd_clique_degree_rejected
            ),
            joint_PD_clique_degree_basis_rejected=int(
                graph.joint_pd_clique_degree_basis_rejected
            ),
            enhanced_joint_rejected=int(graph.enhanced_joint_rejected),
            incremental_enhanced_over_current_joint=int(
                graph.incremental_enhanced_over_current_joint
            ),
        )
        for name, rejected in graph.component_tier_only_rejected.items():
            counts[f"component_B_{name}_only_rejected"] += int(rejected)

    def selected_indices(attribute: str) -> list[int]:
        return [
            graph.index
            for graph in graphs
            if graph.index is not None and bool(getattr(graph, attribute))
        ]

    old_indices = set(selected_indices("old_cover_rejected"))
    cap7_to6_indices = set(selected_indices("cap7_to6_rejected")) - old_indices
    cap6_to5_indices = set(selected_indices("cap6_to5_rejected")) - old_indices
    cap5_to4_indices = set(selected_indices("cap5_to4_rejected")) - old_indices
    cap4_to3_indices = set(selected_indices("cap4_to3_rejected")) - old_indices
    support_indices = set(selected_indices("support_only_rejected")) - old_indices
    subspace_indices = set(
        selected_indices("subspace_k_only_rejected")
    ) - old_indices
    ordinary_indices = set(selected_indices("ordinary_b_only_rejected")) - old_indices
    component_indices = set(selected_indices("component_b_only_rejected")) - old_indices
    joint_indices = set(selected_indices("joint_rejected")) - old_indices
    joint_pd_clique_indices = set(
        selected_indices("joint_pd_clique_rejected")
    ) - old_indices
    joint_pd_clique_degree_indices = set(
        selected_indices("joint_pd_clique_degree_rejected")
    ) - old_indices
    joint_pd_clique_degree_basis_indices = set(
        selected_indices("joint_pd_clique_degree_basis_rejected")
    ) - old_indices
    pd_clique_indices = set(
        selected_indices("pd_clique_only_rejected")
    ) - old_indices
    perpendicular_degree_indices = set(
        selected_indices("perpendicular_degree_only_rejected")
    ) - old_indices
    basis_kernel_indices = set(
        selected_indices("basis_kernel_only_rejected")
    ) - old_indices
    saturating_mask_indices = set(
        selected_indices("saturating_mask_only_rejected")
    ) - old_indices
    enhanced_joint_indices = set(
        selected_indices("enhanced_joint_rejected")
    ) - old_indices
    marginal_cap7_to6 = cap7_to6_indices
    marginal_cap6_to5 = cap6_to5_indices - marginal_cap7_to6
    marginal_cap5_to4 = (
        cap5_to4_indices - marginal_cap7_to6 - marginal_cap6_to5
    )
    marginal_cap4_to3 = (
        cap4_to3_indices
        - marginal_cap7_to6
        - marginal_cap6_to5
        - marginal_cap5_to4
    )
    cap_marginal_union = (
        marginal_cap7_to6
        | marginal_cap6_to5
        | marginal_cap5_to4
        | marginal_cap4_to3
    )
    marginal_support = support_indices - cap_marginal_union
    marginal_subspace = subspace_indices - cap_marginal_union - marginal_support
    marginal_ordinary = (
        ordinary_indices
        - cap_marginal_union
        - marginal_support
        - marginal_subspace
    )
    marginal_component = (
        component_indices
        - cap_marginal_union
        - marginal_support
        - marginal_subspace
        - marginal_ordinary
    )
    marginal_synergy = (
        joint_indices
        - cap_marginal_union
        - marginal_support
        - marginal_subspace
        - marginal_ordinary
        - marginal_component
    )
    marginal_pd_clique = pd_clique_indices - joint_indices
    marginal_perpendicular_degree = (
        perpendicular_degree_indices
        - joint_indices
        - marginal_pd_clique
    )
    marginal_basis_kernel = (
        basis_kernel_indices
        - joint_indices
        - marginal_pd_clique
        - marginal_perpendicular_degree
    )
    marginal_saturating_masks = (
        saturating_mask_indices
        - joint_indices
        - marginal_pd_clique
        - marginal_perpendicular_degree
        - marginal_basis_kernel
    )
    enhanced_synergy = (
        enhanced_joint_indices
        - joint_indices
        - marginal_pd_clique
        - marginal_perpendicular_degree
        - marginal_basis_kernel
        - marginal_saturating_masks
    )
    if not (
        joint_indices
        <= joint_pd_clique_indices
        <= joint_pd_clique_degree_indices
        <= joint_pd_clique_degree_basis_indices
        <= enhanced_joint_indices
    ):
        raise AssertionError("enhanced existential decisions must be nested")
    existential_pd_clique = joint_pd_clique_indices - joint_indices
    existential_perpendicular_degree = (
        joint_pd_clique_degree_indices - joint_pd_clique_indices
    )
    existential_basis_kernel = (
        joint_pd_clique_degree_basis_indices
        - joint_pd_clique_degree_indices
    )
    existential_saturating_masks = (
        enhanced_joint_indices - joint_pd_clique_degree_basis_indices
    )
    component_tier_indices = {
        name: {
            graph.index
            for graph in graphs
            if graph.index is not None
            and graph.component_tier_only_rejected[name]
            and graph.index not in old_indices
        }
        for name in ("trivial", "tier1", "tier2", "tier3", "exact")
    }
    prior_individual_union = (
        cap_marginal_union
        | support_indices
        | subspace_indices
        | ordinary_indices
    )
    component_tier_marginal = {
        name: indices - prior_individual_union
        for name, indices in component_tier_indices.items()
    }
    survivor_indices = sorted(
        graph.index
        for graph in graphs
        if graph.index is not None
        and graph.applicable
        and graph.index not in old_indices
        and graph.index not in joint_indices
    )
    enhanced_survivor_indices = sorted(
        graph.index
        for graph in graphs
        if graph.index is not None
        and graph.applicable
        and graph.index not in old_indices
        and graph.index not in enhanced_joint_indices
    )
    report = {
        "schema": 1,
        "description": (
            "Independent exact K7 existential-cover support, K/B rank, "
            "positive-definite clique, rational basis-kernel, and "
            "saturating-clique Schur-mask reference. Existing pre-extension "
            "decisions are retained separately."
        ),
        "exact_K_structure_rules": {
            "rank_upper": "U_K=min(7-|Z|,nu_N,nu_T-|Z|)",
            "positive_definite_clique": "omega(G[N])<=U_K",
            "perpendicular_neighbourhood": "Delta(complement(G[N]))<=U_K-1",
            "basis_kernel": (
                "for every size-U_K clique C, there is a diagonal D>I "
                "with (D+A_R)ker(K[C,R])=0; tested over exact rationals"
            ),
            "saturating_clique_masks": (
                "every outside basis-neighbour mask is nonempty/proper; "
                "masks are distinct and intersecting; nonedge masks are "
                "incomparable; edge-mask unions are proper"
            ),
            "candidate_nonedge_semantics": (
                "candidate nonedges remain optional unit pairs; the exact "
                "G[N] zero pattern follows only after fixing a cover Z"
            ),
        },
        "graphs": len(graphs),
        "worker_count": workers,
        "applicable_K7_graphs": sum(graph.applicable for graph in graphs),
        "seeds_checked": sum(graph.seeds for graph in graphs),
        "eligible_covers_checked": sum(graph.covers for graph in graphs),
        "cover_size_histogram": dict(
            sorted(aggregate_cover_sizes.items(), key=lambda item: int(item[0]))
        ),
        "individual_cover_failures": dict(aggregate_failures),
        "sequential_marginal_cover_failures": dict(aggregate_marginal),
        "graph_rejections_by_stratum": {
            name: dict(counts) for name, counts in sorted(by_stratum.items())
        },
        "individual_graph_decisions": {
            "old_cover_rejected": sorted(old_indices),
            "cap7_to6_rejected": sorted(cap7_to6_indices),
            "cap6_to5_rejected": sorted(cap6_to5_indices),
            "cap5_to4_rejected": sorted(cap5_to4_indices),
            "cap4_to3_rejected": sorted(cap4_to3_indices),
            "support_only_rejected": sorted(support_indices),
            "subspace_K_only_rejected": sorted(subspace_indices),
            "ordinary_B_only_rejected": sorted(ordinary_indices),
            "component_B_only_rejected": sorted(component_indices),
            "PD_clique_only_rejected": sorted(pd_clique_indices),
            "perpendicular_degree_only_rejected": sorted(
                perpendicular_degree_indices
            ),
            "basis_kernel_only_rejected": sorted(basis_kernel_indices),
            "saturating_mask_only_rejected": sorted(
                saturating_mask_indices
            ),
            "joint_existential_rejected": sorted(joint_indices),
            "joint_plus_PD_clique_rejected": sorted(
                joint_pd_clique_indices
            ),
            "joint_plus_PD_clique_degree_rejected": sorted(
                joint_pd_clique_degree_indices
            ),
            "joint_plus_PD_clique_degree_basis_rejected": sorted(
                joint_pd_clique_degree_basis_indices
            ),
            "enhanced_joint_existential_rejected": sorted(
                enhanced_joint_indices
            ),
            "survivors": survivor_indices,
            "enhanced_survivors": enhanced_survivor_indices,
        },
        "individual_graph_counts": {
            "old_cover_rejected": len(old_indices),
            "cap7_to6_rejected": len(cap7_to6_indices),
            "cap6_to5_rejected": len(cap6_to5_indices),
            "cap5_to4_rejected": len(cap5_to4_indices),
            "cap4_to3_rejected": len(cap4_to3_indices),
            "support_only_rejected": len(support_indices),
            "subspace_K_only_rejected": len(subspace_indices),
            "ordinary_B_only_rejected": len(ordinary_indices),
            "component_B_only_rejected": len(component_indices),
            "PD_clique_only_rejected": len(pd_clique_indices),
            "perpendicular_degree_only_rejected": len(
                perpendicular_degree_indices
            ),
            "basis_kernel_only_rejected": len(basis_kernel_indices),
            "saturating_mask_only_rejected": len(saturating_mask_indices),
            "joint_existential_rejected": len(joint_indices),
            "joint_plus_PD_clique_rejected": len(
                joint_pd_clique_indices
            ),
            "joint_plus_PD_clique_degree_rejected": len(
                joint_pd_clique_degree_indices
            ),
            "joint_plus_PD_clique_degree_basis_rejected": len(
                joint_pd_clique_degree_basis_indices
            ),
            "enhanced_joint_existential_rejected": len(
                enhanced_joint_indices
            ),
            "survivors": len(survivor_indices),
            "enhanced_survivors": len(enhanced_survivor_indices),
        },
        "sequential_marginal_graph_decisions": {
            "cap7_to6_after_old_cover": sorted(marginal_cap7_to6),
            "cap6_to5_after_cap6": sorted(marginal_cap6_to5),
            "cap5_to4_after_cap5": sorted(marginal_cap5_to4),
            "cap4_to3_after_cap4": sorted(marginal_cap4_to3),
            "support_after_cap3": sorted(marginal_support),
            "subspace_K_after_support": sorted(marginal_subspace),
            "ordinary_B_after_support_and_subspace": sorted(marginal_ordinary),
            "component_B_after_support_subspace_and_ordinary": sorted(
                marginal_component
            ),
            "joint_existential_synergy_after_individual_rules": sorted(
                marginal_synergy
            ),
            "PD_clique_existential_after_current_joint": sorted(
                existential_pd_clique
            ),
            "perpendicular_degree_existential_after_PD_clique": sorted(
                existential_perpendicular_degree
            ),
            "basis_kernel_existential_after_degree": sorted(
                existential_basis_kernel
            ),
            "saturating_masks_existential_after_basis_kernel": sorted(
                existential_saturating_masks
            ),
        },
        "sequential_marginal_graph_counts": {
            "cap7_to6_after_old_cover": len(marginal_cap7_to6),
            "cap6_to5_after_cap6": len(marginal_cap6_to5),
            "cap5_to4_after_cap5": len(marginal_cap5_to4),
            "cap4_to3_after_cap4": len(marginal_cap4_to3),
            "support_after_cap3": len(marginal_support),
            "subspace_K_after_support": len(marginal_subspace),
            "ordinary_B_after_support_and_subspace": len(marginal_ordinary),
            "component_B_after_support_subspace_and_ordinary": len(
                marginal_component
            ),
            "joint_existential_synergy_after_individual_rules": len(
                marginal_synergy
            ),
            "PD_clique_existential_after_current_joint": len(
                existential_pd_clique
            ),
            "perpendicular_degree_existential_after_PD_clique": len(
                existential_perpendicular_degree
            ),
            "basis_kernel_existential_after_degree": len(
                existential_basis_kernel
            ),
            "saturating_masks_existential_after_basis_kernel": len(
                existential_saturating_masks
            ),
        },
        "new_rule_individual_marginal_graph_decisions": {
            "PD_clique_standalone_after_current_joint": sorted(
                marginal_pd_clique
            ),
            "perpendicular_degree_standalone_after_prior": sorted(
                marginal_perpendicular_degree
            ),
            "basis_kernel_standalone_after_prior": sorted(
                marginal_basis_kernel
            ),
            "saturating_masks_standalone_after_prior": sorted(
                marginal_saturating_masks
            ),
            "enhanced_existential_synergy_beyond_all_standalone_rules": sorted(
                enhanced_synergy
            ),
        },
        "new_rule_individual_marginal_graph_counts": {
            "PD_clique_standalone_after_current_joint": len(
                marginal_pd_clique
            ),
            "perpendicular_degree_standalone_after_prior": len(
                marginal_perpendicular_degree
            ),
            "basis_kernel_standalone_after_prior": len(
                marginal_basis_kernel
            ),
            "saturating_masks_standalone_after_prior": len(
                marginal_saturating_masks
            ),
            "enhanced_existential_synergy_beyond_all_standalone_rules": len(
                enhanced_synergy
            ),
        },
        "component_B_tier_graph_decisions": {
            name: {
                "standalone_rejected": sorted(component_tier_indices[name]),
                "marginal_after_prior_individual_rules": sorted(
                    component_tier_marginal[name]
                ),
            }
            for name in ("trivial", "tier1", "tier2", "tier3", "exact")
        },
        "component_B_tier_graph_counts": {
            name: {
                "standalone_rejected": len(component_tier_indices[name]),
                "marginal_after_prior_individual_rules": len(
                    component_tier_marginal[name]
                ),
                "fraction_of_exact_marginal": (
                    len(component_tier_marginal[name])
                    / len(component_tier_marginal["exact"])
                    if component_tier_marginal["exact"] else 1.0
                ),
            }
            for name in ("trivial", "tier1", "tier2", "tier3", "exact")
        },
        "runtime": {
            "wall_seconds": elapsed,
            "sum_graph_seconds": sum(graph.elapsed_seconds for graph in graphs),
            "maximum_graph_seconds": max(
                (graph.elapsed_seconds for graph in graphs), default=0.0
            ),
            **dict(counter_totals),
            "saturating_cliques_logically_checked": sum(
                graph.saturating_cliques_checked for graph in graphs
            ),
            "basis_kernel_clique_failures": sum(
                graph.basis_kernel_clique_failures for graph in graphs
            ),
            "saturating_mask_clique_failures": sum(
                graph.saturating_mask_clique_failures for graph in graphs
            ),
            "component_tier_forcing_subset_checks_upper_bound": dict(
                tier_check_totals
            ),
        },
    }
    if include_per_graph:
        report["per_graph"] = [asdict(graph) for graph in graphs]
    return report


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--compact",
        action="store_true",
        help="omit verbose per-graph diagnostics (indices remain in summaries)",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    with args.sample.open(encoding="utf-8") as stream:
        sample = json.load(stream)
    report = sample_report(
        sample, workers=args.workers, include_per_graph=not args.compact
    )
    report["input_sample"] = {
        "file": args.sample.name,
        "sha256": file_sha256(args.sample),
    }
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
        print(f"wrote {args.output}")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
