#!/usr/bin/env python3
"""Exact clique-link caps for propagated K7 defect masks.

For a seed-coordinate set ``S``, every outside point whose exact support is
contained in ``S`` is a common unit neighbour of the complementary seed
clique.  A propagated mask ``E_x`` is an upper bound on exact support, so
``E_x subseteq S`` is already enough to count the point.  Fixed labeled
zero-factor supports can be counted at the same time.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Sequence


# Link bounds minus the k seed vertices in S that are themselves common unit
# neighbours of the complementary K_(7-k).
OUTSIDE_CAPS = {1: 1, 2: 2, 3: 7, 4: 8, 5: 11}
SUBSETS_BY_SIZE = {
    size: tuple(sum(1 << coordinate for coordinate in chosen)
                for chosen in combinations(range(7), size))
    for size in OUTSIDE_CAPS
}


@dataclass(frozen=True)
class LinkCapFailure:
    subset_mask: int
    subset_size: int
    outside_cap: int
    forced_count: int
    zero_factor_count: int
    nonzero_factor_count: int


def first_link_cap_failure(
    zero_factor_supports: Sequence[int],
    propagated_nonzero_masks: Sequence[int],
) -> LinkCapFailure | None:
    """Return the first exact propagated link-cap obstruction, if any.

    ``zero_factor_supports`` are exact actual supports.  Each propagated N
    mask merely contains its unknown actual support; only masks wholly inside
    the tested seed subset are counted.  Ignoring every other point preserves
    the sound direction.
    """

    z_supports = tuple(map(int, zero_factor_supports))
    n_masks = tuple(map(int, propagated_nonzero_masks))
    if any(mask <= 0 or mask >= 128 for mask in (*z_supports, *n_masks)):
        raise ValueError("supports and propagated masks must be nonempty 7-bit masks")

    for size, outside_cap in OUTSIDE_CAPS.items():
        for subset in SUBSETS_BY_SIZE[size]:
            z_count = sum(not (support & ~subset) for support in z_supports)
            n_count = sum(not (mask & ~subset) for mask in n_masks)
            total = z_count + n_count
            if total > outside_cap:
                return LinkCapFailure(
                    subset_mask=subset,
                    subset_size=size,
                    outside_cap=outside_cap,
                    forced_count=total,
                    zero_factor_count=z_count,
                    nonzero_factor_count=n_count,
                )
    return None

