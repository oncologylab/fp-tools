"""Maximum-cardinality matching within a fixed squared-distance caliper."""

import numpy as np
from scipy.optimize import linear_sum_assignment


def caliper_assignment(squared_distance, maximum_distance):
    """Maximize allowed pairs, then minimize their total squared distance.

    Dummy columns allow each row to remain unmatched. The unmatched penalty
    exceeds the sum of every possible normalized real-edge cost, so no distance
    improvement can compensate for losing an eligible pair.
    """
    distance = np.asarray(squared_distance, dtype=float)
    if distance.ndim != 2 or not np.isfinite(distance).all() or np.any(distance < 0):
        raise ValueError("distances must be a finite nonnegative matrix")
    if not np.isfinite(maximum_distance) or maximum_distance <= 0:
        raise ValueError("caliper must be positive and finite")
    n, m = distance.shape
    if not n or not m:
        return np.array([], dtype=int), np.array([], dtype=int)
    limit = maximum_distance ** 2
    allowed = distance <= limit
    penalty = float(n + 1)
    cost = np.full((n, m + n), penalty)
    cost[:, :m] = np.where(allowed, distance / limit, 2 * penalty)
    rows, columns = linear_sum_assignment(cost)
    real = columns < m
    rows, columns = rows[real], columns[real]
    if not allowed[rows, columns].all():
        raise RuntimeError("assignment contains a forbidden edge")
    return rows, columns
