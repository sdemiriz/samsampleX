"""Quality metrics for comparing depth distributions."""

from __future__ import annotations

import sys
from collections import Counter
from dataclasses import dataclass

import numpy as np

from .depth import DepthArray


@dataclass
class MetricsResult:
    wasserstein: float
    tv: float


def _wasserstein(a: np.ndarray, b: np.ndarray) -> float:
    """Normalised Wasserstein-1 distance between two 1-D depth distributions.

    W1 is computed as the mean absolute difference between the empirical CDFs
    of the two depth arrays.
    """
    n = len(a)
    if n == 0:
        return 0.0

    cum_a = np.cumsum(a, dtype=np.float64)
    cum_b = np.cumsum(b, dtype=np.float64)

    total_a = cum_a[-1]
    total_b = cum_b[-1]

    # Prepend 0 so the arrays have length n+1
    cum_a = np.concatenate(([0.0], cum_a))
    cum_b = np.concatenate(([0.0], cum_b))

    if total_a > 0:
        cum_a /= total_a
    if total_b > 0:
        cum_b /= total_b

    return float(np.sum(np.abs(cum_a - cum_b)) / (n + 1))


def _total_variation(a: np.ndarray, b: np.ndarray) -> float:
    """Total Variation distance between empirical distributions of *a* and *b*.

    For PMFs p(x) = count(x)/n over each multiset, TV = (1/2) * sum_x |p_a(x) - p_b(x)|.
    """
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return 0.0
    c1, c2 = Counter(a), Counter(b)
    all_vals = c1.keys() | c2.keys()
    return float(
        sum(abs(c1.get(x, 0) / n1 - c2.get(x, 0) / n2) for x in all_vals) / 2
    )


def metrics_calculate(depth_a: DepthArray, depth_b: DepthArray) -> MetricsResult:
    """Calculate comparison metrics between two DepthArrays."""
    if depth_a.length != depth_b.length:
        raise ValueError(
            f"Depth arrays have different lengths ({depth_a.length} vs {depth_b.length})"
        )

    return MetricsResult(
        wasserstein=_wasserstein(depth_a.depths, depth_b.depths),
        tv=_total_variation(depth_a.depths, depth_b.depths),
    )


def metrics_print(result: MetricsResult, label_a: str = "A", label_b: str = "B") -> None:
    """Print metrics to stderr."""
    print(f"Total Variation:         {result.tv:.4f}", file=sys.stderr)
    print(f"Norm. Wasserstein Dist.: {result.wasserstein:.6f}", file=sys.stderr)
