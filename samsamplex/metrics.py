"""Quality metrics for comparing depth distributions."""

from __future__ import annotations

import sys
from dataclasses import dataclass

import numpy as np
from scipy.stats import wasserstein_distance

from .depth import DepthArray


def depth_from_path(path: str, contig: str, start: int, end: int) -> DepthArray:
    """Load a DepthArray from *path*, auto-detecting format by extension.

    .bed files are read via bed_read_depths; everything else is treated as BAM.
    """
    if path.endswith(".bed"):
        from .bed import bed_read_depths
        return bed_read_depths(path, contig, start, end)
    from .depth import depth_from_bam
    return depth_from_bam(path, contig, start, end)


@dataclass
class MetricsResult:
    mae: float
    wasserstein: float


def _wasserstein(a: np.ndarray, b: np.ndarray) -> float:
    """
    Wasserstein-1 style L1 distance between empirical CDFs of two depth arrays.
    """
    n = min(len(a), len(b))
    assert n > 0, "Wasserstein distance requires non-empty arrays"

    return float(wasserstein_distance(a, b))


def _mean_absolute_error(a: np.ndarray, b: np.ndarray) -> float:
    """Mean absolute per-base depth difference between aligned arrays *a* and *b*."""
    assert len(a) > 0 and len(b) > 0, "MAE requires non-empty arrays"
    return float(np.mean(np.abs(a.astype(np.float64) - b.astype(np.float64))))


def metrics_calculate(depth_a: DepthArray, depth_b: DepthArray) -> MetricsResult:
    """Calculate comparison metrics between two DepthArrays."""
    if depth_a.length != depth_b.length:
        raise ValueError(
            f"Depth arrays have different lengths ({depth_a.length} vs {depth_b.length})"
        )

    return MetricsResult(
        wasserstein=_wasserstein(depth_a.depths, depth_b.depths),
        mae=_mean_absolute_error(depth_a.depths, depth_b.depths),
    )


def metrics_print(result: MetricsResult, label_a: str = "A", label_b: str = "B") -> None:
    """Print metrics to stderr."""
    print(f"Mean Absolute Error:     {result.mae:.4f}", file=sys.stderr)
    print(f"Wasserstein-1 distance: {result.wasserstein:.6f}", file=sys.stderr)
