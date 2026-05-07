"""Tests for metrics.py: Wasserstein-1, mean absolute error, metrics_calculate."""

import numpy as np
import pytest

from samsamplex.depth import DepthArray
from samsamplex.metrics import (
    _mean_absolute_error,
    _mean_signed_depth_difference,
    _wasserstein,
    metrics_calculate,
)


def _make(depths):
    d = np.array(depths, dtype=np.int32)
    return DepthArray(contig="chr1", start=0, end=len(d), depths=d)


# ── _wasserstein ─────────────────────────────────────────────────────────────


class TestWasserstein:
    def test_identical_is_zero(self):
        a = np.array([10, 20, 30], dtype=np.int32)
        assert _wasserstein(a, a) == pytest.approx(0.0)

    def test_empty_is_zero(self):
        a = np.array([], dtype=np.int32)
        assert _wasserstein(a, a) == 0.0

    def test_positive_for_different(self):
        a = np.array([10, 10, 10], dtype=np.int32)
        b = np.array([0, 0, 30], dtype=np.int32)
        assert _wasserstein(a, b) > 0.0

    def test_symmetric(self):
        a = np.array([10, 20, 30], dtype=np.int32)
        b = np.array([30, 20, 10], dtype=np.int32)
        assert _wasserstein(a, b) == pytest.approx(_wasserstein(b, a))

    def test_all_zeros(self):
        a = np.zeros(5, dtype=np.int32)
        assert _wasserstein(a, a) == pytest.approx(0.0)


# ── _mean_absolute_error ─────────────────────────────────────────────────────


class TestMeanAbsoluteError:
    def test_identical_is_zero(self):
        a = np.array([10, 20, 30], dtype=np.int32)
        assert _mean_absolute_error(a, a) == pytest.approx(0.0)

    def test_empty_is_zero(self):
        a = np.array([], dtype=np.int32)
        assert _mean_absolute_error(a, a) == 0.0

    def test_constant_shift(self):
        a = np.array([10, 10, 10], dtype=np.int32)
        b = np.array([20, 20, 20], dtype=np.int32)
        assert _mean_absolute_error(a, b) == pytest.approx(10.0)

    def test_symmetric(self):
        a = np.array([5, 15, 25], dtype=np.int32)
        b = np.array([25, 15, 5], dtype=np.int32)
        assert _mean_absolute_error(a, b) == pytest.approx(_mean_absolute_error(b, a))

    def test_mixed_positions(self):
        a = np.array([10, 10, 10, 10], dtype=np.int32)
        b = np.array([10, 10, 10, 20], dtype=np.int32)
        assert _mean_absolute_error(a, b) == pytest.approx(2.5)


# ── _mean_signed_depth_difference ───────────────────────────────────────────


class TestMeanSignedDepthDifference:
    def test_identical_is_zero(self):
        a = np.array([10, 20, 30], dtype=np.int32)
        assert _mean_signed_depth_difference(a, a) == pytest.approx(0.0)

    def test_empty_is_zero(self):
        a = np.array([], dtype=np.int32)
        assert _mean_signed_depth_difference(a, a) == 0.0

    def test_b_deeper_positive(self):
        a = np.array([10, 10, 10], dtype=np.int32)
        b = np.array([20, 20, 20], dtype=np.int32)
        assert _mean_signed_depth_difference(a, b) == pytest.approx(10.0)

    def test_a_deeper_negative(self):
        a = np.array([20, 20, 20], dtype=np.int32)
        b = np.array([10, 10, 10], dtype=np.int32)
        assert _mean_signed_depth_difference(a, b) == pytest.approx(-10.0)


# ── metrics_calculate ────────────────────────────────────────────────────────


class TestMetricsCalculate:
    def test_identical_arrays(self):
        a = _make([10, 20, 30])
        result = metrics_calculate(a, a)
        assert result.wasserstein == pytest.approx(0.0)
        assert result.mae == pytest.approx(0.0)
        assert result.mean_signed_depth_diff == pytest.approx(0.0)

    def test_different_arrays(self):
        a = _make([10, 10, 10])
        b = _make([20, 20, 20])
        result = metrics_calculate(a, b)
        assert result.mae == pytest.approx(10.0)
        assert result.mean_signed_depth_diff == pytest.approx(10.0)
        assert result.wasserstein >= 0.0

    def test_length_mismatch_raises(self):
        a = _make([10, 20])
        b = _make([10, 20, 30])
        with pytest.raises(ValueError, match="different lengths"):
            metrics_calculate(a, b)
