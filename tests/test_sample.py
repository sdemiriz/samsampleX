"""Tests for sample.py: hashing, ratio computation, ratio lookup helpers."""

import numpy as np
import pytest

from samsamplex.depth import DepthArray
from samsamplex.sample import (
    _compute_ratios,
    _get_max_ratio,
    _get_mean_ratio,
    _get_median_ratio,
    _get_min_ratio,
    _xxh32_fraction,
)


def _make(depths, start=0):
    d = np.array(depths, dtype=np.int32)
    return DepthArray(contig="chr1", start=start, end=start + len(d), depths=d)


# ── _xxh32_fraction ──────────────────────────────────────────────────────────


class TestXxh32Fraction:
    def test_in_range(self):
        f = _xxh32_fraction("read1", seed=42)
        assert 0.0 <= f < 1.0

    def test_deterministic(self):
        assert _xxh32_fraction("readA", 42) == _xxh32_fraction("readA", 42)

    def test_different_names_differ(self):
        assert _xxh32_fraction("readA", 42) != _xxh32_fraction("readB", 42)

    def test_different_seeds_differ(self):
        assert _xxh32_fraction("readA", 1) != _xxh32_fraction("readA", 2)

    def test_many_values_in_range(self):
        fractions = [_xxh32_fraction(f"read_{i}", seed=0) for i in range(1000)]
        assert all(0.0 <= f < 1.0 for f in fractions)


# ── _compute_ratios ──────────────────────────────────────────────────────────


class TestComputeRatios:
    def test_equal_depths(self):
        t = _make([10, 20, 30])
        s = _make([10, 20, 30])
        ratios = _compute_ratios(t, s)
        np.testing.assert_array_almost_equal(ratios, [1.0, 1.0, 1.0])

    def test_template_lower(self):
        t = _make([5, 10, 15])
        s = _make([10, 20, 30])
        ratios = _compute_ratios(t, s)
        np.testing.assert_array_almost_equal(ratios, [0.5, 0.5, 0.5])

    def test_template_higher_capped(self):
        t = _make([100, 100, 100])
        s = _make([10, 20, 30])
        ratios = _compute_ratios(t, s)
        np.testing.assert_array_almost_equal(ratios, [1.0, 1.0, 1.0])

    def test_zero_source(self):
        t = _make([10, 20, 30])
        s = _make([0, 0, 0])
        ratios = _compute_ratios(t, s)
        np.testing.assert_array_almost_equal(ratios, [0.0, 0.0, 0.0])

    def test_mixed(self):
        t = _make([5, 0, 50])
        s = _make([10, 0, 25])
        ratios = _compute_ratios(t, s)
        np.testing.assert_array_almost_equal(ratios, [0.5, 0.0, 1.0])


# ── _get_mean_ratio ──────────────────────────────────────────────────────────


class TestGetMeanRatio:
    def _setup_cumsum(self, ratios, start=0):
        """Build cumsum from a ratio list, return (cumsum, region_start, region_end)."""
        r = np.array(ratios, dtype=np.float64)
        cumsum = np.concatenate(([0.0], np.cumsum(r)))
        return cumsum, start, start + len(ratios)

    def test_full_span(self):
        cumsum, rs, re = self._setup_cumsum([0.5, 0.5, 0.5])
        assert _get_mean_ratio(cumsum, rs, re, rs, re) == pytest.approx(0.5)

    def test_partial_span(self):
        cumsum, rs, re = self._setup_cumsum([0.2, 0.4, 0.6, 0.8])
        # Read covers positions 1-3 (indices 1,2 → values 0.4, 0.6)
        assert _get_mean_ratio(cumsum, rs, re, 1, 3) == pytest.approx(0.5)

    def test_no_overlap_returns_zero(self):
        cumsum, rs, re = self._setup_cumsum([0.5, 0.5])
        assert _get_mean_ratio(cumsum, rs, re, 10, 20) == pytest.approx(0.0)

    def test_read_extends_beyond_region(self):
        """Read that extends past region end should be clipped."""
        cumsum, rs, re = self._setup_cumsum([1.0, 1.0, 1.0])  # region 0-3
        assert _get_mean_ratio(cumsum, rs, re, 1, 100) == pytest.approx(1.0)

    def test_single_position(self):
        cumsum, rs, re = self._setup_cumsum([0.3, 0.7, 0.5])
        assert _get_mean_ratio(cumsum, rs, re, 1, 2) == pytest.approx(0.7)


# ── _get_min_ratio ───────────────────────────────────────────────────────────


class TestGetMinRatio:
    def test_uniform(self):
        ratios = np.array([0.5, 0.5, 0.5])
        assert _get_min_ratio(ratios, 0, 3, 0, 3) == pytest.approx(0.5)

    def test_varying(self):
        ratios = np.array([0.1, 0.9, 0.3, 0.7])
        assert _get_min_ratio(ratios, 0, 4, 0, 4) == pytest.approx(0.1)

    def test_partial(self):
        ratios = np.array([0.1, 0.9, 0.3, 0.7])
        assert _get_min_ratio(ratios, 0, 4, 1, 3) == pytest.approx(0.3)

    def test_no_overlap(self):
        ratios = np.array([0.5, 0.5])
        assert _get_min_ratio(ratios, 0, 2, 5, 10) == pytest.approx(0.0)


# ── _get_max_ratio ───────────────────────────────────────────────────────────


class TestGetMaxRatio:
    def test_uniform(self):
        ratios = np.array([0.5, 0.5, 0.5])
        assert _get_max_ratio(ratios, 0, 3, 0, 3) == pytest.approx(0.5)

    def test_varying(self):
        ratios = np.array([0.1, 0.9, 0.3, 0.7])
        assert _get_max_ratio(ratios, 0, 4, 0, 4) == pytest.approx(0.9)

    def test_partial(self):
        ratios = np.array([0.1, 0.9, 0.3, 0.7])
        assert _get_max_ratio(ratios, 0, 4, 2, 4) == pytest.approx(0.7)

    def test_no_overlap(self):
        ratios = np.array([0.5, 0.5])
        assert _get_max_ratio(ratios, 0, 2, 5, 10) == pytest.approx(0.0)


# ── _get_median_ratio ────────────────────────────────────────────────────────


class TestGetMedianRatio:
    def test_odd_count(self):
        ratios = np.array([0.1, 0.5, 0.9])
        assert _get_median_ratio(ratios, 0, 3, 0, 3) == pytest.approx(0.5)

    def test_even_count(self):
        ratios = np.array([0.1, 0.3, 0.7, 0.9])
        # median of [0.1, 0.3, 0.7, 0.9] = (0.3 + 0.7) / 2 = 0.5
        assert _get_median_ratio(ratios, 0, 4, 0, 4) == pytest.approx(0.5)

    def test_single_position(self):
        ratios = np.array([0.2, 0.8, 0.4])
        assert _get_median_ratio(ratios, 0, 3, 1, 2) == pytest.approx(0.8)

    def test_no_overlap(self):
        ratios = np.array([0.5, 0.5])
        assert _get_median_ratio(ratios, 0, 2, 5, 10) == pytest.approx(0.0)
