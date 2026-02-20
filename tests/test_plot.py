"""Tests for plot.py: _pick_step, _downsample, write_tsv."""

import io

import numpy as np
import pytest

from samsamplex.depth import DepthArray
from samsamplex.plot import _downsample, _pick_step, write_tsv


# ── _pick_step ───────────────────────────────────────────────────────────────


class TestPickStep:
    def test_small_region(self):
        assert _pick_step(100) == 1

    def test_boundary_10k(self):
        assert _pick_step(10_000) == 1

    def test_just_above_10k(self):
        assert _pick_step(10_001) == 10

    def test_100k(self):
        assert _pick_step(100_000) == 10

    def test_1m(self):
        assert _pick_step(1_000_000) == 100

    def test_10m(self):
        assert _pick_step(10_000_000) == 1000

    def test_above_10m(self):
        assert _pick_step(50_000_000) == 5000

    def test_monotonic_increase(self):
        """Larger regions should always get equal or larger step sizes."""
        sizes = [1, 100, 10_000, 10_001, 100_000, 1_000_000, 10_000_000, 100_000_000]
        steps = [_pick_step(s) for s in sizes]
        for i in range(1, len(steps)):
            assert steps[i] >= steps[i - 1]


# ── _downsample ──────────────────────────────────────────────────────────────


class TestDownsample:
    def test_step_1_returns_all(self):
        pos = np.arange(5)
        vals = np.array([10, 20, 30, 40, 50])
        ds_pos, ds_vals = _downsample(pos, vals, step=1)
        np.testing.assert_array_equal(ds_pos, pos)
        np.testing.assert_array_equal(ds_vals, vals)

    def test_step_2(self):
        pos = np.arange(6)
        vals = np.array([10, 20, 30, 40, 50, 60])
        ds_pos, ds_vals = _downsample(pos, vals, step=2)
        # indices 0, 2, 4 — last element (5) is included too
        np.testing.assert_array_equal(ds_pos, [0, 2, 4, 5])
        np.testing.assert_array_equal(ds_vals, [10, 30, 50, 60])

    def test_last_point_always_included(self):
        pos = np.arange(10)
        vals = np.arange(10)
        ds_pos, _ = _downsample(pos, vals, step=3)
        assert ds_pos[-1] == 9

    def test_step_larger_than_array(self):
        pos = np.arange(3)
        vals = np.array([1, 2, 3])
        ds_pos, ds_vals = _downsample(pos, vals, step=100)
        # only index 0 from arange, plus last (2)
        np.testing.assert_array_equal(ds_pos, [0, 2])
        np.testing.assert_array_equal(ds_vals, [1, 3])

    def test_multiple_arrays(self):
        pos = np.arange(6)
        a = np.array([1, 2, 3, 4, 5, 6])
        b = np.array([10, 20, 30, 40, 50, 60])
        ds_pos, ds_a, ds_b = _downsample(pos, a, b, step=3)
        np.testing.assert_array_equal(ds_pos, [0, 3, 5])
        np.testing.assert_array_equal(ds_a, [1, 4, 6])
        np.testing.assert_array_equal(ds_b, [10, 40, 60])


# ── write_tsv ────────────────────────────────────────────────────────────────


class TestWriteTsv:
    def _make(self, depths, start=0):
        d = np.array(depths, dtype=np.int32)
        return DepthArray(contig="chr1", start=start, end=start + len(d), depths=d)

    def test_header_line(self):
        src = self._make([10, 20])
        tpl = self._make([5, 10])
        out = self._make([8, 15])
        buf = io.StringIO()
        write_tsv(buf, 0, src, tpl, out, step=1)
        lines = buf.getvalue().strip().split("\n")
        assert lines[0] == "position\tsource_depth\ttemplate_depth\toutput_depth"

    def test_correct_positions_1based(self):
        """Positions in TSV should be 1-based."""
        src = self._make([10, 20, 30], start=100)
        tpl = self._make([5, 10, 15], start=100)
        out = self._make([8, 15, 25], start=100)
        buf = io.StringIO()
        write_tsv(buf, 100, src, tpl, out, step=1)
        lines = buf.getvalue().strip().split("\n")
        # First data line should start at position 101 (100 + 0 + 1)
        assert lines[1].startswith("101\t")
        assert lines[3].startswith("103\t")

    def test_data_values(self):
        src = self._make([10, 20])
        tpl = self._make([5, 10])
        out = self._make([8, 15])
        buf = io.StringIO()
        write_tsv(buf, 0, src, tpl, out, step=1)
        lines = buf.getvalue().strip().split("\n")
        assert lines[1] == "1\t10\t5\t8"
        assert lines[2] == "2\t20\t10\t15"

    def test_downsampled_output(self):
        src = self._make([1, 2, 3, 4, 5, 6])
        tpl = self._make([10, 20, 30, 40, 50, 60])
        out = self._make([0, 0, 0, 0, 0, 0])
        buf = io.StringIO()
        write_tsv(buf, 0, src, tpl, out, step=3)
        lines = buf.getvalue().strip().split("\n")
        # header + downsampled data rows (indices 0, 3, 5)
        assert len(lines) == 4  # header + 3 data rows
