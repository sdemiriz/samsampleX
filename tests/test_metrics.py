"""Tests for metrics.py: windowed stats, top-by-|value| selection, depth helpers."""

import contextlib
import io
import math

import numpy as np
import pytest

from samsamplex.depth import DepthArray
from samsamplex.metrics import (
    STATS_HEADER,
    WindowMetricsRow,
    _fmt_stderr_metric,
    _mean_absolute_error,
    _mean_signed_depth_difference,
    any_zero_depth_window,
    metrics_window_rows,
    print_stats_stderr,
    top_k_window_indices_by_abs,
    write_stats_tsv,
)


def _make(depths, *, start=0, chrom="chr1"):
    d = np.array(depths, dtype=np.int32)
    return DepthArray(contig=chrom, start=start, end=start + len(d), depths=d)


def _row(
    chrom: str,
    start: int,
    end: int,
    mean_temp: float,
    depth_diff: float,
    rel_depth: float,
    mean_res: float,
) -> WindowMetricsRow:
    return WindowMetricsRow(
        chrom=chrom,
        start=start,
        end=end,
        mean_depth_temp=mean_temp,
        depth_diff=depth_diff,
        rel_depth_diff=rel_depth,
        mean_depth_res=mean_res,
    )


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


# ── _mean_signed_depth_difference ─────────────────────────────────────────────


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


# ── metrics_window_rows ───────────────────────────────────────────────────────


class TestMetricsWindowRows:
    def test_uniform_shift(self):
        temp = _make([10, 10, 10, 10])
        res = _make([20, 20, 20, 20])
        rows = metrics_window_rows(temp, res, window=4)
        assert len(rows) == 1
        r = rows[0]
        assert r.depth_diff == pytest.approx(10.0)
        assert r.mean_depth_temp == pytest.approx(10.0)
        assert r.mean_depth_res == pytest.approx(20.0)
        assert r.rel_depth_diff == pytest.approx(1.0)

    def test_partial_last_window(self):
        temp = _make([10] * 10)
        res = _make([12] * 10)
        rows = metrics_window_rows(temp, res, window=4)
        assert len(rows) == 3
        assert rows[0].start == 0 and rows[0].end == 4
        assert rows[1].start == 4 and rows[1].end == 8
        assert rows[2].start == 8 and rows[2].end == 10
        assert rows[2].depth_diff == pytest.approx(2.0)

    def test_zero_template_relative_nan(self):
        temp = _make([0, 0])
        res = _make([5, 5])
        rows = metrics_window_rows(temp, res, window=2)
        assert math.isnan(rows[0].rel_depth_diff)
        assert rows[0].depth_diff == pytest.approx(5.0)

    def test_length_mismatch_raises(self):
        a = _make([1, 2])
        b = _make([1, 2, 3])
        with pytest.raises(ValueError, match="different lengths"):
            metrics_window_rows(a, b, window=1)

    def test_window_bad_raises(self):
        a = _make([1])
        b = _make([1])
        with pytest.raises(ValueError, match="window"):
            metrics_window_rows(a, b, window=0)


# ── top_k_window_indices_by_abs ───────────────────────────────────────────────


class TestTopKWindowIndicesByAbs:
    def test_top_ten_of_many(self):
        """|i-50| is largest near endpoints; tie-break by start."""
        rows = [
            _row("chr1", i, i + 1, 10.0, float(i - 50), 0.0, 10.0) for i in range(100)
        ]
        idx = top_k_window_indices_by_abs(rows, lambda r: r.depth_diff, k=10)
        assert len(idx) == 10
        assert idx == [0, 1, 99, 2, 98, 3, 97, 4, 96, 5]

    def test_all_when_fewer_than_k(self):
        rows = [_row("chr1", i, i + 1, 1.0, float(i), 0.0, 1.0) for i in range(3)]
        idx = top_k_window_indices_by_abs(rows, lambda r: r.depth_diff, k=10)
        assert len(idx) == 3
        assert idx == [2, 1, 0]

    def test_ties_same_abs_prefers_lower_start(self):
        rows = [
            _row("chr1", 0, 1, 1.0, 5.0, 0.0, 1.0),
            _row("chr1", 1, 2, 1.0, -5.0, 0.0, 1.0),
            _row("chr1", 2, 3, 1.0, 3.0, 0.0, 1.0),
        ]
        idx = top_k_window_indices_by_abs(rows, lambda r: r.depth_diff, k=2)
        assert idx == [0, 1]

    def test_excludes_nan(self):
        rows = [
            _row("chr1", 0, 1, 0.0, 1.0, float("nan"), 1.0),
            _row("chr1", 1, 2, 1.0, 1.0, 0.5, 1.0),
        ]
        idx = top_k_window_indices_by_abs(rows, lambda r: r.rel_depth_diff, k=10)
        assert idx == [1]

    def test_empty_rows(self):
        assert top_k_window_indices_by_abs([], lambda r: r.depth_diff, k=10) == []


# ── _fmt_stderr_metric ─────────────────────────────────────────────────────────


class TestFmtStderrMetric:
    def test_signed_plus_minus(self):
        assert _fmt_stderr_metric(1.5, signed=True).startswith("+")
        assert _fmt_stderr_metric(-2.0, signed=True).startswith("-")
        assert not _fmt_stderr_metric(3.0, signed=False).startswith("+")

    def test_unsigned_no_plus(self):
        s = _fmt_stderr_metric(3.0, signed=False)
        assert not s.startswith("+")


# ── any_zero_depth_window / write_stats_tsv ───────────────────────────────────


class TestHelpers:
    def test_any_zero_depth(self):
        rows = [_row("chr1", 0, 1, 1.0, 0.0, 0.0, 1.0)]
        assert not any_zero_depth_window(rows)
        rows.append(_row("chr1", 1, 2, 0.0, 0.0, float("nan"), 1.0))
        assert any_zero_depth_window(rows)

    def test_write_stats_tsv_header(self):
        rows = [_row("chr1", 10, 20, 10.0, 1.0, 0.1, 11.0)]
        buf = io.StringIO()
        write_stats_tsv(buf, rows)
        out = buf.getvalue()
        assert out.startswith(STATS_HEADER)
        assert "chr1\t10\t20\t" in out


class TestPrintStatsStderr:
    def test_stderr_rows_limits_each_section(self):
        rows = [
            _row("chr1", i, i + 1, 10.0, float(i), float(i) / 10.0, 10.0) for i in range(20)
        ]
        buf = io.StringIO()
        with contextlib.redirect_stderr(buf):
            print_stats_stderr(
                ["prog", "stats"],
                rows,
                warn_zero_depth=False,
                stderr_rows=2,
            )
        text = buf.getvalue()
        assert text.count("Top 2 windows") == 2
