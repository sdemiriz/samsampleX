"""Tests for bed.py: contig matching, BED I/O, depth combining."""

import io

import numpy as np
import pytest

from samsamplex.bed import (
    bed_combine_depths,
    bed_read_depths,
    contig_names_match,
    write_bed_entry,
    write_bed_output,
)
from samsamplex.depth import DepthArray


# ── contig_names_match ───────────────────────────────────────────────────────


class TestContigNamesMatch:
    def test_exact(self):
        assert contig_names_match("chr1", "chr1") is True

    def test_chr_prefix_left(self):
        assert contig_names_match("chr1", "1") is True

    def test_chr_prefix_right(self):
        assert contig_names_match("1", "chr1") is True

    def test_no_match(self):
        assert contig_names_match("chr1", "chr2") is False

    def test_both_bare(self):
        assert contig_names_match("X", "X") is True

    def test_chrX_vs_X(self):
        assert contig_names_match("chrX", "X") is True


# ── write_bed_entry ──────────────────────────────────────────────────────────


class TestWriteBedEntry:
    def test_single_line(self):
        buf = io.StringIO()
        write_bed_entry(buf, "chr1", 100, 200, 30)
        assert buf.getvalue() == "chr1\t100\t200\t30\n"


# ── write_bed_output ─────────────────────────────────────────────────────────


class TestWriteBedOutput:
    def _make_arr(self, depths, contig="chr1", start=0):
        d = np.array(depths, dtype=np.int32)
        return DepthArray(contig=contig, start=start, end=start + len(d), depths=d)

    def test_per_base_output(self):
        arr = self._make_arr([10, 20, 30], start=100)
        buf = io.StringIO()
        write_bed_output(buf, arr, collapse=0)
        lines = buf.getvalue().strip().split("\n")
        assert len(lines) == 3
        assert lines[0] == "chr1\t100\t101\t10"
        assert lines[1] == "chr1\t101\t102\t20"
        assert lines[2] == "chr1\t102\t103\t30"

    def test_empty_array(self):
        arr = self._make_arr([], start=0)
        buf = io.StringIO()
        write_bed_output(buf, arr, collapse=0)
        assert buf.getvalue() == ""

    def test_collapse_no_change(self):
        """All depths identical → single interval."""
        arr = self._make_arr([10, 10, 10, 10], start=0)
        buf = io.StringIO()
        write_bed_output(buf, arr, collapse=1)
        lines = buf.getvalue().strip().split("\n")
        assert len(lines) == 1
        assert lines[0] == "chr1\t0\t4\t10"

    def test_collapse_splits(self):
        """Depth jump exceeds threshold → two intervals."""
        arr = self._make_arr([10, 10, 50, 50], start=0)
        buf = io.StringIO()
        write_bed_output(buf, arr, collapse=5)
        lines = buf.getvalue().strip().split("\n")
        assert len(lines) == 2
        assert lines[0] == "chr1\t0\t2\t10"
        assert lines[1] == "chr1\t2\t4\t50"

    def test_collapse_within_threshold(self):
        """Small depth changes within threshold → merged."""
        arr = self._make_arr([10, 12, 11, 13], start=0)
        buf = io.StringIO()
        write_bed_output(buf, arr, collapse=5)
        lines = buf.getvalue().strip().split("\n")
        assert len(lines) == 1


# ── bed_read_depths ──────────────────────────────────────────────────────────


class TestBedReadDepths:
    def test_basic_read(self, tmp_path):
        bed = tmp_path / "test.bed"
        bed.write_text("chr1\t0\t3\t10\nchr1\t3\t5\t20\n")

        arr = bed_read_depths(str(bed), "chr1", 0, 5)
        np.testing.assert_array_equal(arr.depths, [10, 10, 10, 20, 20])

    def test_region_subset(self, tmp_path):
        bed = tmp_path / "test.bed"
        bed.write_text("chr1\t0\t10\t50\n")

        arr = bed_read_depths(str(bed), "chr1", 3, 7)
        assert arr.length == 4
        np.testing.assert_array_equal(arr.depths, [50, 50, 50, 50])

    def test_skip_other_contig(self, tmp_path):
        bed = tmp_path / "test.bed"
        bed.write_text("chr2\t0\t5\t99\nchr1\t0\t5\t10\n")

        arr = bed_read_depths(str(bed), "chr1", 0, 5)
        np.testing.assert_array_equal(arr.depths, [10, 10, 10, 10, 10])

    def test_chr_prefix_matching(self, tmp_path):
        """BED has 'chr1' but query uses '1'."""
        bed = tmp_path / "test.bed"
        bed.write_text("chr1\t0\t3\t25\n")

        arr = bed_read_depths(str(bed), "1", 0, 3)
        np.testing.assert_array_equal(arr.depths, [25, 25, 25])

    def test_comments_skipped(self, tmp_path):
        bed = tmp_path / "test.bed"
        bed.write_text("# header\nchr1\t0\t3\t10\n\nchr1\t3\t5\t20\n")

        arr = bed_read_depths(str(bed), "chr1", 0, 5)
        np.testing.assert_array_equal(arr.depths, [10, 10, 10, 20, 20])

    def test_no_overlap_returns_zeros(self, tmp_path):
        bed = tmp_path / "test.bed"
        bed.write_text("chr1\t100\t200\t50\n")

        arr = bed_read_depths(str(bed), "chr1", 0, 10)
        np.testing.assert_array_equal(arr.depths, np.zeros(10, dtype=np.int32))


# ── bed_combine_depths ───────────────────────────────────────────────────────


class TestBedCombineDepths:
    def _make(self, depths):
        d = np.array(depths, dtype=np.int32)
        return DepthArray(contig="chr1", start=0, end=len(d), depths=d)

    def test_single_array_returns_copy(self):
        a = self._make([10, 20, 30])
        result = bed_combine_depths([a])
        np.testing.assert_array_equal(result.depths, a.depths)
        assert result.depths is not a.depths  # must be a copy

    def test_mode_min(self):
        a = self._make([10, 50, 30])
        b = self._make([20, 10, 40])
        result = bed_combine_depths([a, b], mode="min")
        np.testing.assert_array_equal(result.depths, [10, 10, 30])

    def test_mode_max(self):
        a = self._make([10, 50, 30])
        b = self._make([20, 10, 40])
        result = bed_combine_depths([a, b], mode="max")
        np.testing.assert_array_equal(result.depths, [20, 50, 40])

    def test_mode_mean(self):
        a = self._make([10, 50, 30])
        b = self._make([20, 10, 40])
        result = bed_combine_depths([a, b], mode="mean")
        # integer division: (10+20)//2=15, (50+10)//2=30, (30+40)//2=35
        np.testing.assert_array_equal(result.depths, [15, 30, 35])

    def test_mode_random_deterministic(self):
        a = self._make([10, 10, 10])
        b = self._make([20, 20, 20])
        r1 = bed_combine_depths([a, b], mode="random", seed=42)
        r2 = bed_combine_depths([a, b], mode="random", seed=42)
        np.testing.assert_array_equal(r1.depths, r2.depths)

    def test_mode_random_in_range(self):
        a = self._make([10, 10, 10])
        b = self._make([20, 20, 20])
        result = bed_combine_depths([a, b], mode="random", seed=1)
        assert all(10 <= v <= 20 for v in result.depths)

    def test_unknown_mode_raises(self):
        a = self._make([1, 2, 3])
        b = self._make([4, 5, 6])
        with pytest.raises(ValueError, match="Unknown combine mode"):
            bed_combine_depths([a, b], mode="bogus")
