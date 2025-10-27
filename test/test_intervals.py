"""
Tests for Intervals class.
"""

import pytest
import numpy as np
from subsample_reads.Intervals import Intervals


class TestIntervals:
    """Test cases for Intervals class."""

    # Test data file paths - all centralized in test/data/
    TEST_BED = "test/data/test-100bp-10count.bed"
    TEST_BED_VALID = "test/data/test-dataframe-dimensions.bed"
    TEST_BED_MULTI_CONTIG = "test/data/test-contigs-not-unique.bed"
    TEST_BED_OVERLAPPING = "test/data/test-overlapping-intervals.bed"

    def test_file_not_found(self):
        """Test Intervals with nonexistent BED file."""
        with pytest.raises(FileNotFoundError):
            Intervals(bed_file="DOES_NOT_EXIST.bed")

    def test_init_valid_bed_file(self):
        """Test Intervals initialization with valid BED file."""
        intervals = Intervals(bed_file=self.TEST_BED)

        # Verify basic properties are set
        assert intervals.contig is not None
        assert intervals.start is not None
        assert intervals.end is not None
        assert intervals.tree is not None
        assert intervals.bed is not None

    def test_init_with_chr6_bed(self):
        """Test Intervals initialization with chr6 test BED."""
        intervals = Intervals(bed_file=self.TEST_BED_VALID)

        # Verify values from the actual test file
        assert intervals.contig == "chr6"
        assert intervals.start == 25000000
        assert intervals.end == 35000000
        assert len(intervals.tree) == 10

    def test_bed_dataframe_columns(self):
        """Test that BED DataFrame has expected columns."""
        intervals = Intervals(bed_file=self.TEST_BED)

        expected_cols = ["contig", "start", "end", "read_count", "fraction"]
        assert list(intervals.bed.columns) == expected_cols

    def test_non_unique_contigs(self):
        """Test that multiple contigs in BED file raises AssertionError."""
        with pytest.raises(AssertionError) as exc_info:
            Intervals(bed_file=self.TEST_BED_MULTI_CONTIG)

        assert "Not all contig values in BED file identical" in str(exc_info.value)

    def test_overlapping_intervals(self):
        """Test that overlapping intervals raise AssertionError."""
        with pytest.raises(AssertionError) as exc_info:
            Intervals(bed_file=self.TEST_BED_OVERLAPPING)

        assert "BED file contains overlapping intervals" in str(exc_info.value)

    def test_no_bed_provided(self):
        """Test that providing no file raises ValueError."""
        with pytest.raises(ValueError) as exc_info:
            Intervals(bed_file=None)

        assert "No BED file provided" in str(exc_info.value)

    def test_len_method(self):
        """Test __len__ method returns number of intervals."""
        intervals = Intervals(bed_file=self.TEST_BED)

        # __len__ should return the number of intervals in the tree
        assert len(intervals) > 0
        assert len(intervals) == len(intervals.tree)

    def test_get_limits(self):
        """Test get_limits returns correct start and end."""
        intervals = Intervals(bed_file=self.TEST_BED_VALID)

        # Verify start is minimum and end is maximum from BED
        assert intervals.start == intervals.bed["start"].min()
        assert intervals.end == intervals.bed["end"].max()

    def test_tree_content(self):
        """Test that IntervalTree contains correct intervals with data."""
        intervals = Intervals(bed_file=self.TEST_BED)

        # Verify tree has the same number of intervals as BED file
        assert len(intervals.tree) == len(intervals.bed)

        # Verify each interval from BED exists in tree with correct data
        for idx, row in intervals.bed.iterrows():
            # Query tree for intervals overlapping this range
            overlaps = intervals.tree[row["start"] : row["end"]]

            # Should find at least one interval (the one we're looking for)
            assert len(overlaps) > 0

            # Find the exact interval
            found = False
            for interval in overlaps:
                if interval.begin == row["start"] and interval.end == row["end"]:
                    # Verify the data (read_count) is stored correctly
                    assert interval.data == row["read_count"]
                    found = True
                    break

            assert found, f"Interval [{row['start']}, {row['end']}) not found in tree"

    def test_bed_data_types(self):
        """Test BED DataFrame has correct data types."""
        intervals = Intervals(bed_file=self.TEST_BED)

        # Verify data types are correct
        assert intervals.bed["contig"].dtype == object  # string type
        assert intervals.bed["start"].dtype == np.int64
        assert intervals.bed["end"].dtype == np.int64
        assert intervals.bed["read_count"].dtype == np.int64
        assert intervals.bed["fraction"].dtype == np.float64

    def test_tree_intervals_sorted(self):
        """Test that IntervalTree intervals are properly sorted."""
        intervals = Intervals(bed_file=self.TEST_BED_VALID)

        # Convert tree to sorted list
        sorted_intervals = sorted(intervals.tree)

        # Verify intervals are in order
        for i in range(len(sorted_intervals) - 1):
            current_interval = sorted_intervals[i]
            next_interval = sorted_intervals[i + 1]

            # Current interval should start before or at the same position as next
            assert current_interval.begin <= next_interval.begin
