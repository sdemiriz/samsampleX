"""
Tests for Comparator class.
"""

import pytest
import pandas as pd
from subsample_reads.Comparator import Comparator


class TestComparator:

    TEST_BAM = "test/data/test-100bp-10count.bam"

    def test_file_not_found_left_bam(self, tmp_path):
        """Test Comparator with nonexistent left BAM file."""
        out_csv = tmp_path / "output.csv"
        with pytest.raises(FileNotFoundError):
            Comparator(
                bam_left_path="DOES_NOT_EXIST.bam",
                bam_right_path=self.TEST_BAM,
                out=str(out_csv),
            )

    def test_file_not_found_right_bam(self, tmp_path):
        """Test Comparator with nonexistent right BAM file."""
        out_csv = tmp_path / "output.csv"
        with pytest.raises(FileNotFoundError):
            Comparator(
                bam_left_path=self.TEST_BAM,
                bam_right_path="DOES_NOT_EXIST.bam",
                out=str(out_csv),
            )

    def test_comparator_creates_output_csv(self, tmp_path):
        """Test that Comparator creates output CSV file."""
        out_csv = tmp_path / "comparison.csv"

        Comparator(
            bam_left_path=self.TEST_BAM,
            bam_right_path=self.TEST_BAM,
            out=str(out_csv),
        )

        assert out_csv.exists()
        assert out_csv.is_file()
        assert out_csv.stat().st_size > 0

    def test_comparator_output_columns(self, tmp_path):
        """Test that output CSV has expected columns."""
        out_csv = tmp_path / "comparison.csv"

        Comparator(
            bam_left_path=self.TEST_BAM,
            bam_right_path=self.TEST_BAM,
            out=str(out_csv),
        )

        df = pd.read_csv(out_csv, sep="\t")
        expected_cols = [
            "name",
            "ref_name_l",
            "ref_name_r",
            "ref_pos_l",
            "ref_pos_r",
            "next_ref_name_l",
            "next_ref_name_r",
            "next_ref_pos_l",
            "next_ref_pos_r",
            "length_l",
            "length_r",
        ]

        for col in expected_cols:
            assert col in df.columns, f"Expected column '{col}' not found in output"

        assert len(df) > 0, "Output CSV should contain data rows"

    def test_comparator_output_has_data(self, tmp_path):
        """Test that Comparator produces rows when comparing identical BAM files."""
        out_csv = tmp_path / "comparison.csv"

        Comparator(
            bam_left_path=self.TEST_BAM,
            bam_right_path=self.TEST_BAM,
            out=str(out_csv),
        )

        df = pd.read_csv(out_csv, sep="\t")
        assert len(df) > 0
        assert (df["ref_name_l"] == df["ref_name_r"]).all()
        assert (df["ref_pos_l"] == df["ref_pos_r"]).all()
