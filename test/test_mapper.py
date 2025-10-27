"""
Tests for Mapper class.
"""

import pytest
import pandas as pd
from pathlib import Path
from samsampleX.Mapper import Mapper


class TestMapper:
    """Test cases for Mapper class."""

    # Test data file paths
    TEST_BAM = "test/data/test-100bp-10count.bam"

    def test_file_not_found(self, tmp_path):
        """Test Mapper with nonexistent BAM file."""
        with pytest.raises(FileNotFoundError):
            Mapper(
                bam_paths=["DOES_NOT_EXIST.bam"],
                contig="chr1",
                start="100",
                end="1000",
                interval_count="10",
                interval_length=None,
                bed_dir=str(tmp_path / "bed"),
            )

    def test_no_interval_parameters(self, tmp_path):
        """Test that providing neither length nor count raises ValueError."""
        with pytest.raises(ValueError) as exc_info:
            Mapper(
                bam_paths=[self.TEST_BAM],
                contig="chr1",
                start="100",
                end="1000",
                interval_count=None,
                interval_length=None,
                bed_dir=str(tmp_path / "bed"),
            )
        assert "No interval length or count provided" in str(exc_info.value)

    def test_init_with_interval_count(self, tmp_path):
        """Test Mapper initialization with interval count."""
        bed_dir = tmp_path / "bed"

        mapper = Mapper(
            bam_paths=[self.TEST_BAM],
            contig="chr1",
            start="100",
            end="1100",
            interval_count="10",
            interval_length=None,
            bed_dir=str(bed_dir),
        )

        # Verify interval settings
        assert mapper.interval_count == 10
        assert mapper.interval_length is None

        # Verify BED file was created
        expected_bed = bed_dir / "test-100bp-10count.bed"
        assert expected_bed.exists()

    def test_init_with_interval_length(self, tmp_path):
        """Test Mapper initialization with interval length."""
        bed_dir = tmp_path / "bed"

        mapper = Mapper(
            bam_paths=[self.TEST_BAM],
            contig="chr1",
            start="100",
            end="1100",
            interval_count=None,
            interval_length="100",
            bed_dir=str(bed_dir),
        )

        # Verify interval settings
        assert mapper.interval_length == 100
        assert mapper.interval_count is None

        # Verify BED file was created
        expected_bed = bed_dir / "test-100bp-10count.bed"
        assert expected_bed.exists()

    def test_bed_output_columns(self, tmp_path):
        """Test that output BED file has correct columns."""
        bed_dir = tmp_path / "bed"

        Mapper(
            bam_paths=[self.TEST_BAM],
            contig="chr1",
            start="100",
            end="1100",
            interval_count="10",
            interval_length=None,
            bed_dir=str(bed_dir),
        )

        # Read the output BED file
        bed_file = bed_dir / "test-100bp-10count.bed"
        df = pd.read_csv(bed_file, sep="\t", header=None)

        # Should have 4 columns: contig, start, end, read_count
        assert df.shape[1] == 4

        # Verify column data types
        assert df[0].dtype == object  # contig (string)
        assert df[1].dtype in [int, "int64"]  # start
        assert df[2].dtype in [int, "int64"]  # end
        assert df[3].dtype in [int, "int64"]  # read_count

    def test_init_with_bed_files(self, tmp_path):
        """Test Mapper initialization with specific BED file paths."""
        output_bed = tmp_path / "custom_output.bed"

        mapper = Mapper(
            bam_paths=[self.TEST_BAM],
            contig="chr1",
            start="100",
            end="1100",
            interval_count="10",
            interval_length=None,
            bed=[str(output_bed)],
        )

        # Verify the custom BED path was used
        assert output_bed.exists()
        assert len(mapper.bed_paths) == 1
        assert mapper.bed_paths[0] == output_bed

    def test_init_with_bed_dir(self, tmp_path):
        """Test Mapper initialization with BED directory."""
        bed_dir = tmp_path / "custom_bed_dir"

        mapper = Mapper(
            bam_paths=[self.TEST_BAM],
            contig="chr1",
            start="100",
            end="1100",
            interval_count="5",
            interval_length=None,
            bed_dir=str(bed_dir),
        )

        # Verify directory was created
        assert bed_dir.exists()
        assert bed_dir.is_dir()

        # Verify BED file was created in the directory
        expected_bed = bed_dir / "test-100bp-10count.bed"
        assert expected_bed.exists()

    def test_make_bed_dir(self, tmp_path):
        """Test make_bed_dir creates directory."""
        test_dir = tmp_path / "test_bed_dir"

        mapper = Mapper(
            bam_paths=[self.TEST_BAM],
            contig="chr1",
            start="100",
            end="1100",
            interval_count="5",
            interval_length=None,
            bed_dir=str(test_dir),
        )

        assert mapper.bed_dir == test_dir
        assert test_dir.exists()

    def test_get_bed_paths(self, tmp_path):
        """Test get_bed_paths generates correct paths."""
        bed_dir = tmp_path / "bed"

        mapper = Mapper(
            bam_paths=[self.TEST_BAM],
            contig="chr1",
            start="100",
            end="1100",
            interval_count="5",
            interval_length=None,
            bed_dir=str(bed_dir),
        )

        # Should have one BED path matching the BAM name
        assert len(mapper.bed_paths) == 1
        assert mapper.bed_paths[0].name == "test-100bp-10count.bed"
        assert mapper.bed_paths[0].parent == bed_dir

    def test_setup_intervals_with_length(self, tmp_path):
        """Test setup_intervals with interval length."""
        mapper = Mapper(
            bam_paths=[self.TEST_BAM],
            contig="chr1",
            start="100",
            end="1100",
            interval_count=None,
            interval_length="200",
            bed_dir=str(tmp_path / "bed"),
        )

        assert mapper.interval_length == 200
        assert mapper.interval_count is None

    def test_setup_intervals_with_count(self, tmp_path):
        """Test setup_intervals with interval count."""
        mapper = Mapper(
            bam_paths=[self.TEST_BAM],
            contig="chr1",
            start="100",
            end="1100",
            interval_count="15",
            interval_length=None,
            bed_dir=str(tmp_path / "bed"),
        )

        assert mapper.interval_count == 15
        assert mapper.interval_length is None

    def test_setup_intervals_no_parameters(self, tmp_path):
        """Test that no interval parameters raises ValueError."""
        with pytest.raises(ValueError):
            Mapper(
                bam_paths=[self.TEST_BAM],
                contig="chr1",
                start="100",
                end="1100",
                interval_count=None,
                interval_length=None,
                bed_dir=str(tmp_path / "bed"),
            )

    def test_load_bams(self, tmp_path):
        """Test load_bams initializes Loaders."""
        mapper = Mapper(
            bam_paths=[self.TEST_BAM],
            contig="chr1",
            start="100",
            end="1100",
            interval_count="5",
            interval_length=None,
            bed_dir=str(tmp_path / "bed"),
        )

        assert len(mapper.bams) == 1
        assert hasattr(mapper.bams[0], "bam")  # Should be a Loader instance

    def test_make_beds_with_length(self, tmp_path):
        """Test make_beds creates correct DataFrames with length."""
        mapper = Mapper(
            bam_paths=[self.TEST_BAM],
            contig="chr1",
            start="100",
            end="1000",
            interval_count=None,
            interval_length="100",
            bed_dir=str(tmp_path / "bed"),
        )

        # Should have one BED DataFrame per BAM
        assert len(mapper.beds) == 1

        # Check DataFrame structure
        bed = mapper.beds[0]
        assert "contig" in bed.columns
        assert "start" in bed.columns
        assert "end" in bed.columns
        assert "read_count" in bed.columns

    def test_make_beds_with_count(self, tmp_path):
        """Test make_beds creates correct DataFrames with count."""
        mapper = Mapper(
            bam_paths=[self.TEST_BAM],
            contig="chr1",
            start="100",
            end="1100",
            interval_count="10",
            interval_length=None,
            bed_dir=str(tmp_path / "bed"),
        )

        # Should have 10 intervals
        bed = mapper.beds[0]
        assert len(bed) == 10

    def test_get_interval_boundaries_with_length(self, tmp_path):
        """Test get_interval_boundaries with interval length."""
        mapper = Mapper(
            bam_paths=[self.TEST_BAM],
            contig="chr1",
            start="0",
            end="1000",
            interval_count=None,
            interval_length="100",
            bed_dir=str(tmp_path / "bed"),
        )

        # With length=100 and range 0-1000, should have 10 intervals
        # boundaries: [0, 100, 200, ..., 900, 1000]
        bed = mapper.beds[0]
        assert len(bed) == 10

    def test_get_interval_boundaries_with_count(self, tmp_path):
        """Test get_interval_boundaries with interval count."""
        mapper = Mapper(
            bam_paths=[self.TEST_BAM],
            contig="chr1",
            start="0",
            end="1000",
            interval_count="5",
            interval_length=None,
            bed_dir=str(tmp_path / "bed"),
        )

        # Should have exactly 5 intervals
        bed = mapper.beds[0]
        assert len(bed) == 5

    def test_populate_read_counts(self, tmp_path):
        """Test that read counts are populated in BED DataFrames."""
        mapper = Mapper(
            bam_paths=[self.TEST_BAM],
            contig="chr1",
            start="100",
            end="1100",
            interval_count="5",
            interval_length=None,
            bed_dir=str(tmp_path / "bed"),
        )

        # Read counts should not be -1 (default) after population
        bed = mapper.beds[0]
        assert all(bed["read_count"] >= 0)

    def test_write_beds(self, tmp_path):
        """Test that BED files are written correctly."""
        bed_dir = tmp_path / "bed"
        output_bed = bed_dir / "test-100bp-10count.bed"

        Mapper(
            bam_paths=[self.TEST_BAM],
            contig="chr1",
            start="100",
            end="1100",
            interval_count="5",
            interval_length=None,
            bed_dir=str(bed_dir),
        )

        # File should exist
        assert output_bed.exists()

        # Should be readable as DataFrame
        df = pd.read_csv(output_bed, sep="\t", header=None)
        assert len(df) == 5  # 5 intervals
