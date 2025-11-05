"""
Tests for Sampler class.
"""

import pytest
import pysam
import numpy as np
import tempfile
import os
import glob
from samsampleX.Sampler import Sampler
from samsampleX.Intervals import Intervals


class TestSampler:
    """Test cases for Sampler class."""

    # Test data file paths
    TEST_BAM = "test/data/test-100bp-10count.bam"
    TEST_BED = "test/data/test-100bp-10count.bed"

    def test_init(self):
        """Test Sampler initialization."""
        sampler = Sampler(bam_path=self.TEST_BAM, bed_file=self.TEST_BED)

        # Should have a Loader instance
        assert hasattr(sampler, "bam")
        assert sampler.bam is not None

        # Should have intervals
        assert hasattr(sampler, "intervals")
        assert isinstance(sampler.intervals, Intervals)

        # Clean up
        sampler.bam.close()

    def test_overlap_true(self):
        """Test overlap method with overlapping coordinates."""
        sampler = Sampler(bam_path=self.TEST_BAM, bed_file=self.TEST_BED)

        # Test various overlapping scenarios
        assert sampler.overlap((0, 10), (0, 10))  # Identical ranges
        assert sampler.overlap((0, 5), (3, 7))  # Partial overlap
        assert sampler.overlap((0, 10), (4, 5))  # One inside other
        assert sampler.overlap((0, 10), (9, 20))  # Overlap at end

        sampler.bam.close()

    def test_overlap_false(self):
        """Test overlap method with non-overlapping coordinates."""
        sampler = Sampler(bam_path=self.TEST_BAM, bed_file=self.TEST_BED)

        # Test non-overlapping scenarios
        assert not sampler.overlap((0, 10), (20, 30))  # Completely separate
        assert not sampler.overlap(
            (0, 10), (10, 20)
        )  # Adjacent (touching but not overlapping)
        assert not sampler.overlap((0, 10), (11, 20))  # Sequential

        sampler.bam.close()

    def test_get_intervals(self):
        """Test get_intervals method."""
        sampler = Sampler(bam_path=self.TEST_BAM, bed_file=self.TEST_BED)

        # Should return an Intervals object
        assert isinstance(sampler.intervals, Intervals)
        assert len(sampler.intervals) > 0

        sampler.bam.close()

    def test_get_interval_seeds(self):
        """Test get_interval_seeds method."""
        sampler = Sampler(bam_path=self.TEST_BAM, bed_file=self.TEST_BED)

        # Get seeds
        seeds = sampler.get_interval_seeds(main_seed=42)

        # Should return a numpy array with one seed per interval
        assert isinstance(seeds, np.ndarray)
        assert len(seeds) == len(sampler.intervals)

        sampler.bam.close()

    def test_setup_buckets(self):
        """Test setup_buckets method."""
        sampler = Sampler(bam_path=self.TEST_BAM, bed_file=self.TEST_BED)

        # Setup buckets - returns a list, not a dict
        buckets = sampler.setup_buckets()

        # Should return a list with one bucket per interval
        assert isinstance(buckets, list)
        assert len(buckets) == len(sampler.intervals)

        # Each bucket should be a list
        for bucket in buckets:
            assert isinstance(bucket, list)

        sampler.bam.close()

    def test_get_mapped_reads(self):
        """Test get_mapped_reads method."""
        sampler = Sampler(bam_path=self.TEST_BAM, bed_file=self.TEST_BED)

        # Get mapped reads in a range
        reads = list(sampler.get_mapped_reads(start=0, end=10000))

        # Should return some reads (generator converted to list)
        # All should be mapped
        for read in reads:
            assert read.is_mapped

        sampler.bam.close()

    def test_normalize_contig(self):
        """Test normalize_contig method."""
        sampler = Sampler(bam_path=self.TEST_BAM, bed_file=self.TEST_BED)

        # Should return the same if already correct and present
        # Our test BAM has chr1
        result = sampler.normalize_contig("chr1")
        assert result == "chr1"

        sampler.bam.close()

    def test_sort_and_index(self):
        """Test sort_and_index method creates indexed BAM file."""
        # Use same filesystem to avoid cross-device link issues
        temp_dir = tempfile.mkdtemp(dir=".")

        try:
            # Create a BAM file first
            from samsampleX.Loader import Loader

            source = Loader(bam_path=self.TEST_BAM)
            output_bam = os.path.join(temp_dir, "output.bam")
            dest = Loader(bam_path=output_bam, template=source.bam)

            # Write some reads
            for i, read in enumerate(source.fetch()):
                dest.bam.write(read)
                if i >= 5:
                    break

            dest.close()
            source.close()

            # Now sort and index using Sampler
            sampler = Sampler(bam_path=self.TEST_BAM, bed_file=self.TEST_BED)
            sampler.write_out = output_bam
            sampler.sort_and_index()

            # The key indicator of success is the .bai index file
            # Check for .bai files in the temp directory
            bai_files = glob.glob(os.path.join(temp_dir, "*.bai"))

            assert (
                len(bai_files) > 0
            ), f"Expected .bai index file after sort_and_index. Found files: {os.listdir(temp_dir)}"

            sampler.bam.close()

        finally:
            # Clean up temp directory
            import shutil

            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)

    def test_run_sampling_regular_mode(self):
        """Test run_sampling in regular mode (end-to-end integration test)."""
        # Use same filesystem to avoid cross-device link issues
        temp_dir = tempfile.mkdtemp(dir=".")

        try:
            # Create sampler
            sampler = Sampler(bam_path=self.TEST_BAM, bed_file=self.TEST_BED)

            # Setup output path
            output_bam = os.path.join(temp_dir, "sampled_output.bam")

            # Run sampling in regular (non-HLA-LA) mode
            sampler.run_sampling(
                main_seed=42,
                out_bam=output_bam,
            )

            # Verify output file was created
            # Note: The output might be .sorted.bam instead of .bam
            sorted_output = output_bam.replace(".bam", ".sorted.bam")
            assert os.path.exists(output_bam) or os.path.exists(
                sorted_output
            ), f"Expected output BAM not found: {output_bam} or {sorted_output}"

        finally:
            # Clean up temp directory
            import shutil

            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
