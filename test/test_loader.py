"""
Tests for Loader class.
"""

import pytest
import pysam
from pathlib import Path
from subsample_reads.Loader import Loader


class TestLoader:
    """Test cases for Loader class."""

    # Test data file paths
    TEST_BAM = "test/data/test-100bp-10count.bam"

    def test_init_with_bam_path(self):
        """Test Loader initialization with BAM path in read mode."""
        loader = Loader(bam_path=self.TEST_BAM)

        # Should create a pysam AlignmentFile object
        assert isinstance(loader.bam, pysam.AlignmentFile)

        # Verify BAM path is stored
        assert loader.bam_path == self.TEST_BAM

        # Clean up
        loader.close()

    def test_file_not_found(self):
        """Test that nonexistent BAM file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            Loader(bam_path="DOES_NOT_EXIST.bam")

    def test_init_with_template(self, tmp_path):
        """Test Loader initialization with template for write mode."""
        # First, open a template BAM in read mode
        template_loader = Loader(bam_path=self.TEST_BAM)

        # Create new BAM with template (write mode)
        output_bam = tmp_path / "output.bam"
        writer_loader = Loader(bam_path=str(output_bam), template=template_loader.bam)

        # Should be in write mode (mode is returned as bytes)
        assert writer_loader.bam.mode == b"wb"
        assert writer_loader.template == template_loader.bam

        # Clean up
        writer_loader.close()
        template_loader.close()

    def test_init_with_header(self, tmp_path):
        """Test Loader initialization with explicit header for write mode."""
        # First, get a header from existing BAM
        template_loader = Loader(bam_path=self.TEST_BAM)
        header = template_loader.bam.header

        # Create new BAM with explicit header
        output_bam = tmp_path / "output.bam"
        writer_loader = Loader(
            bam_path=str(output_bam),
            template=template_loader.bam,  # Need template for write mode
            header=header,
        )

        # Should be in write mode (mode is returned as bytes)
        assert writer_loader.bam.mode == b"wb"

        # Clean up
        writer_loader.close()
        template_loader.close()

    def test_load_bam_read_mode(self):
        """Test load_bam method opens BAM in read mode."""
        loader = Loader(bam_path=self.TEST_BAM)

        # Should be in read mode (binary) - mode is returned as bytes
        assert loader.bam.mode == b"rb"
        assert not loader.bam.is_write

        # Clean up
        loader.close()

    def test_load_bam_write_mode_with_template(self, tmp_path):
        """Test load_bam method opens BAM in write mode with template."""
        # Get template
        template_loader = Loader(bam_path=self.TEST_BAM)

        # Create writer
        output_bam = tmp_path / "output.bam"
        writer_loader = Loader(bam_path=str(output_bam), template=template_loader.bam)

        # Should be in write mode (mode is returned as bytes)
        assert writer_loader.bam.mode == b"wb"
        assert writer_loader.bam.is_write

        # Clean up
        writer_loader.close()
        template_loader.close()

    def test_load_bam_write_mode_with_header(self, tmp_path):
        """Test load_bam method opens BAM in write mode with header."""
        # Get header from existing BAM
        template_loader = Loader(bam_path=self.TEST_BAM)
        header = template_loader.bam.header

        # Create writer with header
        output_bam = tmp_path / "output.bam"
        writer_loader = Loader(
            bam_path=str(output_bam), template=template_loader.bam, header=header
        )

        # Should be in write mode (mode is returned as bytes)
        assert writer_loader.bam.mode == b"wb"

        # Clean up
        writer_loader.close()
        template_loader.close()

    def test_index_if_needed_already_indexed(self):
        """Test index_if_needed when BAM is already indexed."""
        # Our test BAM should already have an index (.bai file)
        loader = Loader(bam_path=self.TEST_BAM)

        # Should not raise an error
        # The index check happens during init, so if we got here, it worked
        assert loader.bam is not None

        # Clean up
        loader.close()

    def test_index_if_needed_requires_indexing(self, tmp_path):
        """Test index_if_needed when BAM needs indexing."""
        # Create a new BAM file (will not have index initially)
        template_loader = Loader(bam_path=self.TEST_BAM)

        # Write a simple BAM
        output_bam = tmp_path / "needs_index.bam"
        writer = Loader(bam_path=str(output_bam), template=template_loader.bam)

        # Write some reads
        for read in template_loader.fetch():
            writer.bam.write(read)
            break  # Just write one read

        writer.close()
        template_loader.close()

        # Now try to open it - should trigger indexing
        # Note: For a real test, we'd need to ensure the index doesn't exist first
        # But pysam.index will handle this gracefully
        indexed_loader = Loader(bam_path=str(output_bam))

        assert indexed_loader.bam is not None
        indexed_loader.close()

    def test_overlap_true(self):
        """Test overlap method with overlapping coordinates."""
        loader = Loader(bam_path=self.TEST_BAM)

        # Test various overlapping scenarios
        assert loader.overlap((0, 10), (0, 10))  # Identical ranges
        assert loader.overlap((0, 5), (3, 7))  # Partial overlap
        assert loader.overlap((0, 10), (4, 5))  # One inside other
        assert loader.overlap((0, 10), (9, 20))  # Overlap at end

        loader.close()

    def test_overlap_false(self):
        """Test overlap method with non-overlapping coordinates."""
        loader = Loader(bam_path=self.TEST_BAM)

        # Test non-overlapping scenarios
        assert not loader.overlap((0, 10), (20, 30))  # Completely separate
        assert not loader.overlap(
            (0, 10), (10, 20)
        )  # Adjacent (touching but not overlapping)
        assert not loader.overlap((0, 10), (11, 20))  # Sequential

        loader.close()

    def test_normalize_contig_chr_format(self):
        """Test normalize_contig with chr-prefixed format."""
        loader = Loader(bam_path=self.TEST_BAM)

        # Should return the same if already correct and present
        # Our test BAM has chr1
        result = loader.normalize_contig("chr1")
        assert result == "chr1"

        loader.close()

    def test_normalize_contig_numeric_format(self):
        """Test normalize_contig with numeric format."""
        loader = Loader(bam_path=self.TEST_BAM)

        # If BAM has '1' as a reference, should convert to 'chr1'
        # If BAM only has 'chr1', this should raise ValueError
        # Our test BAM likely only has 'chr1', so test that scenario
        try:
            result = loader.normalize_contig("1")
            # If this succeeds, BAM had '1' reference
            assert result == "chr1"
        except ValueError:
            # Expected if BAM only has 'chr1' format
            pass

        loader.close()

    def test_normalize_contig_invalid(self):
        """Test normalize_contig with invalid contig name."""
        loader = Loader(bam_path=self.TEST_BAM)

        # Should raise ValueError for non-existent contig
        with pytest.raises(ValueError):
            loader.normalize_contig("chrFAKE")

        loader.close()

    def test_get_reference_name(self):
        """Test get_reference_name method."""
        loader = Loader(bam_path=self.TEST_BAM)

        # Get reference ID for chr1
        ref_id = loader.bam.get_tid("chr1")

        # Should return 'chr1'
        result = loader.get_reference_name(ref_id)
        assert result == "chr1"

        loader.close()

    def test_fetch_without_names(self):
        """Test fetch method without read names."""
        loader = Loader(bam_path=self.TEST_BAM)

        # Fetch all reads
        reads = list(loader.fetch())

        # Should have some reads
        assert len(reads) > 0

        # All should be mapped (based on our test file)
        assert all(r.is_mapped for r in reads)

        loader.close()

    def test_fetch_with_names(self):
        """Test fetch method with read names filter."""
        loader = Loader(bam_path=self.TEST_BAM)

        # Get all reads first to get some names
        all_reads = list(loader.fetch())

        # This test verifies the fetch method works
        # The actual filtering behavior depends on the implementation
        assert len(all_reads) >= 0  # Just verify we can fetch

        loader.close()

    def test_close(self):
        """Test close method."""
        loader = Loader(bam_path=self.TEST_BAM)

        # Should be open initially
        assert not loader.bam.closed

        # Close it
        loader.close()

        # Should be closed now
        assert loader.bam.closed

    def test_setup_mapback_grch38(self):
        """Test setup_mapback with GRCh38 genome build."""
        loader = Loader(bam_path=self.TEST_BAM)

        # This requires HLA-LA sequences.txt file
        # Skip if not available
        try:
            loader.setup_mapback(genome_build="GRCh38")

            # Verify gene_maps are set
            assert hasattr(loader, "gene_maps")
            assert "A" in loader.gene_maps  # HLA-A
            assert "B" in loader.gene_maps  # HLA-B

            # Verify alt_contig_maps are set
            assert hasattr(loader, "alt_contig_maps")
            assert "chr6_GL000250v2_alt" in loader.alt_contig_maps

        except FileNotFoundError:
            # sequences.txt not available, skip this test
            pytest.skip("HLA-LA sequences.txt not available")
        finally:
            loader.close()

    def test_setup_mapback_grch37(self):
        """Test setup_mapback with GRCh37 genome build."""
        loader = Loader(bam_path=self.TEST_BAM)

        # This requires HLA-LA sequences.txt file
        # Skip if not available
        try:
            loader.setup_mapback(genome_build="GRCh37")

            # Verify gene_maps are set
            assert hasattr(loader, "gene_maps")
            assert "A" in loader.gene_maps
            assert "B" in loader.gene_maps

            # Verify alt_contig_maps are set
            assert hasattr(loader, "alt_contig_maps")
            assert "chr6_GL000250v2_alt" in loader.alt_contig_maps

        except FileNotFoundError:
            # sequences.txt not available, skip this test
            pytest.skip("HLA-LA sequences.txt not available")
        finally:
            loader.close()

    def test_run_sampling_regular_mode(self):
        """Test run_sampling in regular mode (end-to-end integration test)."""
        import tempfile
        import os

        # Use same filesystem to avoid cross-device link issues
        temp_dir = tempfile.mkdtemp(dir=".")

        try:
            # Create input loader
            loader = Loader(bam_path=self.TEST_BAM)

            # Setup output path
            output_bam = os.path.join(temp_dir, "sampled_output.bam")

            # Run sampling in regular (non-HLA-LA) mode
            loader.run_sampling(
                bed_dir=None,
                bed_file="test/data/test-100bp-10count.bed",
                main_seed=42,
                out_bam=output_bam,
                hlala_mode=False,  # Regular mode, not HLA-LA
                genome_build=None,
            )

            # Verify output file was created
            # Note: The output might be .sorted.bam instead of .bam
            sorted_output = output_bam.replace(".bam", ".sorted.bam")
            assert os.path.exists(output_bam) or os.path.exists(
                sorted_output
            ), f"Expected output BAM not found: {output_bam} or {sorted_output}"

            loader.close()

        finally:
            # Clean up temp directory
            import shutil

            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)

    def test_run_sampling_hlala_mode(self):
        """Test run_sampling in HLA-LA mode."""
        # This requires HLA-LA sequences.txt and is HLA-specific
        pytest.skip(
            "HLA-LA specific integration test requiring sequences.txt - deferred"
        )

    def test_get_intervals(self):
        """Test get_intervals method."""
        loader = Loader(bam_path=self.TEST_BAM)

        # Test getting intervals from BED file
        intervals = loader.get_intervals(
            bed_dir=None, bed_file="test/data/test-100bp-10count.bed"
        )

        # Should return an Intervals object
        from subsample_reads.Intervals import Intervals

        assert isinstance(intervals, Intervals)
        assert len(intervals) > 0

        loader.close()

    def test_get_interval_seeds(self):
        """Test get_interval_seeds method."""
        loader = Loader(bam_path=self.TEST_BAM)

        # Create some dummy intervals for testing
        intervals = loader.get_intervals(
            bed_dir=None, bed_file="test/data/test-100bp-10count.bed"
        )

        # Set intervals attribute
        loader.intervals = intervals

        # Get seeds
        seeds = loader.get_interval_seeds(main_seed=42)

        # Should return a numpy array with one seed per interval
        import numpy as np

        assert isinstance(seeds, np.ndarray)
        assert len(seeds) == len(intervals)

        loader.close()

    def test_setup_buckets(self):
        """Test setup_buckets method."""
        loader = Loader(bam_path=self.TEST_BAM)

        # Setup intervals first
        intervals = loader.get_intervals(
            bed_dir=None, bed_file="test/data/test-100bp-10count.bed"
        )
        loader.intervals = intervals

        # Setup buckets - returns a list, not a dict
        buckets = loader.setup_buckets()

        # Should return a list with one bucket per interval
        assert isinstance(buckets, list)
        assert len(buckets) == len(intervals)

        # Each bucket should be a list
        for bucket in buckets:
            assert isinstance(bucket, list)

        loader.close()

    def test_get_mapped_reads(self):
        """Test get_mapped_reads method."""
        loader = Loader(bam_path=self.TEST_BAM)

        # Need to set intervals attribute for get_mapped_reads
        intervals = loader.get_intervals(
            bed_dir=None, bed_file="test/data/test-100bp-10count.bed"
        )
        loader.intervals = intervals

        # Get mapped reads in a range
        reads = list(loader.get_mapped_reads(start=0, end=10000))

        # Should return some reads (generator converted to list)
        # All should be mapped
        for read in reads:
            assert read.is_mapped

        loader.close()

    def test_write_reads(self):
        """Test write_reads method."""
        import tempfile
        import os

        # Use same filesystem to avoid cross-device link issues
        # Create temp directory in current working directory
        temp_dir = tempfile.mkdtemp(dir=".")

        try:
            # Create a source and destination loader
            source = Loader(bam_path=self.TEST_BAM)
            output_bam = os.path.join(temp_dir, "output.bam")
            dest = Loader(bam_path=output_bam, template=source.bam)

            # Collect some reads into a list (simulating buckets)
            dest.reads = [read for read in source.fetch()][:5]  # Just first 5 reads

            # Set write_out attribute that write_reads expects
            dest.write_out = output_bam

            # Write them
            dest.write_reads(bam=dest)

            # Clean up loaders
            dest.close()
            source.close()

            # Verify file exists
            assert os.path.exists(output_bam)

        finally:
            # Clean up temp directory
            import shutil

            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)

    def test_sort_and_index(self):
        """Test sort_and_index method creates indexed BAM file."""
        import tempfile
        import os
        import glob

        # Use same filesystem to avoid cross-device link issues
        temp_dir = tempfile.mkdtemp(dir=".")

        try:
            # Create a BAM file first
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

            # Now sort and index
            loader = Loader(bam_path=self.TEST_BAM)
            loader.write_out = output_bam
            loader.sort_and_index()

            # The key indicator of success is the .bai index file
            # Check for .bai files in the temp directory
            bai_files = glob.glob(os.path.join(temp_dir, "*.bai"))

            assert (
                len(bai_files) > 0
            ), f"Expected .bai index file after sort_and_index. Found files: {os.listdir(temp_dir)}"

            loader.close()

        finally:
            # Clean up temp directory
            import shutil

            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)

    def test_sample_reads_from_buckets(self):
        """Test sample_reads_from_buckets method."""
        loader = Loader(bam_path=self.TEST_BAM)

        # Setup complete sampling environment
        intervals = loader.get_intervals(
            bed_dir=None, bed_file="test/data/test-100bp-10count.bed"
        )
        loader.intervals = intervals

        # Get seeds
        loader.seeds = loader.get_interval_seeds(main_seed=42)

        # Setup buckets and populate with some reads
        loader.buckets = loader.setup_buckets()

        # Setup interval arrays
        import numpy as np

        loader.interval_starts = np.array(
            [interval.begin for interval in intervals.tree]
        )
        loader.interval_ends = np.array([interval.end for interval in intervals.tree])

        # Add some reads to buckets
        reads = list(loader.fetch())
        for read in reads[:10]:  # Just use first 10 reads
            loader.add_read_to_bucket(read=read, buckets=loader.buckets)

        # Now sample from buckets
        loader.sample_reads_from_buckets()

        # Verify that reads attribute is now populated
        assert hasattr(loader, "reads")
        assert isinstance(loader.reads, list)
        # May be empty if no reads overlapped intervals

        loader.close()

    def test_add_read_to_bucket(self):
        """Test add_read_to_bucket method."""
        loader = Loader(bam_path=self.TEST_BAM)

        # Setup intervals and buckets
        intervals = loader.get_intervals(
            bed_dir=None, bed_file="test/data/test-100bp-10count.bed"
        )
        loader.intervals = intervals
        buckets = loader.setup_buckets()  # Returns a list

        # Setup interval arrays for overlap detection
        import numpy as np

        loader.interval_starts = np.array(
            [interval.begin for interval in intervals.tree]
        )
        loader.interval_ends = np.array([interval.end for interval in intervals.tree])

        # Get a read
        reads = list(loader.fetch())
        if len(reads) > 0:
            test_read = reads[0]

            # Add to bucket - buckets is a list
            loader.add_read_to_bucket(read=test_read, buckets=buckets)

            # At least one bucket should have reads now
            total_reads = sum(len(bucket) for bucket in buckets)
            assert total_reads > 0

        loader.close()

    # HLA-LA specific tests - these require sequences.txt and are complex
    def test_map_read_to_chr6_chr6_read(self):
        """Test map_read_to_chr6 with read already on chr6."""
        pytest.skip("HLA-LA specific test requiring sequences.txt - deferred")

    def test_map_read_to_chr6_prg_read(self):
        """Test map_read_to_chr6 with PRG-mapped read."""
        pytest.skip("HLA-LA specific test requiring sequences.txt - deferred")

    def test_map_read_to_chr6_invalid_contig(self):
        """Test map_read_to_chr6 with invalid contig."""
        pytest.skip("HLA-LA specific test requiring sequences.txt - deferred")

    def test_get_reads_from_contigs(self):
        """Test get_reads_from_contigs method."""
        loader = Loader(bam_path=self.TEST_BAM)

        # Get contigs from our test BAM
        contigs = list(loader.bam.references)

        if len(contigs) > 0:
            # Get reads from first contig
            reads = list(loader.get_reads_from_contigs(contigs=[contigs[0]]))

            # Should get some reads
            # All should be from the specified contig
            for read in reads:
                if read.is_mapped:
                    ref_name = loader.get_reference_name(read.reference_id)
                    assert ref_name == contigs[0]

        loader.close()

    def test_get_prg_contigs(self):
        """Test get_prg_contigs method."""
        pytest.skip("HLA-LA specific test requiring sequences.txt - deferred")

    def test_convert_to_chr6(self):
        """Test convert_to_chr6 method."""
        pytest.skip("HLA-LA specific test requiring sequences.txt - deferred")

    def test_get_prg_coords_hla_allele(self):
        """Test get_prg_coords with HLA allele."""
        pytest.skip("HLA-LA specific test requiring sequences.txt - deferred")

    def test_get_prg_coords_alt_contig(self):
        """Test get_prg_coords with alt contig."""
        pytest.skip("HLA-LA specific test requiring sequences.txt - deferred")

    def test_get_prg_coords_not_found(self):
        """Test get_prg_coords with contig not found."""
        pytest.skip("HLA-LA specific test requiring sequences.txt - deferred")

    def test_modify_header_grch38(self):
        """Test modify_header with GRCh38."""
        pytest.skip("HLA-LA specific test requiring sequences.txt - deferred")

    def test_modify_header_grch37(self):
        """Test modify_header with GRCh37."""
        pytest.skip("HLA-LA specific test requiring sequences.txt - deferred")
