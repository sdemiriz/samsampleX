import pytest
from subsample_reads.Plotter import Plotter


class TestPlotter:
    """Test cases for Plotter class."""

    # Test data file paths - update these if test files change location
    TEST_BAM = "test/data/test-100bp-10count.bam"
    TEST_BED = "test/data/test-100bp-10count.bed"

    def test_file_not_found_bam(self):
        """Test Plotter initialization with nonexistent BAM file."""
        # This should raise FileNotFoundError when trying to load a BAM that doesn't exist
        with pytest.raises(FileNotFoundError):
            Plotter(
                in_bam="tests/DOES_NOT_EXIST.bam",
                map_bam=None,
                out_bam=None,
                bed=self.TEST_BED,
                out_plt="test_output.png",
            )

    def test_file_not_found_bed(self):
        """Test Plotter initialization with nonexistent BED file."""
        # This should raise FileNotFoundError when trying to load a BED that doesn't exist
        with pytest.raises(FileNotFoundError):
            Plotter(
                in_bam=None,
                map_bam=None,
                out_bam=None,
                bed="tests/DOES_NOT_EXIST.bed",
                out_plt="test_output.png",
            )

    def test_init_no_bams_provided(self):
        """Test Plotter initialization with no BAM files (should work but have empty bams dict)."""
        plotter = Plotter(
            in_bam=None,
            map_bam=None,
            out_bam=None,
            bed=self.TEST_BED,
            out_plt="test_output.png",
        )
        # Should initialize but have no BAMs loaded
        assert len(plotter.bams) == 0
        assert len(plotter.colormap) == 0

    def test_init_with_in_bam_only(self):
        """Test Plotter initialization with only input BAM."""
        plotter = Plotter(
            in_bam=self.TEST_BAM,
            map_bam=None,
            out_bam=None,
            bed=self.TEST_BED,
            out_plt="test_output.png",
        )
        assert len(plotter.bams) == 1
        assert "in" in plotter.bams
        assert "in" in plotter.colormap

    def test_init_with_map_bam_only(self):
        """Test Plotter initialization with only mapping BAM."""
        plotter = Plotter(
            in_bam=None,
            map_bam=self.TEST_BAM,
            out_bam=None,
            bed=self.TEST_BED,
            out_plt="test_output.png",
        )
        assert len(plotter.bams) == 1
        assert "map" in plotter.bams

    def test_init_with_out_bam_only(self):
        """Test Plotter initialization with only output BAM."""
        plotter = Plotter(
            in_bam=None,
            map_bam=None,
            out_bam=self.TEST_BAM,
            bed=self.TEST_BED,
            out_plt="test_output.png",
        )
        assert len(plotter.bams) == 1
        assert "out" in plotter.bams

    def test_init_with_in_and_map_bams(self):
        """Test Plotter initialization with input and mapping BAMs."""
        plotter = Plotter(
            in_bam=self.TEST_BAM,
            map_bam=self.TEST_BAM,
            out_bam=None,
            bed=self.TEST_BED,
            out_plt="test_output.png",
        )
        assert len(plotter.bams) == 2
        assert "in" in plotter.bams

    def test_init_with_in_and_out_bams(self):
        """Test Plotter initialization with input and output BAMs."""
        plotter = Plotter(
            in_bam=self.TEST_BAM,
            map_bam=None,
            out_bam=self.TEST_BAM,
            bed=self.TEST_BED,
            out_plt="test_output.png",
        )
        assert len(plotter.bams) == 2
        assert "in" in plotter.bams

    def test_init_with_map_and_out_bams(self):
        """Test Plotter initialization with mapping and output BAMs."""
        plotter = Plotter(
            in_bam=None,
            map_bam=self.TEST_BAM,
            out_bam=self.TEST_BAM,
            bed=self.TEST_BED,
            out_plt="test_output.png",
        )
        assert len(plotter.bams) == 2
        assert "map" in plotter.bams
        assert "out" in plotter.bams

    def test_init_with_all_three_bams(self):
        """Test Plotter initialization with all three BAM files."""
        plotter = Plotter(
            in_bam=self.TEST_BAM,
            map_bam=self.TEST_BAM,
            out_bam=self.TEST_BAM,
            bed=self.TEST_BED,
            out_plt="test_output.png",
        )
        assert len(plotter.bams) == 3
        assert "in" in plotter.bams
        assert "map" in plotter.bams
        assert "out" in plotter.bams

    def test_init_with_no_det_flag(self):
        """Test Plotter initialization with no_det flag."""
        plotter = Plotter(
            in_bam=self.TEST_BAM,
            map_bam=None,
            out_bam=None,
            bed=self.TEST_BED,
            out_plt="test_output.png",
            no_det=True,
        )
        assert plotter.no_det is True

    def test_setup_plot_static_method(self):
        """Test setup_plot static method returns figure and axes."""
        fig, ax_line, ax_bar = Plotter.setup_plot()

        # Verify we got the expected objects back
        assert fig is not None
        assert ax_line is not None
        assert ax_bar is not None

        # Verify it's a figure with 2 subplots
        assert len(fig.axes) == 2
        assert ax_line in fig.axes
        assert ax_bar in fig.axes

        # Clean up to avoid matplotlib warnings
        import matplotlib.pyplot as plt

        plt.close(fig)

    def test_plotter_creates_output_file(self, tmp_path):
        """Test that plot() creates an output PNG file with content."""
        # Use tmp_path fixture to create a temporary output file path
        out_plt = tmp_path / "test_plot.png"

        # Create plotter with one BAM file
        plotter = Plotter(
            in_bam=self.TEST_BAM,
            map_bam=None,
            out_bam=None,
            bed=self.TEST_BED,
            out_plt=str(out_plt),
        )
        plotter.plot()

        # Verify the file was created
        assert out_plt.exists()
        assert out_plt.is_file()
        assert out_plt.stat().st_size > 0
