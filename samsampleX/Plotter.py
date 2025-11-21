import logging
from typing import Optional

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from samsampleX.Loader import Loader
from samsampleX.Intervals import Intervals
from samsampleX.FileHandler import FileHandler

logger = logging.getLogger(__name__)


class Plotter(FileHandler):

    def __init__(
        self,
        in_bam: Optional[str],
        map_bam: Optional[str],
        out_bam: Optional[str],
        bed: str,
        out_plt: str,
    ) -> None:
        """
        Initialize Plotter with argument names.
        """
        logger.info("[PLOTTER] - Initialize")

        self.intervals = Intervals(bed_path=bed)

        self.bams = {}
        self.colormap = {}
        if in_bam:
            self.in_bam = Loader(bam_path=in_bam)
            self.bams["in"] = self.in_bam
            self.colormap["in"] = "#66c2a5"
        if map_bam:
            self.map_bam = Loader(bam_path=map_bam)
            self.bams["map"] = self.map_bam
            self.colormap["map"] = "#fc8d62"
        if out_bam:
            self.out_bam = Loader(bam_path=out_bam)
            self.bams["out"] = self.out_bam
            self.colormap["out"] = "#8da0cb"

        self.out_plt = out_plt

        logger.info("[PLOTTER] - Complete initialization")

    def get_pileups(self) -> list:
        """Pileup BAMs for the defined region"""
        logger.info("[PLOTTER] - Pileup BAMs")

        pileups = [
            bam.bam.pileup(
                contig=bam.normalize_contig(contig=self.intervals.contig),
                start=self.intervals.start,
                end=self.intervals.end,
            )
            for bam in self.bams.values()
        ]

        logger.info("[PLOTTER] - Complete pileup BAMs")
        return pileups

    def get_counts(self, bam: Loader, start: int, end: int) -> list:
        """Pileup BAMs for the defined region"""
        logger.info("[PLOTTER] - Pileup BAMs")

        counts = [
            bam.bam.count(
                contig=bam.normalize_contig(contig=self.intervals.contig),
                start=start,
                end=end,
            )
        ]

        logger.info("[PLOTTER] - Complete pileup BAMs")
        return counts

    def plot(self) -> None:
        """Plot provided BAM file pileups and read counts in separate subplots"""
        logger.info("[PLOTTER] - Begin plotting")
        fig, ax_depth, ax_count = self.setup_plot()

        # Add content to bar and line plot
        self.add_depth_plot(ax=ax_depth)
        self.add_read_count_plot(ax=ax_count)

        # Calculate and report error metrics
        depth_error, count_error = self.calculate_and_report_error()

        # Add annotations with error values in titles
        self.add_annotations(
            ax_depth=ax_depth,
            ax_count=ax_count,
            depth_error=depth_error,
            count_error=count_error,
        )

        logger.info("[PLOTTER] - Save plot")
        plt.savefig(self.out_plt, dpi=600, bbox_inches="tight")

    @staticmethod
    def setup_plot() -> tuple:
        """Set up the figure with two separate subplots side by side"""
        logger.info("[PLOTTER] - Setup plot")

        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

        fig, (ax_depth, ax_count) = plt.subplots(
            1, 2, figsize=(12, 5), layout="constrained"
        )

        def custom_formatter(x, pos):
            if x >= 1e6:
                return f"{x/1e6:.3f}m"  # Millions
            else:
                return f"{x:.0f}"

        formatter = mticker.FuncFormatter(custom_formatter)

        ax_depth.xaxis.set_major_formatter(formatter)
        ax_count.xaxis.set_major_formatter(formatter)

        return fig, ax_depth, ax_count

    def add_depth_plot(self, ax) -> None:
        """Add line plot representing coverage depth to supplied axis"""
        logger.info("[PLOTTER] - Plot pileups")

        for p, l, c in zip(
            self.get_pileups(), self.bams.keys(), self.colormap.values()
        ):
            pileup = pd.DataFrame(
                [(a.reference_pos, a.nsegments) for a in p], columns=["coord", "depth"]
            )
            ax.plot(pileup["coord"], pileup["depth"], label=l, color=c, linewidth=1)
            ax.fill_between(pileup["coord"], pileup["depth"], alpha=0.3, color=c)

    def add_annotations(
        self,
        ax_depth,
        ax_count,
        depth_error: Optional[float] = None,
        count_error: Optional[float] = None,
    ) -> None:
        """Add various text labels to supplied axes"""
        logger.info("[PLOTTER] - Add axes, title, legend")

        # Annotations for depth plot (left subplot)
        self.remove_spines(ax=ax_depth)
        self.add_horizontal_lines(ax=ax_depth)
        ax_depth.tick_params(axis="y", length=0)
        depth_title = "Depth of coverage"
        if depth_error is not None:
            depth_title += f" (avg. norm. diff. {depth_error:.3f})"
        ax_depth.set_title(depth_title)

        # Annotations for count plot (right subplot)
        self.remove_spines(ax=ax_count)
        self.add_horizontal_lines(ax=ax_count)
        ax_count.tick_params(axis="y", length=0)
        count_title = "Read count"
        if count_error is not None:
            count_title += f" (avg. norm. diff. {count_error:.3f})"
        ax_count.set_title(count_title)

        # Set consistent x-axis ticks for both plots
        num_ticks = 5
        ax_depth.xaxis.set_major_locator(mticker.LinearLocator(numticks=num_ticks))
        ax_count.xaxis.set_major_locator(mticker.LinearLocator(numticks=num_ticks))

        # Get handles and labels for legend
        handles, labels = ax_count.get_legend_handles_labels()

        # Create figure-level legend at bottom center, horizontally arranged
        fig = ax_count.figure
        fig.legend(
            handles,
            labels,
            loc="lower center",
            ncol=len(handles),
            frameon=False,
        )

        # Align y=0 lines for plots
        self.align_zero_lines(ax_depth, ax_count)

        # Add region coordinates text at bottom left of figure
        region_info = (
            f"{self.intervals.contig}:{self.intervals.start:,}-{self.intervals.end:,}"
        )
        fig.text(
            0,
            -0.05,
            region_info,
            fontsize=10,
            va="bottom",
            ha="left",
        )

    @staticmethod
    def align_zero_lines(ax_depth, ax_count) -> None:
        """Align y=0 lines between both plots to the same vertical position"""
        # Get current y-limits from both plots
        y_min_line, y_max_depth = ax_depth.get_ylim()
        y_min_bar, y_max_count = ax_count.get_ylim()

        # Calculate ranges below zero for line plot
        below_0_ratio = -y_min_line / (y_max_depth - y_min_line)

        # Set bar plot to match the same range below the zero line as the line plot
        ax_count.set_ylim(-below_0_ratio * (y_max_count - y_min_bar), y_max_count)

    def add_read_count_plot(self, ax):
        """Add filled line plot of read counts in each interval"""
        logger.info("[PLOTTER] - Add filled line plot")

        # Calculate center of each interval
        self.intervals.bed["center"] = (
            self.intervals.bed["start"] + self.intervals.bed["end"]
        ) / 2

        # Plot read counts for each BAM file
        for bam_key in self.bams.keys():

            # Get temporary dataframe with center and counts
            df = self.intervals.bed.copy()

            # Populate read counts for each interval
            df["counts"] = df.apply(
                lambda row: self.get_counts(
                    bam=self.bams[bam_key], start=row["start"], end=row["end"]
                )[0],
                axis=1,
            )

            ax.plot(
                df["center"],
                df["counts"],
                label=bam_key,
                color=self.colormap[bam_key],
                linewidth=1,
            )
            ax.fill_between(
                df["center"], df["counts"], alpha=0.3, color=self.colormap[bam_key]
            )

    def remove_spines(self, ax) -> None:
        """Remove top, right, and left spines, keep only bottom"""
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

    def add_horizontal_lines(self, ax) -> None:
        """Add horizontal lines to the plot"""
        ax.set_axisbelow(True)
        ax.grid(
            True,
            which="major",
            axis="y",
            color="gray",
            alpha=0.5,
            linestyle="-",
            linewidth=1,
        )

    def calculate_normalized_difference(self, map_val: float, out_val: float) -> float:
        """
        Calculate normalized difference: (map - out) / (map + out)

        This metric ranges from -1 to 1, where:
        - 0 means perfect agreement
        - Positive values mean map > out
        - Negative values mean map < out

        Returns 0.0 if both values are 0 (to avoid division by zero).
        """
        denominator = map_val + out_val
        if denominator == 0:
            return 0.0
        return (map_val - out_val) / denominator

    def calculate_depth_error(self) -> float:
        """
        Calculate average normalized difference for depth of coverage.
        Compares map and out BAM files at each position.

        Uses signed differences, allowing positive and negative differences
        to cancel out, which better represents systematic bias vs balanced differences.

        Returns:
            Average normalized difference across all positions (signed, ranges from -1 to 1)
        """
        if "map" not in self.bams or "out" not in self.bams:
            logger.warning(
                "[PLOTTER] - Cannot calculate depth error: map and/or out BAM not provided"
            )
            return None

        logger.info("[PLOTTER] - Calculate depth error")

        # Get pileups for map and out BAMs
        map_pileup = self.map_bam.bam.pileup(
            contig=self.map_bam.normalize_contig(contig=self.intervals.contig),
            start=self.intervals.start,
            end=self.intervals.end,
        )
        out_pileup = self.out_bam.bam.pileup(
            contig=self.out_bam.normalize_contig(contig=self.intervals.contig),
            start=self.intervals.start,
            end=self.intervals.end,
        )

        # Convert to DataFrames for easier comparison
        map_df = pd.DataFrame(
            [(a.reference_pos, a.nsegments) for a in map_pileup],
            columns=["coord", "depth"],
        )
        out_df = pd.DataFrame(
            [(a.reference_pos, a.nsegments) for a in out_pileup],
            columns=["coord", "depth"],
        )

        # Merge on coordinate to align positions
        merged = pd.merge(
            map_df, out_df, on="coord", how="outer", suffixes=("_map", "_out")
        )
        # Fill NaN with 0 for positions present in only one BAM
        merged = merged.fillna(0)

        # Calculate normalized difference for each position
        merged["normalized_diff"] = merged.apply(
            lambda row: self.calculate_normalized_difference(
                row["depth_map"], row["depth_out"]
            ),
            axis=1,
        )

        # Return signed average (allows positive and negative differences to cancel out)
        avg_error = merged["normalized_diff"].mean()

        logger.info(
            f"[PLOTTER] - Depth error (avg normalized diff, signed): {avg_error:.6f}"
        )
        return avg_error

    def calculate_read_count_error(self) -> float:
        """
        Calculate average normalized difference for read counts per interval.
        Compares map and out BAM files for each interval in the BED file.

        Uses signed differences, allowing positive and negative differences
        to cancel out, which better represents systematic bias vs balanced differences.

        Returns:
            Average normalized difference across all intervals (signed, ranges from -1 to 1)
        """
        if "map" not in self.bams or "out" not in self.bams:
            logger.warning(
                "[PLOTTER] - Cannot calculate read count error: map and/or out BAM not provided"
            )
            return None

        logger.info("[PLOTTER] - Calculate read count error")

        # Get counts for each interval
        df = self.intervals.bed.copy()
        df["counts_map"] = df.apply(
            lambda row: self.get_counts(
                bam=self.map_bam, start=row["start"], end=row["end"]
            )[0],
            axis=1,
        )
        df["counts_out"] = df.apply(
            lambda row: self.get_counts(
                bam=self.out_bam, start=row["start"], end=row["end"]
            )[0],
            axis=1,
        )

        # Calculate normalized difference for each interval
        df["normalized_diff"] = df.apply(
            lambda row: self.calculate_normalized_difference(
                row["counts_map"], row["counts_out"]
            ),
            axis=1,
        )

        # Return signed average (allows positive and negative differences to cancel out)
        avg_error = df["normalized_diff"].mean()

        logger.info(
            f"[PLOTTER] - Read count error (avg normalized diff, signed): {avg_error:.6f}"
        )
        return avg_error

    def calculate_and_report_error(self) -> tuple[Optional[float], Optional[float]]:
        """
        Calculate error metrics for both depth and read counts,
        then report the overall average error.

        Uses signed differences, allowing positive and negative differences
        to cancel out, which better represents systematic bias vs balanced differences.

        Returns:
            Tuple of (depth_error, count_error). Either may be None if BAMs are missing.
            Values are signed (range -1 to 1), where:
            - Positive values indicate map > out on average
            - Negative values indicate map < out on average
            - Values near 0 indicate balanced differences
        """
        depth_error = self.calculate_depth_error()
        count_error = self.calculate_read_count_error()

        if depth_error is None or count_error is None:
            logger.warning(
                "[PLOTTER] - Error calculation incomplete: map and/or out BAM not available"
            )
        else:
            # Calculate overall average error (signed, allows cancellation)
            overall_error = (depth_error + count_error) / 2.0

            logger.info(
                f"[PLOTTER] - Error metrics (signed) - Depth: {depth_error:.6f}, "
                f"Read count: {count_error:.6f}, Overall: {overall_error:.6f}"
            )

        return depth_error, count_error
