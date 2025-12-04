import logging
from typing import Optional

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

from samsampleX.Loader import Loader
from samsampleX.FileHandler import FileHandler

logger = logging.getLogger(__name__)


class Plotter(FileHandler):

    def __init__(
        self,
        in_bam: Optional[str],
        map_bam: Optional[str],
        out_bam: Optional[str],
        region: str,
        out_plt: str,
    ) -> None:

        logger.info("[PLOTTER] - Initialize")

        self.bams = {}
        self.colormap = {}
        self.contig, self.region_start, self.region_end = self.parse_region(region)

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

        self.depths = pd.DataFrame()
        self.counts = pd.DataFrame()

        self.depths["coord"] = np.arange(self.region_start, self.region_end)
        self.counts["coord"] = np.arange(self.region_start, self.region_end)

        self.out_plt = out_plt

    @staticmethod
    def parse_region(region: str) -> tuple[str, int, int]:
        if ":" not in region or "-" not in region:
            raise ValueError("Use samtools-style region notation contig:start-end")
        contig_part, coords = region.split(":", 1)
        start_str, end_str = coords.replace(",", "").split("-", 1)
        start = int(start_str)
        end = int(end_str)

        return contig_part, start, end

    def run_plotting(self) -> None:

        logger.info("[PLOTTER] - Begin plotting")

        self.get_depths()
        self.get_counts()

        fig, ax_depth, ax_count = self.setup_plot()
        self.add_depth_plot(ax=ax_depth)
        self.add_read_count_plot(ax=ax_count)

        self.add_annotations(
            ax_depth=ax_depth,
            ax_count=ax_count,
        )

        logger.info("[PLOTTER] - Save plot")
        plt.savefig(self.out_plt, dpi=600, bbox_inches="tight")

    def get_depths(self) -> None:
        """Get depths per position in region"""
        logger.info("[PLOTTER] - Pileup BAMs")

        for b in self.bams:
            depth_array = np.zeros(self.region_end - self.region_start, dtype=int)
            for column in self.bams[b].bam.pileup(
                contig=self.bams[b].normalize_contig(contig=self.contig),
                start=self.region_start,
                end=self.region_end,
                truncate=True,
            ):
                idx = column.reference_pos - self.region_start
                if 0 <= idx < depth_array.size:
                    depth_array[idx] = column.nsegments

            self.depths[b] = depth_array

    def get_counts(self) -> None:
        """Get read counts per position in region"""
        logger.info("[PLOTTER] - Count reads")

        for b in self.bams:
            coverage = self.bams[b].bam.count_coverage(
                self.bams[b].normalize_contig(contig=self.contig),
                self.region_start,
                self.region_end,
            )

            self.counts[b] = np.sum(np.array(coverage), axis=0)

    @staticmethod
    def setup_plot() -> tuple:
        """Set up the figure with two subplots"""
        logger.info("[PLOTTER] - Setup plot")

        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

        fig, (ax_depth, ax_count) = plt.subplots(
            1, 2, figsize=(12, 5), layout="constrained"
        )

        # Formatting for large x-axis values
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

        for b in self.bams:
            ax.plot(
                self.depths["coord"],
                self.depths[b],
                label=b,
                color=self.colormap[b],
                linewidth=1,
            )
            ax.fill_between(
                self.depths["coord"],
                self.depths[b],
                alpha=0.3,
                color=self.colormap[b],
            )

    def add_annotations(
        self,
        ax_depth,
        ax_count,
    ) -> None:

        logger.info("[PLOTTER] - Add axes, title, legend")

        # Annotations for depth plot
        self.remove_spines(ax=ax_depth)
        self.add_horizontal_lines(ax=ax_depth)

        ax_depth.tick_params(axis="y", length=0)
        ax_depth.set_title("Depth of coverage")

        # Annotations for read count plot
        self.remove_spines(ax=ax_count)
        self.add_horizontal_lines(ax=ax_count)
        ax_count.tick_params(axis="y", length=0)
        ax_count.set_title("Read count")

        # Set consistent x-axis ticks across the two plots
        num_ticks = 5
        ax_depth.xaxis.set_major_locator(mticker.LinearLocator(numticks=num_ticks))
        ax_count.xaxis.set_major_locator(mticker.LinearLocator(numticks=num_ticks))

        # Get handles and labels for legend
        handles, labels = ax_count.get_legend_handles_labels()

        # Create figure-level legend at bottom center
        fig = ax_count.figure
        fig.legend(
            handles,
            labels,
            loc="lower center",
            ncol=len(handles),
            frameon=False,
        )

        # Add region coordinates text at bottom left of figure
        region_info = f"{self.contig}:{self.region_start:,}-{self.region_end:,}"
        fig.text(
            0,
            -0.05,
            region_info,
            fontsize=10,
            va="bottom",
            ha="left",
        )

    def add_read_count_plot(self, ax):
        """Add filled line plot of read counts in each interval"""
        logger.info("[PLOTTER] - Add filled line plot")

        # Plot read counts for each BAM file
        for b in self.bams:

            ax.plot(
                np.arange(self.region_start, self.region_end),
                self.counts[b],
                label=b,
                color=self.colormap[b],
                linewidth=1,
            )
            ax.fill_between(
                np.arange(self.region_start, self.region_end),
                self.counts[b],
                alpha=0.3,
                color=self.colormap[b],
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
