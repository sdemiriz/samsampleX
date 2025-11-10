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
        no_det: bool = False,
    ) -> None:
        """
        Initialize Plotter with argument names.

        Args:
            in_bam:   Path to input/original BAM file.
            map_bam:  Path to mapping BAM file.
            out_bam:  Path to downsampled BAM file.
            bed:      Path to BED file to plot.
            out_plt:  Path for the output plot.
            no_det:   Whether to include details in plots (disable for large regions)
        """
        logger.info("Plotter - Initialize")

        self.intervals = Intervals(bed_path=bed)
        self.boundaries = set(
            list(self.intervals.bed["start"]) + list(self.intervals.bed["end"])
        )

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
        self.no_det = no_det

        logger.info("Plotter - Complete initialization")

    def get_pileups(self) -> list:
        """Pileup BAMs for the defined region"""
        logger.info("Plotter - Pileup BAMs")

        pileups = [
            bam.bam.pileup(
                contig=bam.normalize_contig(contig=self.intervals.contig),
                start=self.intervals.start,
                end=self.intervals.end,
            )
            for bam in self.bams.values()
        ]

        logger.info("Plotter - Complete pileup BAMs")
        return pileups

    def get_counts(self, start: int, end: int) -> list:
        """Get read counts from BAMs for the defined region"""
        logger.info("Plotter - Pileup BAMs")

        counts = [
            bam.bam.count(
                contig=bam.normalize_contig(self.intervals.contig),
                start=start,
                end=end,
            )
            for bam in self.bams.values()
        ]

        logger.info("Plotter - Complete pileup BAMs")
        return list(counts)

    def plot(self) -> None:
        """Plot provided BAM file pileups and read counts in separate subplots"""
        logger.info("Plotter - Begin plotting")
        fig, ax_line, ax_bar = self.setup_plot()

        # Add content to bar and line plot
        self.add_lineplot(ax=ax_line)
        self.add_barplot(ax=ax_bar)
        self.add_annotations(ax_line=ax_line, ax_bar=ax_bar)

        if not self.no_det:

            self.add_fractions(ax=ax_line)
            self.add_fractions(ax=ax_bar)

        logger.info("Plotter - Save plot")
        plt.savefig(self.out_plt, dpi=600, bbox_inches="tight")

    @staticmethod
    def setup_plot() -> tuple:
        """Set up the figure with two separate subplots side by side"""
        logger.info("Plotter - Setup plot")

        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

        fig, (ax_line, ax_bar) = plt.subplots(
            1, 2, figsize=(12, 5), layout="constrained"
        )

        def custom_formatter(x, pos):
            if x >= 1e6:
                return f"{x/1e6:.3f}m"  # Millions
            else:
                return f"{x:.0f}"

        formatter = mticker.FuncFormatter(custom_formatter)

        # Apply formatting to both subplots
        ax_line.yaxis.set_major_formatter(formatter)
        ax_line.xaxis.set_major_formatter(formatter)
        ax_bar.yaxis.set_major_formatter(formatter)
        ax_bar.xaxis.set_major_formatter(formatter)

        return fig, ax_line, ax_bar

    def add_lineplot(self, ax) -> None:
        """Add line plot representing coverage depth to supplied axis"""
        logger.info("Plotter - Plot pileups")

        for p, l, c in zip(
            self.get_pileups(), self.bams.keys(), self.colormap.values()
        ):
            pileup = pd.DataFrame(
                [(a.reference_pos, a.nsegments) for a in p], columns=["coord", "depth"]
            )
            ax.plot(pileup["coord"], pileup["depth"], label=l, color=c)

    def add_annotations(self, ax_line, ax_bar) -> None:
        """Add various text labels to supplied axes"""
        logger.info("Plotter - Add axes, title, legend")

        # Annotations for line plot (left subplot)
        self.remove_spines(ax=ax_line)
        self.add_horizontal_lines(ax=ax_line)
        self.add_minor_ticks(ax=ax_line)
        ax_line.tick_params(axis="y", length=0)
        ax_line.set_title("Depth of coverage")

        # Annotations for bar plot (right subplot)
        self.remove_spines(ax=ax_bar)
        self.add_horizontal_lines(ax=ax_bar)
        self.add_minor_ticks(ax=ax_bar)
        ax_bar.tick_params(axis="y", length=0)
        ax_bar.set_title("Per-interval read count")

        # Get handles and labels from the bar plot for the legend
        handles, labels = ax_bar.get_legend_handles_labels()

        # Create figure-level legend at bottom center, horizontally arranged
        fig = ax_bar.figure
        fig.legend(
            handles,
            labels,
            loc="lower center",
            ncol=len(handles),
            frameon=False,
        )

        # Align y=0 lines between both plots
        self.align_zero_lines(ax_line, ax_bar)

        # Add region coordinates text at bottom left of figure
        region_info = f"Region {self.intervals.contig}:{self.intervals.start:,}-{self.intervals.end:,}"
        fig.text(
            0,
            -0.05,
            region_info,
            fontsize=10,
            va="bottom",
            ha="left",
        )

    @staticmethod
    def align_zero_lines(ax_line, ax_bar) -> None:
        """Align y=0 lines between both plots to the same vertical position, preserving line plot max"""
        # Get current y-limits from both plots
        y_min_line, y_max_line = ax_line.get_ylim()
        y_min_bar, y_max_bar = ax_bar.get_ylim()

        # Calculate ranges below zero for line plot
        below_0_ratio = -y_min_line / (y_max_line - y_min_line)

        # Set bar plot to match the same range below the zero line as the line plot
        ax_bar.set_ylim(-below_0_ratio * (y_max_bar - y_min_bar), y_max_bar)

    def add_fractions(self, ax):
        """Add fractions to represent distribution values to supplied axis"""
        logger.info("Plotter - Add fraction annotations")

        # Position at bottom of plot, slightly above the bottom line
        y_min, y_max = ax.get_ylim()
        y = y_min + (y_max - y_min) * 0.01
        ha = "center"
        color = "black"
        size = "x-small"

        for row in self.intervals.bed.iterrows():

            x = (row[1]["start"] + row[1]["end"]) / 2
            try:
                ax.text(
                    x=(row[1]["start"] + row[1]["end"]) / 2,
                    y=y,
                    s=str(row[1]["read_count"] / sum(self.intervals.bed["read_count"]))[
                        :5
                    ],
                    ha=ha,
                    color=color,
                    size=size,
                )
            except ZeroDivisionError:
                logger.info(f"Plotter - An interval contains zero reads")
                ax.text(
                    x=x,
                    y=y,
                    s="0",
                    ha=ha,
                    color=color,
                    size=size,
                )

    def add_barplot(self, ax):
        """Add barplot signifying number of reads in each interval"""
        logger.info("Plotter - Add barplot")

        for row_idx, row in enumerate(self.intervals.bed.iterrows()):
            start, end = row[1]["start"], row[1]["end"]
            width = end - start
            counts = self.get_counts(start=start, end=end)

            bam_count = len(self.bams)
            if bam_count == 1:
                offset = 0
            elif bam_count == 2:
                offset = 0.5
            elif bam_count == 3:
                offset = 1
            else:
                raise ValueError(f"Invalid number of BAMs: {bam_count}")

            for i, (bam_key, c) in enumerate(
                zip(self.bams.keys(), self.colormap.values())
            ):
                # Add label only for the first interval to avoid duplicate legend entries
                label = bam_key if row_idx == 0 else None
                ax.bar(
                    x=(start + end) / 2 - (offset - i) * width / bam_count,
                    height=counts[i],
                    width=width / bam_count,
                    color=c,
                    label=label,
                )

    def remove_spines(self, ax) -> None:
        """Remove top, right, and left spines, keep only bottom"""
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

    def add_horizontal_lines(self, ax) -> None:
        """Add horizontal lines to the plot"""
        # Ensure grid lines are drawn behind all plot elements
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

    def add_minor_ticks(self, ax) -> None:
        ax.set_xticks(list(self.boundaries), minor=True)
        ax.tick_params(axis="x", which="minor", length=5, direction="in")
