import logging
from typing import Optional

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from subsample_reads.Loader import Loader
from subsample_reads.Intervals import Intervals
from subsample_reads.FileHandler import FileHandler

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

        self.intervals = Intervals(bed_file=bed)
        self.boundaries = set(
            list(self.intervals.bed["start"]) + list(self.intervals.bed["end"])
        )

        self.bams = {}
        self.colormap = {}
        if in_bam:
            self.in_bam = Loader(bam_path=in_bam)
            self.bams["in"] = self.in_bam
            self.colormap["in"] = "#648FFF"
        if map_bam:
            self.map_bam = Loader(bam_path=map_bam)
            self.bams["map"] = self.map_bam
            self.colormap["map"] = "#DC267F"
        if out_bam:
            self.out_bam = Loader(bam_path=out_bam)
            self.bams["out"] = self.out_bam
            self.colormap["out"] = "#FFB000"

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
            self.add_boundaries(ax=ax_line)
            self.add_fractions(ax=ax_line)
            self.add_boundaries(ax=ax_bar)

        logger.info("Plotter - Save plot")
        plt.savefig(self.out_plt, dpi=600)

    @staticmethod
    def setup_plot() -> tuple:
        """Set up the figure with two separate subplots side by side"""
        logger.info("Plotter - Setup plot")
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
            ax.plot(
                pileup["coord"], pileup["depth"], label=l, alpha=0.8, color=c, zorder=1
            )

    def add_boundaries(self, ax) -> None:
        """Add vertical lines signifying interval boundaries to supplied axis"""
        logger.info("Plotter - Draw boundaries")

        for b in self.boundaries:
            ax.axvline(
                x=b, linestyle="--", linewidth=1.5, color="gray", alpha=0.5, zorder=1
            )

    def add_annotations(self, ax_line, ax_bar) -> None:
        """Add various text labels to supplied axes"""
        logger.info("Plotter - Add axes, title, legend")

        # Annotations for line plot (left subplot)
        ax_line.set_title("Coverage Depth")
        ax_line.set_xlabel("Chromosomal coordinate")
        ax_line.set_ylabel("Depth of coverage")
        ax_line.margins(y=0.1)

        # Annotations for bar plot (right subplot)
        ax_bar.set_title("Read Counts per Interval")
        ax_bar.set_xlabel("Chromosomal coordinate")
        ax_bar.set_ylabel("Read count")
        ax_bar.margins(y=0.1)

        ax_line.legend(loc="upper right")

        # Add super title to the figure
        region_info = (
            f"{self.intervals.contig}:{self.intervals.start}-{self.intervals.end}"
        )
        ax_line.figure.suptitle(f"{region_info}", fontsize=14)

    def add_fractions(self, ax):
        """Add fractions to represent distribution values to supplied axis"""
        logger.info("Plotter - Add fraction annotations")

        y = ax.get_ylim()[1] * 0.97
        ha = "center"
        color = "black"
        size = "x-small"
        zorder = 2

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
                    zorder=zorder,
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
                    zorder=zorder,
                )

    def add_barplot(self, ax):
        """Add barplot signifying number of reads in each interval"""
        logger.info("Plotter - Add barplot")

        for row in self.intervals.bed.iterrows():
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

            for i, c in zip(range(bam_count), self.colormap.values()):
                ax.bar(
                    x=(start + end) / 2 - (offset - i) * width / bam_count,
                    height=counts[i],
                    width=width / bam_count,
                    color=c,
                    alpha=0.5,
                    zorder=0,
                )
