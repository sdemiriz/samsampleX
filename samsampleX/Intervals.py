from logging import info

import numpy as np
import pandas as pd
from intervaltree import Interval, IntervalTree

from samsampleX.FileHandler import FileHandler


class Intervals(FileHandler):

    def __init__(self, bed_file: str) -> None:
        """
        Class constructor: read, validate BED and populate IntervalTree
        """
        info(f"Intervals - Initialize Intervals")

        # Read BED file
        self.bed_path = self.get_bed_path(bed_file=bed_file)
        self.bed = self.get_bed(bed_path=self.bed_path)

        # Get contig, start, and end from BED
        self.contig = self.get_contig()
        self.start, self.end = self.get_limits()

        # Build IntervalTree from BED
        self.populate_tree()
        self.validate()

        # Pre-compute interval arrays for overlap detection
        self.interval_starts = np.array([interval.begin for interval in self.tree])
        self.interval_ends = np.array([interval.end for interval in self.tree])

        info(f"Intervals - Complete initialize Intervals")

    def get_bed_path(self, bed_file: str) -> str:
        """
        Read and validate the provided BED file path
        """
        info("Intervals - Handle BED file")

        if not bed_file:
            raise ValueError("No BED file provided")

        super().check_file_exists(path=bed_file)
        info(f"Intervals - Received BED file path {bed_file}")
        return bed_file

    def __len__(self) -> int:
        return len(self.tree)

    def get_contig(self) -> str:
        """
        Gets contig from first BED file (assumes BED files have been validated)
        """
        info(f"Intervals - Set contig value")
        return str(self.bed["contig"][0])

    def get_bed(self, bed_path: str) -> pd.DataFrame:
        """
        Read BED file from supplied filename
        """
        info(f"Intervals - Read BED from {bed_path}")

        return pd.read_csv(  # type: ignore
            filepath_or_buffer=bed_path,
            sep="\t",
            header=None,
            names=["contig", "start", "end", "read_count", "fraction"],
            dtype={
                "contig": str,
                "start": int,
                "end": int,
                "read_count": int,
                "fraction": float,
            },
        )

    def populate_tree(self) -> None:
        """
        Populate IntervalTree using rows from BED file
        """
        info("Intervals - Populate interval tree using mean read_counts")

        tree = IntervalTree()
        for row in self.bed.itertuples(index=False):
            tree.add(Interval(begin=row[1], end=row[2], data=row[3]))

        self.tree = IntervalTree(sorted(tree))
        info(f"Intervals - Created IntervalTree from {len(self.tree)} intervals")

    def get_limits(self) -> tuple[int, int]:
        """
        Return start and end of region from the first BED file (assumes BED files have been validated)
        """
        return min(self.bed["start"]), max(self.bed["end"])

    def validate(self) -> None:
        """
        Validate BED files when
        """
        info(f"Intervals - Validate last BED file and tree")

        contig = ""
        start, end = None, None

        # Check for different contigs within BED files
        assert (
            len(pd.unique(self.bed["contig"])) == 1
        ), f"Not all contig values in BED file identical"

        # Check for different contigs between BED files
        if contig:
            assert contig == self.bed["contig"][0], "Contig differs from other contigs"
            contig = self.bed["contig"][0]

        # Check all start and end cooridnates match between BED files
        if start and end:
            assert start == min(
                self.bed["start"]
            ), "Start coordinate differs from other start coordinates"
            assert end == max(
                self.bed["end"]
            ), "End coordinate differs from other end coordinates"
            start, end = min(self.bed["start"]), max(self.bed["end"])

        # Check for interval overlaps
        before = self.tree.copy()
        self.tree.split_overlaps()
        after = self.tree
        assert before == after, "BED file contains overlapping intervals"
