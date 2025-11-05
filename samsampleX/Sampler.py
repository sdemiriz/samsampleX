import os
import tempfile
import logging
from typing import Optional, Generator

import pysam
import numpy as np
import pandas as pd
from tqdm import tqdm

from samsampleX.Loader import Loader
from samsampleX.Intervals import Intervals
from samsampleX.FileHandler import FileHandler

logger = logging.getLogger(__name__)


class Sampler(FileHandler):
    """
    Sampler class to sample reads from a BAM file
    """

    def __init__(self, bam_path: str, bed_file: str) -> None:
        """
        Constructor for Sampler class
        """
        logger.info("Sampler - Initialize Sampler")
        self.bam = Loader(bam_path=bam_path, template=None)
        self.intervals = self.get_intervals(bed_file=bed_file)

    def run_sampling(
        self,
        main_seed: int,
        out_bam: str,
    ) -> None:
        """Sampling method for both regular and HLA*LA modes."""
        logger.info(f"Sampler - Begin sampling")

        # Set up random seed and corresponding output BAM file
        self.main_seed = int(main_seed)
        self.write_out = out_bam
        self.out_bam = Loader(bam_path=out_bam, template=self.bam.bam)

        # Generate user-provided intervals, and a bucket and seed per interval
        self.seeds = self.get_interval_seeds(main_seed=self.main_seed)
        self.buckets = self.setup_buckets()

        self.interval_starts = self.intervals.interval_starts
        self.interval_ends = self.intervals.interval_ends

        # # Pad the interval range by 1000bp to account for overhangs
        # overhang = 1000
        # region_start, region_end = self.intervals.start, self.intervals.end
        # interval_range = (
        #     region_start - overhang,
        #     region_end + overhang,
        # )

        prev_reads = []
        for i, interval in enumerate(self.intervals.tree):

            logger.info(f"Sampler - Processing interval {i}")

            overhang_read_count = sum(
                1
                for r in prev_reads
                if self.overlap(
                    read_coords=(r.reference_start, r.reference_end),
                    int_coords=(self.interval_starts[i], self.interval_ends[i]),
                )
            )

            logger.info(f"Sampler - {overhang_read_count} overhang reads")

            # Get mapped reads have not already been accounted for
            mapped_reads = [
                r
                for r in self.get_mapped_reads(
                    start=self.interval_starts[i], end=self.interval_ends[i]
                )
                if r not in prev_reads
            ]

            # Add current interval reads to previous reads
            prev_reads.extend(mapped_reads)

            logger.info(
                f"Sampler - Sampling {interval.data} reads from {len(mapped_reads)} mapped reads"
            )

            # Sample and hold indices for current intervalreads
            self.buckets[i] = np.random.choice(
                a=np.arange(len(mapped_reads) - overhang_read_count),
                size=interval.data,
                replace=False,
            )

            # Get current interval reads and write out if their index is in the bucket
            for j, r in enumerate(mapped_reads):
                if j in self.buckets[i]:
                    self.out_bam.bam.write(read=r)

        self.bam.close()
        self.out_bam.close()
        self.sort_and_index()

    def get_intervals(self, bed_file: str) -> Intervals:
        """
        Set up Interval instances based on BED-provided coordinates
        """
        logger.info("Sampler - Ingest Intervals from BED files")
        return Intervals(bed_file=bed_file)

    def get_interval_seeds(self, main_seed: int) -> np.ndarray:
        """
        Generate a seed per interval provided
        """
        logger.info(f"Sampler - Generate random seeds")
        np.random.seed(seed=main_seed)
        return np.random.randint(low=0, high=1_000_000, size=len(self.intervals))

    def setup_buckets(self) -> list[list[pysam.AlignedSegment]]:
        """
        Get an empty read bucket to sort reads from per interval provided
        """
        logger.info("Loader - Set up an empty bucket per interval")
        return [[] for i in range(len(self.intervals))]

    def get_mapped_reads(
        self, start: int, end: int
    ) -> Generator[pysam.AlignedSegment, None, None]:
        """
        Yield all mapped reads within limits of BED file
        """
        logger.info("Loader - Fetch mapped reads from supplied region")
        contig = self.normalize_contig(contig=self.intervals.contig)
        for r in self.bam.bam.fetch(
            contig=contig,
            start=start,
            end=end,
        ):
            if r.is_mapped:
                yield r

    def normalize_contig(self, contig: str) -> str:
        """
        Handle both chrN and N contig names (other formats not supported)
        """
        logger.info(f"Sampler - Normalize contig name {contig}")

        if contig in self.bam.bam.references:
            logger.info(f"Sampler - Contig {contig} is as expected")
            return contig

        if contig.startswith("chr"):
            normalized_contig = contig[3:]
        else:
            normalized_contig = "chr" + contig

        if normalized_contig not in self.bam.bam.references:
            raise ValueError(f"Sampler - Cannot auto-correct contig name {contig}")

        return normalized_contig

    def overlap(
        self, read_coords: tuple[int, int], int_coords: tuple[int, int]
    ) -> bool:
        """
        Determine whether the read coordinates overlap the interval coordinates (start-end)
        """
        read_start, read_end = read_coords
        int_start, int_end = int_coords

        return max(read_start, int_start) < min(read_end, int_end)

    def sort_and_index(self) -> None:
        """
        Sort, then index file
        """
        logger.info(f"Sampler - Sort and index output BAM file")

        with tempfile.NamedTemporaryFile(dir=".", delete=False) as temp_file:
            pysam.sort(self.write_out, "-o", temp_file.name)
            os.rename(src=temp_file.name, dst=self.write_out)
        pysam.index(self.write_out)
