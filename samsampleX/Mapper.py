import bisect
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from tqdm import tqdm

from samsampleX.Loader import Loader
from samsampleX.FileHandler import FileHandler

logger = logging.getLogger(__name__)


class Mapper(FileHandler):

    def __init__(
        self,
        bam_paths: list[str],
        contig: str,
        start: str,
        end: str,
        interval_length: str = str(10),
        bed_dir: str = "bed/",
        bed: Optional[list[str]] = None,
    ) -> None:
        logger.info("[MAPPER] - Initialize Mapper")

        self.bam_paths = bam_paths
        self.bed = bed

        if bed is not None:
            self.bed_paths = [Path(bed_file) for bed_file in bed]
            for bed_path in self.bed_paths:
                bed_path.parent.mkdir(parents=True, exist_ok=True)
        else:
            self.make_bed_dir(path=bed_dir)
            self.bed_paths = self.get_bed_paths(bed_dir=self.bed_dir)

        self.bams = self.load_bams(bam_paths=self.bam_paths)
        self.contig = contig
        self.beds = self.make_beds(
            contig=self.contig,
            start=start,
            end=end,
            interval_length=interval_length,
        )

        self.populate_read_counts(contig=self.contig, beds=self.beds, bams=self.bams)
        self.write_beds(beds=self.beds, bed_paths=self.bed_paths)

    def make_bed_dir(self, path: str) -> None:
        """
        Make the target directory to place BED files into
        """
        self.bed_dir = Path(path)
        self.bed_dir.mkdir(parents=True, exist_ok=True)

    def get_bed_paths(self, bed_dir: Path) -> list[Path]:
        """
        Generate output BED filenames from BAM filenames
        """
        return [
            bed_dir / Path(bam_path).with_suffix(".bed").name
            for bam_path in self.bam_paths
        ]

    def load_bams(self, bam_paths: list[str]) -> list[Loader]:
        """
        Initialize Loaders for all supplied BAMs
        """
        logger.info(f"[MAPPER] - Initialize Loaders for supplied BAM files")
        return [Loader(bam_path=path) for path in bam_paths]

    def make_beds(
        self,
        contig: str,
        start: str,
        end: str,
        interval_length: str,
    ) -> list[pd.DataFrame]:
        """
        Construct BED-formatted DataFrame
        """
        logger.info("[MAPPER] - Form intervals for BED file")
        interval_boundaries = self.get_interval_boundaries(
            start=int(start),
            end=int(end),
            interval_length=int(interval_length),
        )

        bed_columns = ["contig", "start", "end", "read_count"]

        beds = []
        for bam in self.bams:
            intervals = []
            for s, e in zip(interval_boundaries[:-1], interval_boundaries[1:]):
                intervals.append(
                    {
                        bed_columns[0]: contig,
                        bed_columns[1]: s,
                        bed_columns[2]: e,
                        bed_columns[3]: -1,
                    }
                )

            bed = pd.DataFrame.from_records(
                intervals,
                columns=bed_columns,
            )

            beds.append(bed)

        return beds

    def get_interval_boundaries(
        self,
        start: int,
        end: int,
        interval_length: int,
    ) -> list[int]:
        """
        Divide region based on interval size
        """
        region_length = end - start

        logger.info("[MAPPER] - Use interval size to set up intervals")

        interval_boundaries = [
            i + start for i in range(0, region_length, interval_length)
        ]
        if interval_boundaries[-1] != end:
            interval_boundaries.append(end)

        return interval_boundaries

    def populate_read_counts(
        self, contig: str, beds: list[pd.DataFrame], bams: list[Loader]
    ) -> None:
        """
        Fill read counts in all BED DataFrames using binary search for efficient
        overlap detection. For each read, uses binary search to find candidate
        intervals that could overlap, then checks only those candidates.
        This reduces complexity from O(reads × intervals) to O(reads × log(intervals)).
        """
        logger.info(f"[MAPPER] - Populate read counts in BED files")
        for bed, bam in zip(beds, bams):
            # Pre-compute interval boundaries as lists for binary search
            # Intervals are already sorted by construction
            interval_starts = bed["start"].tolist()
            interval_ends = bed["end"].tolist()
            num_intervals = len(interval_starts)

            # Get the full region to fetch (first interval start to last interval end)
            region_start = int(bed["start"].min())
            region_end = int(bed["end"].max())

            # Initialize read counts to zero
            read_counts = np.zeros(num_intervals, dtype=int)

            # Fetch all reads in the region once
            logger.info(
                f"[MAPPER] - Fetching reads from {contig}:{region_start}-{region_end}"
            )
            for read in tqdm(
                bam.fetch(contig=contig, start=region_start, end=region_end),
                desc=f"Processing reads from {Path(bam.bam_path).name}",
                unit=" reads",
            ):
                # Get read coordinates
                read_start = read.reference_start
                read_end = read.reference_end

                # Use binary search to find candidate intervals efficiently
                # A read overlaps an interval if: read_end > interval_start AND read_start < interval_end
                # Since intervals are sorted and non-overlapping, we can narrow the search space

                # Find first interval that cannot overlap (starts at or after read_end)
                # All intervals before this are candidates (they start before read_end)
                first_non_candidate = bisect.bisect_left(interval_starts, read_end)

                # Find first interval that could overlap (ends after read_start)
                # All intervals from this point forward are candidates (they end after read_start)
                first_candidate = bisect.bisect_right(interval_ends, read_start)

                # Candidate intervals are in range [first_candidate, first_non_candidate)
                # This narrows from checking all intervals to typically 1-2 intervals
                start_idx = first_candidate
                end_idx = min(first_non_candidate, num_intervals)

                # Check overlap for candidate intervals only
                # Two intervals overlap if: max(start1, start2) < min(end1, end2)
                for i in range(start_idx, end_idx):
                    if max(read_start, interval_starts[i]) < min(
                        read_end, interval_ends[i]
                    ):
                        read_counts[i] += 1

            # Assign counts back to DataFrame
            bed["read_count"] = read_counts

    def write_beds(self, beds: list[pd.DataFrame], bed_paths: list[Path]) -> None:
        """
        Write output to BED file
        """
        logger.info("[MAPPER] - Write BED contents to file")
        for bed, path in zip(beds, bed_paths):
            bed.to_csv(path, sep="\t", index=False, header=False)
            super().check_file_exists(path=str(path))
