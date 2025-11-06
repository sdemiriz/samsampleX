import os
import tempfile
import logging
from typing import Optional, Generator

import pysam
import numpy as np
import pandas as pd
from tqdm import tqdm

from samsampleX.Loader import Loader
from samsampleX.Bucket import Bucket
from samsampleX.Intervals import Intervals
from samsampleX.FileHandler import FileHandler

logger = logging.getLogger(__name__)


class Sampler(FileHandler):
    """
    Sample reads from Loader and Intervals instances
    """

    def __init__(
        self, in_bam_path: str, bed_path: str, seed: int, out_bam_path: str
    ) -> None:
        """
        Initialize components necessary for sampling
        """
        logger.info("Sampler - Initialize Sampler")

        logger.info(f"[SAMPLER] - Set up inputs")
        self.target: Loader = Loader(bam_path=in_bam_path, template=None)
        self.template: Intervals = Intervals(bed_path=bed_path)
        self.contig: str = self.normalize_contig(contig=self.template.contig)
        self.interval_count: int = len(self.template)
        self.main_seed: int = seed

        logger.info(f"[SAMPLER] - Set up outputs")
        self.write_out: str = out_bam_path
        self.result: Loader = Loader(bam_path=out_bam_path, template=self.target.bam)

    def run_sampling(self) -> None:
        """Iterate through intervals and sample reads"""
        logger.info(f"[SAMPLER] - Begin sampling")

        self.seeds: np.ndarray = self.get_interval_seeds()
        self.buckets: list[Bucket] = self.get_empty_buckets()

        self.template_starts: np.ndarray = self.template.interval_starts
        self.template_ends: np.ndarray = self.template.interval_ends

        seen_reads: set[tuple[str, int, int]] = set()
        for i, interval in tqdm(
            enumerate(self.template.tree),
            total=len(self.template),
            desc="Sampling intervals",
            unit=" intervals",
        ):

            # Get reads without duplicates
            mapped_reads = [
                r
                for r in self.target.fetch(
                    contig=self.contig,
                    start=self.template_starts[i],
                    end=self.template_ends[i],
                )
                if (r.query_name, r.reference_start, r.reference_end) not in seen_reads
            ]

            read_count = len(mapped_reads)
            if read_count > 0:
                np.random.seed(self.seeds[i])
                chosen_reads = np.random.choice(
                    a=np.arange(len(mapped_reads)),
                    size=interval.data,
                    replace=False,
                )

                reads = [mapped_reads[j] for j in chosen_reads]
                self.buckets[i].add_read(reads)
                seen_reads.update(
                    [(r.query_name, r.reference_start, r.reference_end) for r in reads]
                )

            logger.info(f"[SAMPLER] - Interval {i}: {len(mapped_reads)} reads found")
            logger.info(f"[SAMPLER] - Interval {i}: {interval.data} reads requested")
            logger.info(
                f"[SAMPLER] - Interval {i}: {len(self.buckets[i])} reads sampled"
            )

            seen_reads.clear()
            self.buckets[i].write_reads()

        # Clean up file I/O
        self.target.close()
        self.result.close()
        self.result.sort_and_index()

    def get_interval_seeds(self) -> np.ndarray:
        """
        Use main seed value to generate independent random seeds per interval
        """
        logger.info(f"[SAMPLER] - Generate random seeds")
        np.random.seed(self.main_seed)
        return np.random.randint(low=0, high=100_000_000, size=self.interval_count)

    def get_empty_buckets(self) -> list[Bucket]:
        """
        Generate an empty bucket per interval
        """
        logger.info("[SAMPLER] - Set up an empty bucket per interval")
        return [Bucket(out_bam=self.result) for _ in range(self.interval_count)]

    def normalize_contig(self, contig: str) -> str:
        """
        Handle both chrN and N contig names (other formats not supported)
        """
        logger.info(f"[SAMPLER] - Normalize contig name {contig}")

        if contig in self.target.bam.references:
            logger.info(f"[SAMPLER] - Contig {contig} found in BAM header")
            return contig

        if contig.startswith("chr"):
            norm_contig: str = contig[3:]
        else:
            norm_contig: str = "chr" + contig

        if norm_contig not in self.target.bam.references:
            raise ValueError(f"Cannot normalize contig name {contig}")

        logger.info(f"[SAMPLER] - Normalized contig name {contig} to {norm_contig}")
        return norm_contig

    def overlap(
        self, read_coords: tuple[int, int], int_coords: tuple[int, int]
    ) -> bool:
        """
        Check if read overlaps interval
        """
        read_start, read_end = read_coords
        int_start, int_end = int_coords

        return max(read_start, int_start) < min(read_end, int_end)
