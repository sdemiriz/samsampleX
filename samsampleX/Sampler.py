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
        self.seeds: np.ndarray = self.get_interval_seeds()

        logger.info(f"[SAMPLER] - Set up output")
        self.result: Loader = Loader(bam_path=out_bam_path, template=self.target.bam)

    def run_sampling(self) -> None:
        """Iterate through intervals and sample reads"""
        logger.info(f"[SAMPLER] - Begin sampling")

        self.buckets = sorted(
            [
                {
                    "reads": [],
                    "target": t.data,
                    "start": t.begin,
                    "end": t.end,
                    "interval_idx": i,
                }
                for i, t in enumerate(self.template.tree)
                if t.data > 0
            ],
            key=lambda x: x["target"],
        )

        # Pass 1: Sort read indices into buckets
        for i, r in enumerate(
            self.target.fetch(
                contig=self.contig, start=self.template.start, end=self.template.end
            )
        ):
            for bucket in self.buckets:
                if self.overlap(
                    read_coords=(r.reference_start, r.reference_end),
                    int_coords=(bucket["start"], bucket["end"]),
                ):
                    bucket["reads"].append(i)

        self.selected_reads = []
        logger.info(f"[SAMPLER] - Pass 1 completed")

        # Process buckets in ascending order of target value
        while len(self.buckets) > 0:
            # Sort buckets, smallest to largest target value
            self.buckets = sorted(self.buckets, key=lambda x: x["target"])
            self.buckets = [b for b in self.buckets if len(b["reads"]) > 0]

            # Check if the first bucket has any reads and target > 0
            if len(self.buckets[0]["reads"]) == 0 or self.buckets[0]["target"] <= 0:

                # Remove leftover reads because they should not be sampled again
                for b in self.buckets:
                    b["reads"] = [
                        r for r in b["reads"] if r not in self.buckets.pop(0)["reads"]
                    ]
                continue

            # Sample a read from the smallest target bucket
            np.random.seed(self.seeds[self.buckets[0]["interval_idx"]])
            sampled_read_idx = int(
                np.random.choice(a=self.buckets[0]["reads"], size=1, replace=False)[0]
            )

            # Find all buckets with this read index
            containing_buckets = [
                b for b in self.buckets if sampled_read_idx in b["reads"]
            ]

            # Drop sampled index from all buckets and decrement targets
            buckets_to_remove = []
            for b in containing_buckets:
                b["reads"].remove(sampled_read_idx)
                b["target"] = b["target"] - 1
                # Mark buckets that should be removed (empty or target reached zero)
                if len(b["reads"]) == 0 or b["target"] <= 0:
                    buckets_to_remove.append(b)

            # Remove exhausted buckets and clean up their reads from remaining buckets
            for b_to_remove in buckets_to_remove:
                leftover_reads = b_to_remove["reads"]
                self.buckets.remove(b_to_remove)
                # Remove leftover reads from all remaining buckets
                for b in self.buckets:
                    b["reads"] = [r for r in b["reads"] if r not in leftover_reads]

            # Record the sampled read
            self.selected_reads.append(sampled_read_idx)

        logger.info(f"[SAMPLER] - Finish sampling")

        # Pass 2: Write reads with selected indices to output BAM
        for i, r in enumerate(
            self.target.fetch(
                contig=self.contig, start=self.template.start, end=self.template.end
            )
        ):
            if i in self.selected_reads:
                self.result.bam.write(r)

        logger.info(f"[SAMPLER] - Pass 2 completed")

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
