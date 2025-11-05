import logging
from typing import Optional, Generator

import pysam
import numpy as np

from samsampleX.FileHandler import FileHandler


logger = logging.getLogger(__name__)


class Loader(FileHandler):

    def __init__(
        self,
        bam_path: str,
        template: Optional[pysam.AlignmentFile] = None,
        header: Optional[pysam.AlignmentHeader] = None,
    ) -> None:
        """
        Constructor
        """
        logger.info("Loader - Initialize")
        self.SEQUENCES = "HLA-LA/graphs/PRG_MHC_GRCh38_withIMGT/sequences.txt"

        # Set up paths for BAM file initialization
        self.bam_path = bam_path
        self.template = template
        self.bam = self.load_bam(
            bam_path=self.bam_path, template=self.template, header=header
        )

        # Only index if this is a read-mode BAM file (no template provided)
        if not template:
            self.index_if_needed(bam=self.bam, bam_path=self.bam_path)

        logger.info("Loader - Initialized")

    def load_bam(
        self,
        bam_path: str,
        template: Optional[pysam.AlignmentFile] = None,
        header: Optional[pysam.AlignmentHeader] = None,
    ) -> pysam.AlignmentFile:
        """
        Open in "w" mode if a template has been provided, otherwise open in "r" mode
        """

        # If a template has been provided, open in write mode. Header being provided overrides template
        if template:
            logger.info(f"Loader - Template found: load BAM {bam_path} (writing)")
            if header:
                bam = pysam.AlignmentFile(bam_path, mode="wb", header=header)
            else:
                bam = pysam.AlignmentFile(bam_path, mode="wb", template=template)

        # If no template has been provided, open in read mode
        else:
            logger.info(f"Loader - No template: load BAM {self.bam_path} (reading)")
            bam = pysam.AlignmentFile(self.bam_path, mode="rb")

        return bam

    def index_if_needed(self, bam: pysam.AlignmentFile, bam_path: str) -> None:
        """
        Index the BAM file if it is not already indexed
        """
        # This throws an error for unindexed files so
        # we except and index the file at path instead
        try:
            bam.check_index()
        except ValueError as e:
            logger.warning(f"Loader - Indexing BAM file {bam_path}")
            pysam.samtools.index(bam_path)

    def write_reads(self, bam) -> None:
        """
        Open the output BAM file and write all kept reads
        Sort and index for file for random access in following steps
        """
        logger.info(f"Loader - Write reads to file")

        for r in self.reads:
            bam.bam.write(read=r)

        bam.close()
        self.sort_and_index()

    def normalize_contig(self, contig: str) -> str:
        """
        Handle both chrN and N contig names (other formats not supported)
        """
        logger.info(f"Loader - Normalize contig name {contig}")

        if contig in self.bam.references:
            logger.info(f"Loader - Contig {contig} is as expected")
            return contig

        if contig.startswith("chr"):
            normalized_contig = contig[3:]
        else:
            normalized_contig = "chr" + contig

        if normalized_contig not in self.bam.references:
            raise ValueError(f"Loader - Cannot auto-correct contig name {contig}")

        return normalized_contig

    def add_read_to_bucket(
        self, read: pysam.AlignedSegment, buckets: list[list[pysam.AlignedSegment]]
    ):
        """
        Pre-computed interval arrays are used to find reads that overlap with intervals
        """
        # Read coordinates
        read_start = read.reference_start
        read_end = read.reference_end

        # Vectorized overlap detection using pre-computed arrays
        # Two intervals overlap if: max(start1, start2) < min(end1, end2)
        overlap_condition = np.maximum(read_start, self.interval_starts) < np.minimum(
            read_end, self.interval_ends
        )

        # Get indices where overlap is True
        candidate_bucket_indices = np.where(overlap_condition)[0]

        # Randomly select one bucket (index) to assign the read
        if len(candidate_bucket_indices) > 0:
            selected_bucket_index = np.random.choice(a=candidate_bucket_indices)
            buckets[selected_bucket_index].append(read)

    def fetch(
        self, names: Optional[list[str]] = None
    ) -> Generator[pysam.AlignedSegment, None, None]:
        """
        Yield all mapped reads within limits of BED file
        """
        logger.info("Loader - Fetch mapped reads from supplied region")
        if names is not None:
            for read in self.bam.fetch():
                yield read
        else:
            yield from self.bam.fetch()

    def get_reference_name(self, reference_id: int) -> str:
        """
        Get reference name from reference ID
        """
        return self.bam.get_reference_name(reference_id)

    def close(self) -> None:
        """
        Close BAM file using pysam's internal method and clear caches
        """
        logger.info(f"Loader - Close BAM file")
        self.bam.close()
