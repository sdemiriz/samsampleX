import os
import tempfile
import logging
from typing import Optional, Generator

import pysam
import numpy as np
import pandas as pd
from tqdm import tqdm

from samsampleX.Intervals import Intervals
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
        # # we except and index the file at path instead
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

    # def sample_reads_from_buckets(self):
    #     """
    #     Sort reads that overlap with BED intervals into buckets
    #     """
    #     logger.info("Loader - Sort reads into buckets")

    #     for i, bucket in enumerate(self.buckets):
    #         logger.info(f"Loader - {len(bucket)} read names in bucket {i}")

    #     self.reads = []
    #     for bucket, interval, seed in zip(
    #         self.buckets, self.intervals.tree, self.seeds
    #     ):
    #         logger.info(
    #             f"Loader - Sampling reads in interval [{interval.begin}-{interval.end}]"
    #         )

    #         # Count reads that overhang from previous intervals
    #         overhang_read_count = sum(
    #             1
    #             for prev_read in self.reads
    #             if self.overlap(
    #                 read_coords=(prev_read.reference_start, prev_read.reference_end),
    #                 int_coords=(interval.begin, interval.end),
    #             )
    #         )

    #         # Calculate actual amount of reads to sample
    #         count = int(interval.data) - overhang_read_count

    #         logger.info(f"Loader - {overhang_read_count} overhang reads")
    #         logger.info(f"Loader - {count} reads need to be sampled")

    #         try:
    #             np.random.seed(seed=seed)
    #             self.reads.extend(np.random.choice(a=bucket, size=count, replace=False))
    #         except ValueError as e:
    #             logger.warning(
    #                 f"Loader - Not enough reads in bucket for interval [{interval.begin}-{interval.end}]"
    #             )

    def map_read_to_chr6(self, read: pysam.AlignedSegment) -> pysam.AlignedSegment:
        """
        Map provided PRG-mapped read back to chr6 and corresponding coordinates
        """
        # Convert to dictionary to modify
        r_dict = read.to_dict()

        old_ref_name = r_dict["ref_name"]

        # Set up reference names and TIDs as variables
        chr6_tid = len(self.out_bam.bam.header["SQ"]) - 1

        # If read is already mapped to chr6, only switch the TID to chr6's index
        if r_dict["ref_name"] in ["6", "chr6"]:
            r_dict["ref_name"] = "chr6"
            r_dict["tid"] = chr6_tid

        # If read is mapped to PRG contig, convert coordinates to chr6
        elif r_dict["ref_name"].startswith("PRG"):

            # Convert read position to chr6 coordinate, and update TID
            r_dict["ref_pos"] = self.out_bam.convert_to_chr6(
                prg_contig=r_dict["ref_name"], prg_coord=r_dict["ref_pos"]
            )
            r_dict["ref_name"] = "chr6"
            r_dict["tid"] = chr6_tid

        else:
            raise ValueError(f"Loader - Cannot map contig {old_ref_name} to chr6")

        if r_dict["next_ref_name"] == "=":
            r_dict["next_ref_name"] = old_ref_name

        if r_dict["next_ref_name"] in ["6", "chr6"]:
            r_dict["next_ref_name"] = "chr6"
            r_dict["next_tid"] = chr6_tid

        elif r_dict["next_ref_name"].startswith("PRG"):
            r_dict["next_ref_pos"] = self.out_bam.convert_to_chr6(
                prg_contig=r_dict["next_ref_name"], prg_coord=r_dict["next_ref_pos"]
            )
            r_dict["next_ref_name"] = "chr6"
            r_dict["next_tid"] = chr6_tid

        else:
            raise ValueError(f"Loader - Cannot map contig {r_dict['ref_name']} to chr6")

        return pysam.AlignedSegment.from_dict(
            sam_dict=r_dict, header=self.out_bam.bam.header
        )

    def get_reads_from_contigs(
        self, contigs: list[str]
    ) -> Generator[pysam.AlignedSegment, None, None]:
        """
        Yield reads mapped to PRG contigs
        """
        for contig in contigs:
            for r in self.bam.fetch(contig=contig):
                if r.is_mapped:
                    yield r

    def get_prg_contigs(self) -> list[str]:
        """
        Get contigs with names that begin with PRG
        """
        contigs = [
            contig
            for contig in self.bam.references
            if contig.startswith(("PRG", "chr6", "6"))
        ]
        return contigs

    def convert_to_chr6(self, prg_contig: str, prg_coord: int) -> str:
        """
        Calculate chr6-based start coordinate of a PRG-mapped read
        """
        chr6_coords = self.get_prg_coords(prg_contig=prg_contig)
        return str((chr6_coords[0]) + int(prg_coord))

    def get_prg_coords(self, prg_contig: str) -> tuple[int, int]:
        """
        Return and cache (start, end) coordinates of PRG contig
        """
        try:
            prg_contig_mask = self.sequence_txt["FASTAID"] == prg_contig
            sequence_id = self.sequence_txt[prg_contig_mask]["Name"].values[0]
            gene_name = sequence_id.split("*")[0]
        except IndexError:
            raise ValueError(f"Loader - Contig {prg_contig} not found in sequences.txt")

        # If PRG corresponds to HLA allele (with * notation in name)
        if gene_name in self.gene_maps:
            chrom, start, end = self.gene_maps[gene_name]
        # If PRG corresponds to alt contig
        elif sequence_id in self.contig_names:
            alt_contig = self.contig_names[sequence_id]
            chrom, start, end = self.alt_contig_maps[alt_contig]

        return (int(start), int(end))

    def modify_header(self, genome_build: str):
        """
        Add chr6 entry to the SQ section of the input BAM file for read back-mapping
        """
        # Get dictionary of BAM header and add chr6 entry, based on genome build
        header_sq = self.bam.header.to_dict()["SQ"]

        if "chr6" not in [h["SN"] for h in header_sq]:
            if genome_build == "GRCh38":
                header_sq.append({"SN": "chr6", "LN": 170805979})
            elif genome_build == "GRCh37":
                header_sq.append({"SN": "chr6", "LN": 171115067})

        # Replace the original header from BAM with modified one
        header = self.bam.header.to_dict()
        header["SQ"] = header_sq

        # Return header as object
        return pysam.AlignmentHeader.from_dict(header_dict=header)

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
