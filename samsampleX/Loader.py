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

    def run_sampling(
        self,
        bed_file: str,
        main_seed: int,
        out_bam: str,
        hlala_mode: bool,
        genome_build: Optional[str] = None,
    ) -> None:
        """Sampling method for both regular and HLA*LA modes."""
        logger.info(f"Loader - Begin {'HLA*LA ' if hlala_mode else ' '}sampling")

        # Set up random seed and corresponding output BAM file
        self.main_seed = int(main_seed)
        self.write_out = out_bam
        self.out_bam = Loader(bam_path=out_bam, template=self.bam)

        # Generate user-provided intervals, and a bucket and seed per interval
        self.intervals = self.get_intervals(bed_file=bed_file)
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

        # HLA-LA remaps reads to PRG contigs, here we map them back to chr6
        if hlala_mode:

            header_with_chr6 = self.modify_header(genome_build=genome_build)
            self.out_bam = Loader(
                bam_path=out_bam, template=self.bam, header=header_with_chr6
            )
            self.out_bam.setup_mapback(genome_build=genome_build)

            # PRG-mapped reads from HLA*LA need to be mapped back to chr6
            prg_reads = self.get_reads_from_contigs(contigs=self.get_prg_contigs())
            prg_read_count = sum(
                1 for _ in self.get_reads_from_contigs(contigs=self.get_prg_contigs())
            )

            for r in tqdm(
                iterable=prg_reads,
                desc="Processed",
                unit=" PRG-mapped reads",
                total=prg_read_count,
            ):
                r_chr6 = self.map_read_to_chr6(read=r)
                read_coords = (
                    r_chr6.reference_start,
                    r_chr6.reference_end,
                )

                if self.overlap(
                    read_coords=read_coords,
                    int_coords=(self.interval_starts[i], self.interval_ends[i]),
                ):
                    self.add_read_to_bucket(read=r_chr6, buckets=self.buckets)

        # If we are not dealing with HLA-LA output, no need for mapback
        else:

            for i, interval in enumerate(self.intervals.tree):
                logger.info(
                    f"Loader - Counting reads in interval [{interval.begin}-{interval.end}]"
                )
                mapped_read_count = sum(
                    1
                    for _ in self.get_mapped_reads(
                        start=self.interval_starts[i], end=self.interval_ends[i]
                    )
                )
                logger.info(
                    f"Loader - Sampling {interval.data} reads from interval [{interval.begin}-{interval.end}] with {mapped_read_count} reads"
                )
                self.buckets[i] = np.random.choice(
                    a=np.arange(mapped_read_count), size=interval.data, replace=False
                )

        # After reads have been sorted into buckets, sample from them and write to file
        for i, interval in enumerate(self.intervals.tree):
            logger.info(
                f"Loader - Fetching reads from interval [{interval.begin}-{interval.end}]"
            )
            reads = []
            for j, r in enumerate(
                self.get_mapped_reads(
                    start=self.interval_starts[i], end=self.interval_ends[i]
                )
            ):
                if j in self.buckets[i]:
                    reads.append(r)

            for r in reads:
                self.out_bam.bam.write(read=r)

        self.out_bam.close()
        self.sort_and_index()

    def setup_mapback(self, genome_build: str) -> None:
        """
        Set up HLA*LA-specific variables
        """

        # https://github.com/DiltheyLab/ContigAnalysisScripts/blob/master/fasta2svg.py
        # Define contig names and corresponding alt contig names
        self.contig_names = {
            "apd": "chr6_GL000250v2_alt",
            "cox": "chr6_GL000251v2_alt",
            "dbb": "chr6_GL000252v2_alt",
            "mann": "chr6_GL000253v2_alt",
            "mcf": "chr6_GL000254v2_alt",
            "qbl": "chr6_GL000255v2_alt",
            "ssto": "chr6_GL000256v2_alt",
            "chr6": "chr6",
            "6": "chr6",
        }

        # Coordinates change between genome builds so we handle the two main ones
        if genome_build == "GRCh38":

            # Boundaries of HLA alleles
            # GENCODE v48 (non-lncRNA ones selected)
            self.gene_maps = {
                "A": ("chr6", 29941260, 29949572),
                "B": ("chr6", 31353872, 31367067),
                "C": ("chr6", 31268749, 31272130),
                "DMA": ("chr6", 32948613, 32969094),
                "DMB": ("chr6", 32934629, 32941028),
                "DOA": ("chr6", 33004182, 33009591),
                "DPA1": ("chr6", 33064569, 33080775),
                "DPB1": ("chr6", 33075936, 33089696),
                "DQA1": ("chr6", 32628179, 32647062),
                "DQB1": ("chr6", 32659467, 32668383),
                "DRA": ("chr6", 32439878, 32445046),
                "DRB1": ("chr6", 32577902, 32589848),
                "DRB3": ("chr6_GL000250v2_alt", 3824514, 3837642),
                "DRB4": ("chr6_GL000253v2_alt", 3840435, 3855431),
                "E": ("chr6", 30489509, 30494194),
                "F": ("chr6", 29722775, 29738528),
                "G": ("chr6", 29826967, 29831125),
                "H": ("chr6", 29887752, 29890482),
                "K": ("chr6", 29926459, 29929232),
                "L": ("chr6", 30259625, 30261703),
                "MICA": ("chr6", 31399784, 31415315),
                "MICB": ("chr6", 31494881, 31511124),
                "P": ("chr6", 29800415, 29802425),
                "TAP1": ("chr6", 32845209, 32853816),
                "TAP2": ("chr6", 32821833, 32838739),
                "V": ("chr6", 29792234, 29793136),
            }

            # How alt contigs map to chr6
            # From UCSC Browser
            self.alt_contig_maps = {
                "chr6_GL000250v2_alt": ("chr6", 28734408, 33367716),
                "chr6_GL000251v2_alt": ("chr6", 28510120, 33383765),
                "chr6_GL000252v2_alt": ("chr6", 28734408, 33361299),
                "chr6_GL000253v2_alt": ("chr6", 28734408, 33258200),
                "chr6_GL000254v2_alt": ("chr6", 28734408, 33391865),
                "chr6_GL000255v2_alt": ("chr6", 28734408, 33411973),
                "chr6_GL000256v2_alt": ("chr6", 28691466, 33480577),
                "chr6": ("chr6", 1, 33480577),
            }

        elif genome_build == "GRCh37":

            # Liftover from GRCh38 (DRB3/4 alt contigs not in GRCh37, kept as is)
            self.gene_maps = {
                "A": ("chr6", 29909037, 29917349),
                "B": ("chr6", 31321649, 31334844),
                "C": ("chr6", 31236526, 31239907),
                "DMA": ("chr6", 32916390, 32936871),
                "DMB": ("chr6", 32902406, 32908805),
                "DOA": ("chr6", 32971959, 32977368),
                "DPA1": ("chr6", 33032346, 33048552),
                "DPB1": ("chr6", 33043713, 33057473),
                "DQA1": ("chr6", 32595956, 32614839),
                "DQB1": ("chr6", 32627244, 32636160),
                "DRA": ("chr6", 32407655, 32412823),
                "DRB1": ("chr6", 32545679, 32557625),
                "DRB3": ("chr6_GL000250v2_alt", 3824514, 3837642),
                "DRB4": ("chr6_GL000253v2_alt", 3840435, 3855431),
                "E": ("chr6", 30457286, 30461971),
                "F": ("chr6", 29690552, 29706305),
                "G": ("chr6", 29794744, 29798902),
                "H": ("chr6", 29855529, 29858259),
                "K": ("chr6", 29894236, 29897009),
                "L": ("chr6", 30227402, 30229480),
                "MICA": ("chr6", 31367561, 31383092),
                "MICB": ("chr6", 31462658, 31478901),
                "P": ("chr6", 29768192, 29770202),
                "TAP1": ("chr6", 32812986, 32821593),
                "TAP2": ("chr6", 32789610, 32806516),
                "V": ("chr6", 29760011, 29760913),
            }

            # Liftover from GRCh38
            self.alt_contig_maps = {
                "chr6_GL000250v2_alt": ("chr6", 28702185, 33335493),
                "chr6_GL000251v2_alt": ("chr6", 28477897, 33351542),
                "chr6_GL000252v2_alt": ("chr6", 28702185, 33329076),
                "chr6_GL000253v2_alt": ("chr6", 28702185, 33225977),
                "chr6_GL000254v2_alt": ("chr6", 28702185, 33359642),
                "chr6_GL000255v2_alt": ("chr6", 28702185, 33379750),
                "chr6_GL000256v2_alt": ("chr6", 28659243, 33448354),
                "chr6": ("chr6", 60000, 33448354),
            }

        # If HLA-LA is installed in the same directory, source the sequences.txt file to
        # map PRG-mapped reads back to chr6
        super().check_file_exists(path=self.SEQUENCES)
        self.sequence_txt = pd.read_csv(
            filepath_or_buffer=self.SEQUENCES,
            sep="\t",
            usecols=["Name", "FASTAID"],  # type: ignore
        )

    def get_intervals(self, bed_file: str) -> Intervals:
        """
        Set up Interval instances based on BED-provided coordinates
        """
        logger.info("Loader - Ingest Intervals from BED files")
        return Intervals(bed_file=bed_file)

    def get_interval_seeds(self, main_seed: int) -> np.ndarray:
        """
        Generate a seed per interval provided
        """
        logger.info(f"Loader - Generate random seeds")
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
        for r in self.bam.fetch(
            contig=contig,
            start=start,
            end=end,
        ):
            if r.is_mapped:
                yield r

    def overlap(
        self, read_coords: tuple[int, int], int_coords: tuple[int, int]
    ) -> bool:
        """
        Determine whether the read coordinates overlap the interval coordinates (start-end)
        """
        read_start, read_end = read_coords
        int_start, int_end = int_coords

        return max(read_start, int_start) < min(read_end, int_end)

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

    def sort_and_index(self) -> None:
        """
        Sort, then index file
        """
        logger.info(f"Loader - Sort and index output BAM file")

        with tempfile.NamedTemporaryFile(dir=".", delete=False) as temp_file:
            pysam.sort(self.write_out, "-o", temp_file.name)
            os.rename(src=temp_file.name, dst=self.write_out)
        pysam.index(self.write_out)

    def sample_reads_from_buckets(self):
        """
        Sort reads that overlap with BED intervals into buckets
        """
        logger.info("Loader - Sort reads into buckets")

        for i, bucket in enumerate(self.buckets):
            logger.info(f"Loader - {len(bucket)} read names in bucket {i}")

        self.reads = []
        for bucket, interval, seed in zip(
            self.buckets, self.intervals.tree, self.seeds
        ):
            logger.info(
                f"Loader - Sampling reads in interval [{interval.begin}-{interval.end}]"
            )

            # Count reads that overhang from previous intervals
            overhang_read_count = sum(
                1
                for prev_read in self.reads
                if self.overlap(
                    read_coords=(prev_read.reference_start, prev_read.reference_end),
                    int_coords=(interval.begin, interval.end),
                )
            )

            # Calculate actual amount of reads to sample
            count = int(interval.data) - overhang_read_count

            logger.info(f"Loader - {overhang_read_count} overhang reads")
            logger.info(f"Loader - {count} reads need to be sampled")

            try:
                np.random.seed(seed=seed)
                self.reads.extend(np.random.choice(a=bucket, size=count, replace=False))
            except ValueError as e:
                logger.warning(
                    f"Loader - Not enough reads in bucket for interval [{interval.begin}-{interval.end}]"
                )

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

    # def add_read_to_bucket(
    #     self, read: pysam.AlignedSegment, buckets: list[list[pysam.AlignedSegment]]
    # ):
    #     """
    #     Use IntervalTree to efficiently find overlapping intervals for read assignment
    #     """
    #     # Empty list to store bucket indices instead of whole buckets
    #     candidate_bucket_indices = []
    #     for i, interval in enumerate(self.intervals.tree):
    #         if self.overlap(
    #             read_coords=(read.reference_start, read.reference_end),
    #             int_coords=(interval.begin, interval.end),
    #         ):
    #             candidate_bucket_indices.append(i)

    #     # Randomly select one bucket (index) to assign the read
    #     np.random.seed(seed=self.main_seed)
    #     if candidate_bucket_indices:
    #         selected_bucket_index = np.random.choice(a=candidate_bucket_indices)
    #         buckets[selected_bucket_index].append(read)

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
