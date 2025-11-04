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
        hlala_mode: bool,
        genome_build: Optional[str] = None,
    ) -> None:
        """Sampling method for both regular and HLA*LA modes."""
        logger.info(f"Sampler - Begin {'HLA*LA ' if hlala_mode else ' '}sampling")

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
