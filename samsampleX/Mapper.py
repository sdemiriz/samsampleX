import logging
from pathlib import Path
from typing import Optional

import pandas as pd

from samsampleX.Loader import Loader
from samsampleX.FileHandler import FileHandler

logger = logging.getLogger(__name__)


class Mapper(FileHandler):
    """
    Chart the distribution of the provided BAM file
    """

    def __init__(
        self,
        bam_paths: list[str],
        contig: str,
        start: str,
        end: str,
        interval_length: Optional[str] = None,
        interval_count: Optional[str] = None,
        bed_dir: str = "bed/",
        bed: Optional[list[str]] = None,
    ) -> None:
        """
        Constructor for class
        """
        logger.info("Mapper - Initialize Mapper")

        self.bam_paths = bam_paths
        self.bed = bed

        if bed is not None:
            self.bed_paths = [Path(bed_file) for bed_file in bed]
            for bed_path in self.bed_paths:
                bed_path.parent.mkdir(parents=True, exist_ok=True)
        else:
            self.make_bed_dir(path=bed_dir)
            self.bed_paths = self.get_bed_paths(bed_dir=self.bed_dir)

        self.interval_length, self.interval_count = self.setup_intervals(
            interval_length=interval_length,
            interval_count=interval_count,
        )

        self.bams = self.load_bams(bam_paths=self.bam_paths)

        self.contig = contig
        self.beds = self.make_beds(
            contig=self.contig,
            start=start,
            end=end,
            interval_length=self.interval_length,
            interval_count=self.interval_count,
        )

        self.populate_read_counts(contig=self.contig, beds=self.beds, bams=self.bams)
        # self.populate_fractions(contig=self.contig, beds=self.beds, bams=self.bams)

        self.write_beds(beds=self.beds, bed_paths=self.bed_paths)

    def make_bed_dir(self, path: str) -> None:
        """
        Make the target directory to place BED files into
        """
        self.bed_dir = Path(path)
        self.bed_dir.mkdir(parents=True, exist_ok=True)

    def get_bed_paths(self, bed_dir: Path) -> list[Path]:
        """
        Generate output BED filenames from  BAM filenames
        """
        return [
            bed_dir / Path(bam_path).with_suffix(".bed").name
            for bam_path in self.bam_paths
        ]

    def setup_intervals(
        self, interval_length: Optional[str], interval_count: Optional[str]
    ) -> tuple[Optional[int], Optional[int]]:
        """
        Settle whether to use interval length or count when mapping
        """
        if interval_length:
            length = int(interval_length)
            count = None
        elif interval_count:
            length = None
            count = int(interval_count)
        else:
            e = "Mapper - No interval length or count provided"
            logger.error(e)
            raise ValueError(e)

        return length, count

    def load_bams(self, bam_paths: list[str]) -> list[Loader]:
        """
        Initialize Loaders for all supplied BAMs
        """
        logger.info(f"Mapper - Initialize Loaders for supplied BAM files")
        return [Loader(bam_path=path) for path in bam_paths]

    def make_beds(
        self,
        contig: str,
        start: str,
        end: str,
        interval_length: Optional[int] = None,
        interval_count: Optional[int] = None,
    ) -> list[pd.DataFrame]:
        """
        Construct BED-formatted DataFrame
        """
        logger.info("Mapper - Form intervals for BED file")
        interval_boundaries = self.get_interval_boundaries(
            start=int(start),
            end=int(end),
            interval_length=interval_length,
            interval_count=interval_count,
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
        interval_length: Optional[int],
        interval_count: Optional[int],
    ) -> list[int]:
        """
        Divide region based on interval size or count
        """
        region_length = end - start

        if interval_length:
            logger.info("Mapper - Use interval size to set up intervals")

            interval_boundaries = [
                i + start for i in range(0, region_length, interval_length)
            ]
            if interval_boundaries[-1] != region_length:
                interval_boundaries.append(end)

        elif interval_count:
            logger.info("Mapper - Use interval count to set up intervals")

            interval_size = round(region_length / interval_count)
            interval_boundaries = [
                start + (i * interval_size) for i in range(0, interval_count + 1)
            ]
            interval_boundaries[-1] = end

        return interval_boundaries

    def populate_read_counts(
        self, contig: str, beds: list[pd.DataFrame], bams: list[Loader]
    ) -> None:
        """
        Fill read counts in all BED DataFrames
        """
        logger.info(f"Mapper - Populate read counts in BED files")
        for bed, bam in zip(beds, bams):
            read_counts = [
                bam.bam.count(contig=contig, start=start, end=end)
                for start, end in zip(bed["start"], bed["end"])
            ]
            bed["read_count"] = read_counts

    # def populate_fractions(
    #     self, contig: str, beds: list[pd.DataFrame], bams: list[Loader]
    # ) -> None:
    #     """
    #     Fill read fractions in all BED DataFrames
    #     """
    #     logger.info(f"Mapper - Populate read fractions in BED files")
    #     for bed, bam in zip(beds, bams):
    #         bed["fraction"] = [
    #             bam.bam.count(contig=contig, start=row[1], end=row[2])
    #             / sum(bed["read_count"])
    #             for row in bed.itertuples(index=False)
    #         ]

    def write_beds(self, beds: list[pd.DataFrame], bed_paths: list[Path]) -> None:
        """
        Write output to BED file using filename specified
        """
        logger.info("Mapper - Write BED contents to file")
        for bed, path in zip(beds, bed_paths):
            bed.to_csv(path, sep="\t", index=False, header=False)
            super().check_file_exists(path=str(path))
