import heapq
import logging
from collections import deque
from pathlib import Path
from typing import Generator

from samsampleX.FileHandler import FileHandler
from samsampleX.Loader import Loader

logger = logging.getLogger(__name__)


class BoxcarMapper(FileHandler):
    """Generate a template of read start counts via boxcar deconvolution."""

    def __init__(
        self,
        template_bam: str,
        region: str,
        read_length: int,
        out_template: str,
    ) -> None:
        if read_length <= 0:
            raise ValueError("Read length must be greater than 0")

        self.template_bam = template_bam
        self.contig, self.region_start, self.region_end = self.parse_region(region)
        if self.region_start >= self.region_end:
            raise ValueError("Region end cannot be less than start")

        self.read_length = int(read_length)

        # Calculate start coordinate for proper initial sampling
        self.template_start = max(0, self.region_start - self.read_length)

        # Make output file and parent directory
        self.out_template = Path(out_template)
        self.out_template.parent.mkdir(parents=True, exist_ok=True)

    def run_mapping(self) -> None:
        logger.info("[MAPPER] - Begin mapping")

        loader = Loader(bam_path=self.template_bam)
        try:
            self.write_template(loader)
        finally:
            loader.close()

        logger.info(f"[MAPPER] Finished writing {self.out_template}")
        self.check_file_exists(self.out_template)

    def write_template(self, loader: Loader) -> None:
        logger.info(
            f"[MAPPER] Fetching {self.contig}:{self.template_start}-{self.region_end}"
        )

        # Get samtools read iterator, and first read from it
        read_iter = loader.fetch(
            contig=self.contig,
            start=self.template_start,
            end=self.region_end,
        )
        next_read = self.next_read(read_iter)

        # Track reads contributing depth to current position
        active_ends: list[int] = []
        start_history = deque(maxlen=self.read_length)

        # Get starting depth
        prev_depth, next_read = self.depth_at(
            position=self.template_start,
            active_ends=active_ends,
            read_iter=read_iter,
            next_read=next_read,
        )

        with self.out_template.open("w") as out_handle:

            # For all coordinates in region
            for pos in range(self.template_start, self.region_end):

                # Get depth at current position
                depth, next_read = self.depth_at(
                    position=pos,
                    active_ends=active_ends,
                    read_iter=read_iter,
                    next_read=next_read,
                )

                # Get depth one read length before current position
                prev_start = (
                    start_history[0] if len(start_history) == self.read_length else 0
                )

                # Count how many reads start at current position
                start_count = depth - prev_depth + prev_start

                # Write to output if any reads start at current position
                if start_count <= 0:
                    start_count = 0
                else:
                    out_handle.write(f"{self.contig}\t{pos}\t{start_count}\n")

                # Add current start count as most recent entry
                start_history.append(start_count)

                # Set up previous depth for next iteration
                prev_depth = depth

    def depth_at(
        self,
        position: int,
        active_ends: list[int],
        read_iter: Generator,
        next_read,
    ) -> tuple[int, object | None]:
        """Advance sweep-line to the requested position and return depth."""
        self.expire_reads(position, active_ends)
        next_read = self.advance_reads(
            upto=position,
            active_ends=active_ends,
            read_iter=read_iter,
            next_read=next_read,
        )
        return len(active_ends), next_read

    def advance_reads(
        self,
        upto: int,
        active_ends: list[int],
        read_iter: Generator,
        next_read,
    ):
        while next_read:
            # Fetch reads until end of requested position
            read_start = next_read.reference_start
            if read_start > upto:
                break

            # If read passes end of requested position, add to buffer
            read_end = read_start + self.read_length
            if read_end > upto:
                heapq.heappush(active_ends, read_end)

            # Fetch next read
            next_read = self.next_read(read_iter)
        return next_read

    @staticmethod
    def expire_reads(position: int, active_ends: list[int]) -> None:
        # Remove past reads
        while active_ends and active_ends[0] <= position:
            heapq.heappop(active_ends)

    @staticmethod
    def next_read(read_iter: Generator):
        """Return mapped reads, end with None"""
        # Skip unmapped and None reads
        for read in read_iter:
            if read.reference_start is None:
                continue
            if read.is_unmapped:
                continue

            return read
        return None

    @staticmethod
    def parse_region(region: str) -> tuple[str, int, int]:
        if ":" not in region or "-" not in region:
            raise ValueError("--region must be formatted as contig:start-end")

        # Split string into contig, start, end based on symbols
        contig_part, coords = region.split(":", 1)
        start_str, end_str = coords.replace(",", "").split("-", 1)

        start = int(start_str)
        end = int(end_str)

        if start < 0 or end < 0:
            raise ValueError("--region coordinates must be >= 0")

        if end <= start:
            raise ValueError("--region end must be greater than start")

        return contig_part, start, end
