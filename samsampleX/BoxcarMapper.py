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
        tolerance: int,
        out_chunks: str,
    ) -> None:
        if read_length <= 0:
            raise ValueError("Read length must be greater than 0")

        self.template_bam = template_bam
        self.contig, self.region_start, self.region_end = self.parse_region(region)
        if self.region_start >= self.region_end:
            raise ValueError("Region end cannot be less than start")

        self.read_length = int(read_length)
        self.tolerance = int(tolerance)
        self.template_start = max(0, self.region_start - self.read_length)
        self.fetch_start = max(0, self.template_start - self.read_length)
        self.out_chunks = Path(out_chunks)
        self.out_chunks.parent.mkdir(parents=True, exist_ok=True)

    def run_mapping(self) -> None:
        logger.info(
            "[BOXCAR-MAP] Template=%s Region=%s:%d-%d "
            "(template span=%d-%d, fetch=%d-%d) ReadLength=%d Var=%d",
            self.template_bam,
            self.contig,
            self.region_start,
            self.region_end,
            self.template_start,
            self.region_end,
            self.fetch_start,
            self.region_end,
            self.read_length,
            self.tolerance,
        )

        loader = Loader(bam_path=self.template_bam)
        try:
            self._write_template(loader)
        finally:
            loader.close()

        logger.info("[BOXCAR-MAP] Finished writing %s", self.out_chunks)
        self.check_file_exists(self.out_chunks)

    def _write_template(self, loader: Loader) -> None:
        fetch_start = self.fetch_start
        logger.info(
            "[BOXCAR-MAP] Fetching %s:%d-%d",
            self.contig,
            fetch_start,
            self.region_end,
        )

        read_iter = loader.fetch(
            contig=self.contig,
            start=fetch_start,
            end=self.region_end,
        )
        next_read = self.next_read(read_iter)
        active_ends: list[int] = []
        start_history = deque(maxlen=self.read_length)

        if self.template_start > 0:
            prev_depth, next_read = self._depth_at(
                position=self.template_start - 1,
                active_ends=active_ends,
                read_iter=read_iter,
                next_read=next_read,
            )
        else:
            prev_depth = 0

        total_positions = 0
        total_template_reads = 0
        with self.out_chunks.open("w") as out_handle:
            for pos in range(self.template_start, self.region_end):
                depth, next_read = self._depth_at(
                    position=pos,
                    active_ends=active_ends,
                    read_iter=read_iter,
                    next_read=next_read,
                )

                prev_start = (
                    start_history[0] if len(start_history) == self.read_length else 0
                )
                start_count = depth - prev_depth + prev_start
                if start_count < 0:
                    start_count = 0

                if start_count > 0:
                    out_handle.write(f"{self.contig}\t{pos}\t{start_count}\n")
                    total_template_reads += start_count

                start_history.append(start_count)
                prev_depth = depth
                total_positions += 1

        logger.info(
            "[BOXCAR-MAP] Processed %d positions (%d-%d); template starts=%d",
            total_positions,
            self.template_start,
            self.region_end,
            total_template_reads,
        )

    def _depth_at(
        self,
        position: int,
        active_ends: list[int],
        read_iter: Generator,
        next_read,
    ) -> tuple[int, object | None]:
        """Advance sweep-line to the requested position and return depth."""
        self._expire_reads(position, active_ends)
        next_read = self._advance_reads(
            upto=position,
            active_ends=active_ends,
            read_iter=read_iter,
            next_read=next_read,
        )
        return len(active_ends), next_read

    def _advance_reads(
        self,
        upto: int,
        active_ends: list[int],
        read_iter: Generator,
        next_read,
    ):
        while next_read is not None and next_read.reference_start is not None:
            read_start = next_read.reference_start
            if read_start > upto:
                break

            read_end = read_start + self.read_length
            if read_end > upto:
                heapq.heappush(active_ends, read_end)
            next_read = self.next_read(read_iter)
        return next_read

    @staticmethod
    def _expire_reads(position: int, active_ends: list[int]) -> None:
        while active_ends and active_ends[0] <= position:
            heapq.heappop(active_ends)

    @staticmethod
    def next_read(read_iter: Generator):
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
        contig_part, coords = region.split(":", 1)
        start_str, end_str = coords.replace(",", "").split("-", 1)
        start = int(start_str)
        end = int(end_str)
        if start < 0 or end < 0:
            raise ValueError("--region coordinates must be >= 0")
        return contig_part, start, end
