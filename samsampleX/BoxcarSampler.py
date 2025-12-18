import logging
from pathlib import Path
from collections import deque
from typing import Deque, Dict, Generator, List, Optional

import numpy as np

from samsampleX.Loader import Loader
from samsampleX.FileHandler import FileHandler

logger = logging.getLogger(__name__)


class BoxcarSampler(FileHandler):
    """Replay template start counts on a new BAM using streaming fetches."""

    def __init__(
        self,
        in_bam: str,
        template: str,
        region: str,
        read_length: int,
        window_size: int,
        out_bam: str,
        tolerance: int,
        seed: Optional[int] = None,
    ) -> None:
        if read_length <= 0:
            raise ValueError("Read length must be greater than 0")
        if window_size < 0:
            raise ValueError("Window size must be >0")
        if tolerance < 0:
            raise ValueError("Read length tolerance must be >0")

        # File IO
        self.in_bam = in_bam
        self.template = Path(template)
        self.out_bam = out_bam
        self.check_file_exists(self.in_bam)
        self.check_file_exists(self.template)
        self.check_file_exists(self.out_bam)

        # Parameters
        self.read_length = int(read_length)
        self.window_size = int(window_size)
        self.length_variation = int(tolerance)
        self.contig, self.region_start, self.region_end = self._parse_region(region)
        if self.region_start >= self.region_end:
            raise ValueError("--region start must be less than end")
        self.template_start = max(0, self.region_start - self.read_length)

        # Randomness and sampling
        self.rng = np.random.default_rng(seed)
        self.read_cache: Dict[int, Deque] = {}
        self.start_order: Deque[int] = deque()

    def run_sampling(self) -> None:
        logger.info("[BOXCAR-SAMPLE] - Begin sampling")

        # Set up BAM IO
        read_source = Loader(bam_path=self.in_bam)
        read_out = Loader(bam_path=self.out_bam, template=read_source.bam)

        # Fetch reads from defined region
        fetch_start = max(self.template_start - self.window_size, 0)
        fetch_end = self.region_end + self.window_size
        read_iterator = read_source.fetch(
            contig=self.contig,
            start=fetch_start,
            end=fetch_end,
        )
        # First read
        next_read = self.get_next_read(read_iterator)

        # Read tally
        kept_reads = 0
        requested_reads = 0

        with self.template.open() as template_handle:
            for contig, position, count in self.iterate_template(template_handle):
                # if contig != self.contig:
                #     continue
                # if position < self.template_start or position >= self.region_end:
                #     continue

                requested_reads += count

                # Add reads until end of requested window
                max_needed_pos = position + self.window_size
                next_read = self.fill_cache(
                    upto=max_needed_pos,
                    read_iterator=read_iterator,
                    next_read=next_read,
                )

                if count > 0:
                    selected = self.collect_reads(position=position, needed=count)
                    for read in selected:
                        read_out.bam.write(read)
                        kept_reads += 1

                self.clear_position(position)
                self.cache_cleanup(position)

        logger.info(
            "[BOXCAR-SAMPLE] - %.2f%% completion: fetched=%d/requested=%d",
            (kept_reads / requested_reads * 100) if requested_reads else 0.0,
            kept_reads,
            requested_reads,
        )

        read_source.close()
        read_out.close()
        read_out.sort_and_index()

    def iterate_template(self, handle) -> Generator[tuple[str, int, int], None, None]:
        for line in handle:
            line = line.strip()

            # Skip header and empty lines
            if not line or line.startswith("#"):
                continue

            # Extract three values per line
            parts = line.split()
            if len(parts) != 3:
                raise ValueError(
                    f"Template line '{line}' is not formatted as <contig> <pos> <count>"
                )
            contig, pos_str, count_str = parts
            yield contig, int(pos_str), int(count_str)

    def fill_cache(
        self,
        upto: int,
        read_iterator: Generator,
        next_read,
    ):
        while next_read:
            start_pos = next_read.reference_start
            if start_pos > upto:
                break

            # Add start_pos entry to cache if not present
            if start_pos not in self.read_cache:
                self.read_cache[start_pos] = deque()
                self.start_order.append(start_pos)

            # Else, add read to the same index
            self.read_cache[start_pos].append(next_read)
            next_read = self.get_next_read(read_iterator)

        return next_read

    def collect_reads(self, position: int, needed: int) -> List:
        selected = []
        for offset in self._offset_sequence():
            start_pos = position + offset
            bucket = self.read_cache.get(start_pos)
            if not bucket:
                continue

            # while bucket and needed > 0:
            #     self._consume_invalid_reads(bucket)
            #     if not bucket:
            #         break
            #     read = self._pop_random(bucket)
            #     selected.append(read)
            #     needed -= 1

            try:
                selected = self.rng.choice(bucket, size=needed, replace=False)
            except ValueError:
                logger.warning(
                    f"[BOXCAR-SAMPLE] - Position {position + offset} sampling deficit: {len(bucket)}/{needed} reads available"
                )
                selected = list(bucket)

            if not bucket:
                self._remove_start(start_pos)

        return selected

    def _offset_sequence(self) -> Generator[int, None, None]:
        yield 0
        for delta in range(1, self.window_size + 1):
            yield delta
            yield -delta

    def _consume_invalid_reads(self, bucket: Deque) -> None:
        while bucket:
            read = bucket[0]
            if read.query_length is None:
                bucket.popleft()
                continue
            delta = abs(read.query_length - self.read_length)
            if delta > self.length_variation:
                bucket.popleft()
                continue
            break

    def clear_position(self, position: int) -> None:
        if position in self.read_cache:
            self._remove_start(position)

    def cache_cleanup(self, current_position: int) -> None:
        min_allowed = current_position - self.window_size
        while self.start_order and self.start_order[0] < min_allowed:
            start_pos = self.start_order[0]
            self._remove_start(start_pos)

    def _remove_start(self, start_pos: int) -> None:
        bucket = self.read_cache.get(start_pos)
        if bucket:
            bucket.clear()
        self.read_cache.pop(start_pos, None)
        while self.start_order and self.start_order[0] == start_pos:
            self.start_order.popleft()

    def _pop_random(self, bucket: Deque):
        """Remove and return a random read from the deque using RNG."""
        if len(bucket) == 1:
            return bucket.popleft()
        idx = int(self.rng.integers(len(bucket)))
        bucket.rotate(-idx)
        read = bucket.popleft()
        bucket.rotate(idx)
        return read

    @staticmethod
    def get_next_read(read_iterator: Generator):
        for read in read_iterator:
            if read.is_unmapped or read.reference_start is None:
                continue
            if read.is_secondary or read.is_supplementary:
                continue
            return read
        return None

    @staticmethod
    def _parse_region(region: str):
        if ":" not in region or "-" not in region:
            raise ValueError("--region must be formatted as contig:start-end")
        contig_part, coords = region.split(":", 1)
        start_str, end_str = coords.replace(",", "").split("-", 1)
        start = int(start_str)
        end = int(end_str)
        if start < 0 or end < 0:
            raise ValueError("--region coordinates must be >= 0")
        return contig_part, start, end
