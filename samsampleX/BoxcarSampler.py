import logging
from collections import deque
from pathlib import Path
from typing import Deque, Dict, Generator, List

from samsampleX.FileHandler import FileHandler
from samsampleX.Loader import Loader

logger = logging.getLogger(__name__)


class BoxcarSampler(FileHandler):
    """Replay template start counts on a new BAM using streaming fetches."""

    def __init__(
        self,
        in_bam: str,
        chunk_starts: str,
        region: str,
        read_length: int,
        window_size: int,
        out_bam: str,
        tolerance: int,
    ) -> None:
        if read_length <= 0:
            raise ValueError("Read length must be greater than 0")
        if window_size < 0:
            raise ValueError("Window size must be non-negative")
        if tolerance < 0:
            raise ValueError("Read length variation must be non-negative")

        self.in_bam = in_bam
        self.chunk_starts = Path(chunk_starts)
        self.contig, self.region_start, self.region_end = self._parse_region(region)
        if self.region_start >= self.region_end:
            raise ValueError("--region start must be less than end")

        self.read_length = int(read_length)
        self.window_size = int(window_size)
        self.length_variation = int(tolerance)
        self.out_bam = out_bam
        self.template_start = max(0, self.region_start - self.read_length)

        self.read_cache: Dict[int, Deque] = {}
        self.start_order: Deque[int] = deque()

    def run_sampling(self) -> None:
        logger.info(
            "[BOXCAR-SAMPLE] Target=%s Template=%s Region=%s:%d-%d "
            "(template start=%d) ReadLength=%d Var=%d Window=%d",
            self.in_bam,
            self.chunk_starts,
            self.contig,
            self.region_start,
            self.region_end,
            self.template_start,
            self.read_length,
            self.length_variation,
            self.window_size,
        )
        self.check_file_exists(self.chunk_starts)

        target_loader = Loader(bam_path=self.in_bam)
        writer = Loader(bam_path=self.out_bam, template=target_loader.bam)

        fetch_start = max(self.template_start - self.window_size, 0)
        fetch_end = self.region_end + self.window_size
        read_iter = target_loader.fetch(
            contig=self.contig,
            start=fetch_start,
            end=fetch_end,
        )
        next_read = self._next_mapped_read(read_iter)

        kept_reads = 0
        requested_reads = 0
        with self.chunk_starts.open() as template_handle:
            for contig, position, count in self._template_records(template_handle):
                if contig != self.contig:
                    continue
                if position < self.template_start or position >= self.region_end:
                    continue

                requested_reads += count
                max_needed_pos = position + self.window_size
                next_read = self._fill_cache(
                    upto=max_needed_pos,
                    read_iter=read_iter,
                    next_read=next_read,
                )

                selected = self._collect_reads(position=position, needed=count)
                if len(selected) < count:
                    logger.warning(
                        "[BOXCAR-SAMPLE] Requested %d reads at %s:%d but only found %d",
                        count,
                        contig,
                        position,
                        len(selected),
                    )

                for read in selected:
                    writer.bam.write(read)
                    kept_reads += 1

                self._discard_coordinate(position)
                self._prune_old_coordinates(position)

        logger.info(
            "[BOXCAR-SAMPLE] Requested=%d Kept=%d (%.2f%% success)",
            requested_reads,
            kept_reads,
            (kept_reads / requested_reads * 100) if requested_reads else 0.0,
        )

        target_loader.close()
        writer.close()
        writer.sort_and_index()

    def _template_records(self, handle) -> Generator[tuple[str, int, int], None, None]:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 3:
                raise ValueError(
                    f"Template line '{line}' is not formatted as <contig> <pos> <count>"
                )
            contig, pos_str, count_str = parts
            yield contig, int(pos_str), int(count_str)

    def _fill_cache(
        self,
        upto: int,
        read_iter: Generator,
        next_read,
    ):
        while next_read is not None and next_read.reference_start is not None:
            start_pos = next_read.reference_start
            if start_pos > upto:
                break
            if start_pos not in self.read_cache:
                self.read_cache[start_pos] = deque()
                self.start_order.append(start_pos)
            self.read_cache[start_pos].append(next_read)
            next_read = self._next_mapped_read(read_iter)
        return next_read

    def _collect_reads(self, position: int, needed: int) -> List:
        selected = []
        for offset in self._offset_sequence():
            if needed == 0:
                break
            start_pos = position + offset
            bucket = self.read_cache.get(start_pos)
            if not bucket:
                continue

            while bucket and needed > 0:
                self._consume_invalid_reads(bucket)
                if not bucket:
                    break
                read = bucket.popleft()
                selected.append(read)
                needed -= 1

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

    def _discard_coordinate(self, position: int) -> None:
        if position in self.read_cache:
            self._remove_start(position)

    def _prune_old_coordinates(self, current_position: int) -> None:
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

    @staticmethod
    def _next_mapped_read(read_iter: Generator):
        for read in read_iter:
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


import logging
from collections import deque
from pathlib import Path
from typing import Deque, Dict, Generator, List

from samsampleX.FileHandler import FileHandler
from samsampleX.Loader import Loader

logger = logging.getLogger(__name__)


class BoxcarSampler(FileHandler):
    """Replay a template of read starts on a new BAM via streaming."""

    def __init__(
        self,
        in_bam: str,
        chunk_starts: str,
        region: str,
        read_length: int,
        window_size: int,
        out_bam: str,
        tolerance: int,
    ) -> None:
        if read_length <= 0:
            raise ValueError("BoxcarSampler requires --read-length > 0")
        if window_size < 0:
            raise ValueError("--window-size must be non-negative")
        if tolerance < 0:
            raise ValueError("--read-length-variation must be non-negative")

        self.in_bam = in_bam
        self.chunk_starts = Path(chunk_starts)
        self.contig, self.region_start, self.region_end = self._parse_region(region)
        if self.region_start >= self.region_end:
            raise ValueError("--region start must be less than end")

        self.read_length = int(read_length)
        self.window_size = int(window_size)
        self.length_variation = int(tolerance)
        self.out_bam = out_bam

        self.read_cache: Dict[int, Deque] = {}
        self.start_order: Deque[int] = deque()

    def run_sampling(self) -> None:
        logger.info(
            "[BOXCAR-SAMPLE] Target=%s Template=%s Region=%s:%d-%d "
            "ReadLength=%d Var=%d Window=%d",
            self.in_bam,
            self.chunk_starts,
            self.contig,
            self.region_start,
            self.region_end,
            self.read_length,
            self.length_variation,
            self.window_size,
        )
        self.check_file_exists(self.chunk_starts)

        target_loader = Loader(bam_path=self.in_bam)
        writer = Loader(bam_path=self.out_bam, template=target_loader.bam)

        fetch_start = max(self.region_start - self.window_size, 0)
        fetch_end = self.region_end + self.window_size
        read_iter = target_loader.fetch(
            contig=self.contig,
            start=fetch_start,
            end=fetch_end,
        )
        next_read = self._next_mapped_read(read_iter)

        kept_reads = 0
        requested_reads = 0
        with self.chunk_starts.open() as template_handle:
            for contig, position, count in self._template_records(template_handle):
                if contig != self.contig:
                    continue
                if position < self.region_start or position >= self.region_end:
                    continue

                requested_reads += count
                max_needed_pos = position + self.window_size
                next_read = self._fill_cache(
                    upto=max_needed_pos,
                    read_iter=read_iter,
                    next_read=next_read,
                )

                selected = self._collect_reads(position=position, needed=count)
                if len(selected) < count:
                    logger.warning(
                        "[BOXCAR-SAMPLE] Requested %d reads at %s:%d but only found %d",
                        count,
                        contig,
                        position,
                        len(selected),
                    )

                for read in selected:
                    writer.bam.write(read)
                    kept_reads += 1

                # Ensure reads starting at this coordinate are never reused
                self._discard_coordinate(position)
                self._prune_old_coordinates(position)

        logger.info(
            "[BOXCAR-SAMPLE] Requested=%d Kept=%d (%.2f%% success)",
            requested_reads,
            kept_reads,
            (kept_reads / requested_reads * 100) if requested_reads else 0.0,
        )

        target_loader.close()
        writer.close()
        writer.sort_and_index()

    def _template_records(self, handle) -> Generator[tuple[str, int, int], None, None]:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 3:
                raise ValueError(
                    f"Template line '{line}' is not formatted as <contig> <pos> <count>"
                )
            contig, pos_str, count_str = parts
            yield contig, int(pos_str), int(count_str)

    def _fill_cache(
        self,
        upto: int,
        read_iter: Generator,
        next_read,
    ):
        while next_read is not None and next_read.reference_start is not None:
            start_pos = next_read.reference_start
            if start_pos > upto:
                break
            if start_pos not in self.read_cache:
                self.read_cache[start_pos] = deque()
                self.start_order.append(start_pos)
            self.read_cache[start_pos].append(next_read)
            next_read = self._next_mapped_read(read_iter)
        return next_read

    def _collect_reads(self, position: int, needed: int) -> List:
        selected = []
        for offset in self._offset_sequence():
            if needed == 0:
                break
            start_pos = position + offset
            bucket = self.read_cache.get(start_pos)
            if not bucket:
                continue

            while bucket and needed > 0:
                self._consume_invalid_reads(bucket)
                if not bucket:
                    break
                read = bucket.popleft()
                selected.append(read)
                needed -= 1

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

    def _discard_coordinate(self, position: int) -> None:
        if position in self.read_cache:
            self._remove_start(position)

    def _prune_old_coordinates(self, current_position: int) -> None:
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

    @staticmethod
    def _next_mapped_read(read_iter: Generator):
        for read in read_iter:
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
