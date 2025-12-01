import logging

import numpy as np
from tqdm import tqdm

from samsampleX.Loader import Loader
from samsampleX.FileHandler import FileHandler

logger = logging.getLogger(__name__)


class Sampler(FileHandler):
    """
    Sample reads from a target BAM so that its depth matches a template BAM.
    """

    def __init__(
        self,
        source_path: str,
        template_path: str,
        region: str,
        seed: int,
        out_bam_path: str,
    ) -> None:
        """
        Initialize loaders and metadata needed for sampling.
        """
        logger.info("Sampler - Initialize Sampler")

        logger.info("[SAMPLER] - Set up inputs")
        self.target: Loader = Loader(bam_path=source_path, template=None)
        self.template: Loader = Loader(bam_path=template_path, template=None)

        region_contig, region_start, region_end = self._parse_region(region)
        self.contig: str = self.normalize_contig(contig=region_contig)

        if self.contig not in self.template.bam.references:
            raise ValueError(
                f"Contig {self.contig} not found in template BAM header "
                f"{self.template.bam_path}"
            )

        self.region_start: int = region_start
        self.region_end: int = region_end
        if self.region_start >= self.region_end:
            raise ValueError("Region start must be less than end coordinate")

        self.main_seed: int = seed

        logger.info("[SAMPLER] - Set up output")
        self.result: Loader = Loader(bam_path=out_bam_path, template=self.target.bam)

    def run_sampling(self) -> None:
        logger.info("[SAMPLER] - Begin read sampling")

        rng = np.random.default_rng(self.main_seed)
        weight_iter = self._weight_stream(
            chrom=self.contig, start=self.region_start, end=self.region_end
        )
        weight_buffer = _WeightBuffer(start=self.region_start, weight_iter=weight_iter)

        processed_reads = 0
        kept_reads = 0

        total_reads = self.target.bam.count(
            contig=self.contig, start=self.region_start, end=self.region_end
        )
        read_iter = self.target.fetch(
            contig=self.contig, start=self.region_start, end=self.region_end
        )
        for read in tqdm(
            read_iter,
            total=total_reads if total_reads >= 0 else None,
            desc="Sampling reads",
            unit=" reads",
        ):
            # if read.is_unmapped or read.is_secondary or read.is_supplementary:
            #     continue

            if read.reference_end is None:
                continue

            overlap_start = max(read.reference_start, self.region_start)
            overlap_end = min(read.reference_end, self.region_end)

            if overlap_end <= overlap_start:
                continue

            avg_weight = weight_buffer.mean(overlap_start, overlap_end)
            avg_weight = max(0.0, min(1.0, avg_weight))

            if rng.random() < avg_weight:
                self.result.bam.write(read)
                kept_reads += 1

            weight_buffer.discard_before(overlap_start)
            processed_reads += 1

        logger.info(
            "[SAMPLER] - Processed %d reads, retained %d", processed_reads, kept_reads
        )
        self.target.close()
        self.template.close()
        self.result.close()
        self.result.sort_and_index()

    def _weight_stream(self, chrom: str, start: int, end: int):
        """
        Yield per-base weights computed from template and target piles.
        """
        template_pile = self.template.bam.pileup(
            chrom,
            start,
            end,
            truncate=True,
            stepper="all",
            max_depth=0,
            multiple_iterators=True,
        )
        target_pile = self.target.bam.pileup(
            chrom,
            start,
            end,
            truncate=True,
            stepper="all",
            max_depth=0,
            multiple_iterators=True,
        )

        template_iter = iter(template_pile)
        target_iter = iter(target_pile)
        template_col = next(template_iter, None)
        target_col = next(target_iter, None)

        pos = start
        while pos < end:
            template_depth, template_col = self._depth_at_position(
                template_col, template_iter, pos
            )
            target_depth, target_col = self._depth_at_position(
                target_col, target_iter, pos
            )

            if target_depth <= 0:
                yield 0.0
            else:
                ratio = min(template_depth / target_depth, 1.0)
                yield float(ratio)

            pos += 1

    @staticmethod
    def _depth_at_position(current_col, iterator, pos: int):
        """
        Return depth for a specific genomic position from a pileup iterator.
        """
        while current_col is not None and current_col.reference_pos < pos:
            current_col = next(iterator, None)

        if current_col is not None and current_col.reference_pos == pos:
            depth = current_col.nsegments
            current_col = next(iterator, None)
        else:
            depth = 0

        return depth, current_col

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

    @staticmethod
    def _parse_region(region: str) -> tuple[str, int, int]:
        """
        Parse a region string formatted as contig:start-end.
        """
        if ":" not in region or "-" not in region:
            raise ValueError(f"Region '{region}' must be formatted as contig:start-end")

        contig_part, coords = region.split(":", 1)
        start_str, end_str = coords.replace(",", "").split("-", 1)

        try:
            start = int(start_str)
            end = int(end_str)
        except ValueError as exc:
            raise ValueError(f"Region coordinates must be integers: {region}") from exc

        if start < 0 or end < 0:
            raise ValueError("Region coordinates must be non-negative integers")

        return contig_part, start, end


class _WeightBuffer:
    """
    Maintain a sliding window of per-base weights produced by a generator.
    """

    def __init__(self, start: int, weight_iter):
        self._iter = weight_iter
        self._buffer_start = start
        self._weights: list[float] = []
        self._exhausted = False

    def _extend_to(self, upto: int) -> None:
        while self._buffer_start + len(self._weights) < upto:
            if self._exhausted:
                self._weights.append(0.0)
                continue

            try:
                weight = next(self._iter)
            except StopIteration:
                self._exhausted = True
                weight = 0.0

            self._weights.append(weight)

    def mean(self, start: int, end: int) -> float:
        if end <= start:
            return 0.0

        self._extend_to(end)

        rel_start = start - self._buffer_start
        rel_end = end - self._buffer_start

        if rel_start < 0:
            rel_start = 0
        rel_end = min(rel_end, len(self._weights))

        if rel_end <= rel_start:
            return 0.0

        total = 0.0
        for idx in range(rel_start, rel_end):
            total += self._weights[idx]

        length = rel_end - rel_start
        return total / length if length else 0.0

    def discard_before(self, pos: int) -> None:
        if pos <= self._buffer_start:
            return

        drop = min(pos - self._buffer_start, len(self._weights))
        if drop <= 0:
            return

        self._weights = self._weights[drop:]
        self._buffer_start += drop
