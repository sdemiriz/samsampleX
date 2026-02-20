"""Depth array computation from BAM files and region parsing utilities."""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass, field

import numpy as np
import pysam


@dataclass
class Region:
    """Genomic region parsed from a samtools-style string."""

    contig: str
    start: int = -1  # 0-based inclusive; -1 means whole contig
    end: int = -1  # 0-based exclusive; -1 means whole contig


@dataclass
class DepthArray:
    """Per-position depth values over a contiguous genomic region."""

    contig: str
    start: int
    end: int
    depths: np.ndarray = field(repr=False)

    @property
    def length(self) -> int:
        return self.end - self.start


def region_parse(region_str: str) -> Region:
    """Parse a samtools-style region string (e.g. chr1, chr1:1000-2000).

    Input coordinates are 1-based; the returned Region uses 0-based half-open
    coordinates.
    """
    m = re.match(r"^([^:]+)(?::(\d+)(?:-(\d+))?)?$", region_str)
    if not m:
        raise ValueError(f"Invalid region format: {region_str}")

    contig = m.group(1)
    start = int(m.group(2)) - 1 if m.group(2) else -1  # 1-based → 0-based
    end = int(m.group(3)) if m.group(3) else -1  # exclusive, stays as-is

    return Region(contig=contig, start=start, end=end)


def resolve_contig_name(header: pysam.AlignmentHeader, contig: str) -> str | None:
    """Resolve *contig* against a BAM header, trying chr-prefix variants.

    Returns the matching name from the header, or None.
    """
    references = set(header.references)

    if contig in references:
        return contig

    if contig.startswith("chr"):
        alt = contig[3:]
        if alt in references:
            print(f"Note: Using contig '{alt}' (matched from '{contig}')", file=sys.stderr)
            return alt
    else:
        alt = f"chr{contig}"
        if alt in references:
            print(f"Note: Using contig '{alt}' (matched from '{contig}')", file=sys.stderr)
            return alt

    return None


def get_contig_length(header: pysam.AlignmentHeader, contig: str) -> int | None:
    """Return contig length from a BAM header, handling chr-prefix mismatch."""
    resolved = resolve_contig_name(header, contig)
    if resolved is None:
        return None
    return header.get_reference_length(resolved)


def depth_from_bam(
    bam_path: str,
    contig: str,
    start: int,
    end: int,
) -> DepthArray:
    """Compute per-position depth for a region from an indexed BAM file.

    Uses a simplified read-iteration approach:
    for each read overlapping the region, increment depth at every covered
    position (CIGAR-unaware).
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        resolved = resolve_contig_name(bam.header, contig)
        if resolved is None:
            raise ValueError(f"Contig '{contig}' not found in BAM header")

        tid = bam.get_tid(resolved)
        contig_len = bam.get_reference_length(resolved)

        if start < 0:
            start = 0
        if end < 0 or end > contig_len:
            end = contig_len

        depths = np.zeros(end - start, dtype=np.int32)

        for read in bam.fetch(resolved, start, end):
            if read.flag & (0x4 | 0x100 | 0x200 | 0x400):
                # FUNMAP | FSECONDARY | FQCFAIL | FDUP
                continue

            r_start = read.reference_start
            r_end = read.reference_end
            if r_end is None:
                continue

            ov_start = max(r_start, start)
            ov_end = min(r_end, end)
            if ov_start >= ov_end:
                continue

            depths[ov_start - start : ov_end - start] += 1

    return DepthArray(contig=resolved, start=start, end=end, depths=depths)
