"""BED file I/O and depth-array combination utilities."""

from __future__ import annotations

import random
from pathlib import Path
from typing import TextIO

import numpy as np

from .depth import DepthArray


def contig_names_match(name1: str, name2: str) -> bool:
    """Check whether two contig names match, ignoring a leading 'chr' prefix."""
    if name1 == name2:
        return True
    strip = lambda s: s[3:] if s.startswith("chr") else s
    return strip(name1) == strip(name2)


# ── Writing ──────────────────────────────────────────────────────────────────


def write_bed_entry(fp: TextIO, chrom: str, start: int, end: int, depth: int) -> None:
    """Write a single BED4 line."""
    fp.write(f"{chrom}\t{start}\t{end}\t{depth}\n")


def write_bed_output(fp: TextIO, arr: DepthArray, collapse: int = 0) -> None:
    """Write a DepthArray to BED format with optional collapsing.

    When *collapse* is 0 every position gets its own line.  When > 0,
    consecutive positions whose depth differs by <= *collapse* are merged
    into a single interval (using the depth of the first position in the
    interval.
    """
    if arr.length == 0:
        return

    if collapse == 0:
        for i in range(arr.length):
            pos = arr.start + i
            write_bed_entry(fp, arr.contig, pos, pos + 1, int(arr.depths[i]))
    else:
        interval_start = arr.start
        interval_depth = int(arr.depths[0])

        for i in range(1, arr.length):
            current_depth = int(arr.depths[i])
            if abs(current_depth - interval_depth) > collapse:
                write_bed_entry(fp, arr.contig, interval_start, arr.start + i, interval_depth)
                interval_start = arr.start + i
                interval_depth = current_depth

        write_bed_entry(fp, arr.contig, interval_start, arr.end, interval_depth)


# ── Reading ──────────────────────────────────────────────────────────────────


def bed_read_depths(
    bed_path: str,
    contig: str,
    region_start: int,
    region_end: int,
) -> DepthArray:
    """Read depth values from a BED4 file into a DepthArray for a region."""
    depths = np.zeros(region_end - region_start, dtype=np.int32)

    with open(bed_path) as fp:
        for line in fp:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 4:
                parts = line.split()
            if len(parts) < 4:
                continue

            chrom = parts[0]
            bed_start = int(parts[1])
            bed_end = int(parts[2])
            depth = int(parts[3])

            if not contig_names_match(chrom, contig):
                continue

            ov_start = max(bed_start, region_start)
            ov_end = min(bed_end, region_end)
            if ov_start >= ov_end:
                continue

            depths[ov_start - region_start : ov_end - region_start] = depth

    return DepthArray(contig=contig, start=region_start, end=region_end, depths=depths)


# ── Combining multiple templates ─────────────────────────────────────────────


def bed_combine_depths(
    arrays: list[DepthArray],
    mode: str = "min",
    seed: int = 42,
) -> DepthArray:
    """Combine multiple DepthArrays position-by-position.

    Supported *mode* values: ``"min"``, ``"max"``, ``"mean"``, ``"random"``.
    """
    if len(arrays) == 1:
        return DepthArray(
            contig=arrays[0].contig,
            start=arrays[0].start,
            end=arrays[0].end,
            depths=arrays[0].depths.copy(),
        )

    stacked = np.stack([a.depths for a in arrays], axis=0)  # (n_arrays, length)

    if mode == "min":
        combined = stacked.min(axis=0)
    elif mode == "max":
        combined = stacked.max(axis=0)
    elif mode == "mean":
        combined = (stacked.sum(axis=0) // len(arrays)).astype(np.int32)
    elif mode == "random":
        rng = random.Random(seed)
        mins = stacked.min(axis=0)
        maxs = stacked.max(axis=0)
        combined = np.array(
            [
                v if mn == mx else rng.randint(mn, mx)
                for mn, mx, v in zip(mins, maxs, mins)
            ],
            dtype=np.int32,
        )
    else:
        raise ValueError(f"Unknown combine mode: {mode}")

    ref = arrays[0]
    return DepthArray(contig=ref.contig, start=ref.start, end=ref.end, depths=combined)
