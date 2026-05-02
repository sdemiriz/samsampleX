"""Hash-based deterministic BAM sampling to match a template depth profile."""

from __future__ import annotations

import os
import sys
from typing import Sequence

import numpy as np
import pysam
import xxhash

from .bed import bed_combine_depths, bed_read_depths
from .depth import DepthArray, depth_from_bam, region_parse, resolve_contig_name
from .modes import DEPTH_MODES

UINT32_MAX = 0xFFFFFFFF

# Backwards-compatible aliases (same canonical set as map/sample --mode and --stat).
VALID_STAT_MODES = DEPTH_MODES
VALID_COMBINE_MODES = DEPTH_MODES


def _xxh32_fraction(qname: str, seed: int) -> float:
    """Hash a read name to a float in [0, 1)."""
    h = xxhash.xxh32(qname.encode(), seed=seed).intdigest()
    return h / UINT32_MAX


# ── Ratio helpers ────────────────────────────────────────────────────────────

# CIGAR op codes (BAM spec): only M(0), =(7), X(8) are "aligned" positions where the
# read actually has a base aligned to the reference. D(2), N(3) consume reference
# positions but the read has no base there, so we skip them when aggregating ratios.
# I(1), S(4), H(5), P(6) don't consume reference at all.
_ALIGNED_OPS = {0, 7, 8}        # M, =, X (count toward the stat)
_REF_OPS = {0, 2, 3, 7, 8}      # M, D, N, =, X (consume reference; advance ref_pos)


def _compute_ratios(template: DepthArray, source: DepthArray) -> np.ndarray:
    """ratio[i] = min(1.0, template[i] / source[i]), 0 where source is 0."""
    with np.errstate(divide="ignore", invalid="ignore"):
        ratios = np.where(
            source.depths == 0,
            0.0,
            np.minimum(1.0, template.depths.astype(np.float64) / source.depths.astype(np.float64)),
        )
    return ratios


def _aligned_ref_runs(read) -> list[tuple[int, int]]:
    """Walk a read's CIGAR and return [(ref_start, ref_end), ...] for M/=/X spans.

    D and N ops still advance the reference position (so subsequent aligned ops
    map to the right place) but they themselves are not emitted. I/S/H/P ops
    don't consume reference and are silently skipped.
    """
    runs: list[tuple[int, int]] = []
    cigar = read.cigartuples
    if not cigar:
        return runs
    pos = read.reference_start
    for op, length in cigar:
        if op in _ALIGNED_OPS:
            runs.append((pos, pos + length))
            pos += length
        elif op in _REF_OPS:
            pos += length
    return runs


def _collect_aligned_ratios(
    ratios: np.ndarray,
    region_start: int,
    region_end: int,
    runs: list[tuple[int, int]],
) -> np.ndarray | None:
    """Concatenate the slices of *ratios* corresponding to the clipped aligned runs.

    Returns None if no run overlaps the region (so callers can short-circuit to 0.0).
    """
    parts: list[np.ndarray] = []
    for rs, re in runs:
        cs = max(rs, region_start)
        ce = min(re, region_end)
        if cs < ce:
            parts.append(ratios[cs - region_start : ce - region_start])
    if not parts:
        return None
    return parts[0] if len(parts) == 1 else np.concatenate(parts)


def _get_mean_ratio(
    cumsum: np.ndarray,
    region_start: int,
    region_end: int,
    runs: list[tuple[int, int]],
) -> float:
    """Mean ratio over aligned reference positions, using cumsum per run for speed."""
    total = 0.0
    count = 0
    for rs, re in runs:
        cs = max(rs, region_start)
        ce = min(re, region_end)
        if cs >= ce:
            continue
        total += float(cumsum[ce - region_start] - cumsum[cs - region_start])
        count += ce - cs
    if count == 0:
        return 0.0
    return total / count


def _get_min_ratio(
    ratios: np.ndarray,
    region_start: int,
    region_end: int,
    runs: list[tuple[int, int]],
) -> float:
    arr = _collect_aligned_ratios(ratios, region_start, region_end, runs)
    if arr is None or arr.size == 0:
        return 0.0
    return float(arr.min())


def _get_max_ratio(
    ratios: np.ndarray,
    region_start: int,
    region_end: int,
    runs: list[tuple[int, int]],
) -> float:
    arr = _collect_aligned_ratios(ratios, region_start, region_end, runs)
    if arr is None or arr.size == 0:
        return 0.0
    return float(arr.max())


def _get_median_ratio(
    ratios: np.ndarray,
    region_start: int,
    region_end: int,
    runs: list[tuple[int, int]],
) -> float:
    arr = _collect_aligned_ratios(ratios, region_start, region_end, runs)
    if arr is None or arr.size == 0:
        return 0.0
    return float(np.median(arr))


def _get_random_ratio(
    ratios: np.ndarray,
    region_start: int,
    region_end: int,
    runs: list[tuple[int, int]],
    seed: int,
    read_start: int,
    read_end: int,
) -> float:
    """Pick one ratio from the read's aligned positions; index is deterministic in (span, seed).

    The hash key is the original (read_start, read_end) span so determinism per
    (span, seed) is preserved across the CIGAR-aware change; only the array
    being indexed differs.
    """
    arr = _collect_aligned_ratios(ratios, region_start, region_end, runs)
    if arr is None or arr.size == 0:
        return 0.0
    h = xxhash.xxh32(f"{read_start}:{read_end}".encode(), seed=seed).intdigest()
    return float(arr[h % arr.size])


# ── Main sampling routine ───────────────────────────────────────────────────


def sample_run(
    source_bam: str,
    template_beds: Sequence[str],
    region_str: str,
    out_bam: str = "out.bam",
    mode: str = "random",
    stat: str = "mean",
    seed: int = 42,
    # no_sort: bool = False,
    uniform_fraction: float | None = None,
    threads: int = 2,
) -> int:
    """Run the sample subcommand. Returns 0 on success."""
    log = lambda msg: print(msg, file=sys.stderr)

    log(f"[sample] Source BAM: {source_bam}")
    if uniform_fraction is not None:
        log(f"[sample] Uniform fraction: {uniform_fraction}")
    else:
        log(f"[sample] Template BEDs: {len(template_beds)} file(s)")
        for i, b in enumerate(template_beds):
            log(f"[sample]   {i + 1}: {b}")
        log(f"[sample] Stat: {stat}")
        log(f"[sample] Mode: {mode}")
    log(f"[sample] Region: {region_str}")
    log(f"[sample] Seed: {seed}")
    log(f"[sample] Output BAM: {out_bam}")

    region = region_parse(region_str)

    # Open source BAM to resolve region bounds
    with pysam.AlignmentFile(source_bam, "rb") as src:
        resolved = resolve_contig_name(src.header, region.contig)
        if resolved is None:
            log(f"Error: Contig '{region.contig}' not found in BAM")
            return 1
        region.contig = resolved

        if region.start < 0:
            region.start = 0
        if region.end < 0:
            region.end = src.get_reference_length(resolved)

    log(f"[sample] Parsed region: {region.contig}:{region.start}-{region.end}")

    if uniform_fraction is None:
        if mode not in DEPTH_MODES:
            log(
                f"Error: Unknown combine mode {mode!r} "
                f"(expected one of: {', '.join(DEPTH_MODES)})"
            )
            return 1
        if stat not in DEPTH_MODES:
            log(
                f"Error: Unknown stat mode {stat!r} "
                f"(expected one of: {', '.join(DEPTH_MODES)})"
            )
            return 1
        if not template_beds:
            log("Error: Template BED(s) required when not using --uniform")
            return 1

        # Load template depth(s)
        log("[sample] Loading template BED file(s)...")
        template_arrays = [
            bed_read_depths(bp, region.contig, region.start, region.end)
            for bp in template_beds
        ]

        if len(template_arrays) == 1:
            template_depth = template_arrays[0]
        else:
            log(f"[sample] Combining {len(template_arrays)} templates using '{mode}' mode...")
            template_depth = bed_combine_depths(template_arrays, mode=mode, seed=seed)

        # Compute source depth
        log("[sample] Computing source depth array...")
        source_depth = depth_from_bam(source_bam, region.contig, region.start, region.end)

        # Compute ratios
        log("[sample] Computing sampling ratios...")
        ratios = _compute_ratios(template_depth, source_depth)

        # Precompute cumulative sum (used by mean stat mode, always needed for default)
        cumsum = np.concatenate(([0.0], np.cumsum(ratios)))

        # Choose ratio-lookup function based on stat mode. Each lambda walks the
        # read's CIGAR (via _aligned_ref_runs) so only M/=/X reference positions
        # contribute; D/N/I/S/H/P are excluded.
        if stat == "mean":
            get_ratio = lambda r: _get_mean_ratio(
                cumsum, region.start, region.end, _aligned_ref_runs(r)
            )
        elif stat == "min":
            get_ratio = lambda r: _get_min_ratio(
                ratios, region.start, region.end, _aligned_ref_runs(r)
            )
        elif stat == "max":
            get_ratio = lambda r: _get_max_ratio(
                ratios, region.start, region.end, _aligned_ref_runs(r)
            )
        elif stat == "median":
            get_ratio = lambda r: _get_median_ratio(
                ratios, region.start, region.end, _aligned_ref_runs(r)
            )
        elif stat == "random":
            get_ratio = lambda r: _get_random_ratio(
                ratios, region.start, region.end, _aligned_ref_runs(r),
                seed, r.reference_start, r.reference_end,
            )
        else:
            log(f"Error: Unknown stat mode '{stat}'")
            return 1

    # Sampling loop
    log("[sample] Sampling reads...")
    total_reads = 0
    kept_reads = 0

    with pysam.AlignmentFile(source_bam, "rb") as src:
        with pysam.AlignmentFile(out_bam, "wb", header=src.header, threads=threads) as out:
            for read in src.fetch(region.contig, region.start, region.end):
                if read.is_unmapped:
                    continue

                total_reads += 1

                hash_frac = _xxh32_fraction(read.query_name, seed)
                if uniform_fraction is not None:
                    read_ratio = uniform_fraction
                else:
                    read_ratio = get_ratio(read)

                if hash_frac < read_ratio:
                    out.write(read)
                    kept_reads += 1

                if total_reads % 1_000_000 == 0:
                    pct = 100.0 * kept_reads / total_reads
                    log(f"[sample]   Processed {total_reads} reads, kept {kept_reads} ({pct:.1f}%)")

    pct = 100.0 * kept_reads / total_reads if total_reads else 0.0
    log(f"[sample] Processed {total_reads} reads, kept {kept_reads} ({pct:.1f}%)")

    # Sort and index
    # if not no_sort:
    #     log("[sample] Sorting output BAM...")
    #     tmp_sorted = out_bam + ".tmp.sorted.bam"
    #     try:
    #         pysam.sort("-o", tmp_sorted, out_bam)
    #         os.replace(tmp_sorted, out_bam)
    #         pysam.index(out_bam)
    #         log("[sample] Sorting and indexing complete.")
    #     except Exception as exc:
    #         log(f"Warning: Failed to sort/index output BAM: {exc}")

    log(f"[sample] Done. Output written to: {out_bam}")
    return 0
