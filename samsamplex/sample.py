"""Hash-based deterministic BAM sampling to match a template depth profile."""

from __future__ import annotations

import os
import sys
import tempfile
from typing import Sequence

import numpy as np
import pysam
import xxhash

from .bed import bed_combine_depths, bed_read_depths
from .depth import DepthArray, depth_from_bam, region_parse, resolve_contig_name
from .metrics import metrics_calculate, metrics_print

UINT32_MAX = 0xFFFFFFFF

VALID_STAT_MODES = ("mean", "min", "max", "median")
VALID_COMBINE_MODES = ("min", "max", "mean", "random")


def _xxh32_fraction(qname: str, seed: int) -> float:
    """Hash a read name to a float in [0, 1)."""
    h = xxhash.xxh32(qname.encode(), seed=seed).intdigest()
    return h / UINT32_MAX


# ── Ratio helpers ────────────────────────────────────────────────────────────


def _compute_ratios(template: DepthArray, source: DepthArray) -> np.ndarray:
    """ratio[i] = min(1.0, template[i] / source[i]), 0 where source is 0."""
    with np.errstate(divide="ignore", invalid="ignore"):
        ratios = np.where(
            source.depths == 0,
            0.0,
            np.minimum(1.0, template.depths.astype(np.float64) / source.depths.astype(np.float64)),
        )
    return ratios


def _get_mean_ratio(
    cumsum: np.ndarray, region_start: int, region_end: int, read_start: int, read_end: int,
) -> float:
    """O(1) mean ratio over a read span using the precomputed cumulative sum."""
    cs = max(read_start, region_start)
    ce = min(read_end, region_end)
    if cs >= ce:
        return 0.0
    i1 = cs - region_start
    i2 = ce - region_start
    return float(cumsum[i2] - cumsum[i1]) / (ce - cs)


def _get_min_ratio(
    ratios: np.ndarray, region_start: int, region_end: int, read_start: int, read_end: int,
) -> float:
    cs = max(read_start, region_start)
    ce = min(read_end, region_end)
    if cs >= ce:
        return 0.0
    return float(ratios[cs - region_start : ce - region_start].min())


def _get_max_ratio(
    ratios: np.ndarray, region_start: int, region_end: int, read_start: int, read_end: int,
) -> float:
    cs = max(read_start, region_start)
    ce = min(read_end, region_end)
    if cs >= ce:
        return 0.0
    return float(ratios[cs - region_start : ce - region_start].max())


def _get_median_ratio(
    ratios: np.ndarray, region_start: int, region_end: int, read_start: int, read_end: int,
) -> float:
    cs = max(read_start, region_start)
    ce = min(read_end, region_end)
    if cs >= ce:
        return 0.0
    return float(np.median(ratios[cs - region_start : ce - region_start]))


# ── Main sampling routine ───────────────────────────────────────────────────


def sample_run(
    source_bam: str,
    template_beds: Sequence[str],
    region_str: str,
    out_bam: str = "out.bam",
    mode: str = "min",
    stat: str = "mean",
    seed: int = 42,
    no_sort: bool = False,
    no_metrics: bool = False,
) -> int:
    """Run the sample subcommand. Returns 0 on success."""
    log = lambda msg: print(msg, file=sys.stderr)

    log(f"[sample] Source BAM: {source_bam}")
    log(f"[sample] Template BEDs: {len(template_beds)} file(s)")
    for i, b in enumerate(template_beds):
        log(f"[sample]   {i + 1}: {b}")
    log(f"[sample] Region: {region_str}")
    log(f"[sample] Stat: {stat}")
    log(f"[sample] Mode: {mode}")
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

    # Choose ratio-lookup function based on stat mode
    if stat == "mean":
        get_ratio = lambda rs, re: _get_mean_ratio(cumsum, region.start, region.end, rs, re)
    elif stat == "min":
        get_ratio = lambda rs, re: _get_min_ratio(ratios, region.start, region.end, rs, re)
    elif stat == "max":
        get_ratio = lambda rs, re: _get_max_ratio(ratios, region.start, region.end, rs, re)
    elif stat == "median":
        get_ratio = lambda rs, re: _get_median_ratio(ratios, region.start, region.end, rs, re)
    else:
        log(f"Error: Unknown stat mode '{stat}'")
        return 1

    # Sampling loop
    log("[sample] Sampling reads...")
    total_reads = 0
    kept_reads = 0

    with pysam.AlignmentFile(source_bam, "rb") as src:
        with pysam.AlignmentFile(out_bam, "wb", header=src.header) as out:
            for read in src.fetch(region.contig, region.start, region.end):
                if read.is_unmapped:
                    continue

                total_reads += 1

                hash_frac = _xxh32_fraction(read.query_name, seed)
                ratio = get_ratio(read.reference_start, read.reference_end)

                if hash_frac < ratio:
                    out.write(read)
                    kept_reads += 1

                if total_reads % 1_000_000 == 0:
                    pct = 100.0 * kept_reads / total_reads
                    log(f"[sample]   Processed {total_reads} reads, kept {kept_reads} ({pct:.1f}%)")

    pct = 100.0 * kept_reads / total_reads if total_reads else 0.0
    log(f"[sample] Processed {total_reads} reads, kept {kept_reads} ({pct:.1f}%)")

    # Sort and index
    if not no_sort:
        log("[sample] Sorting output BAM...")
        tmp_sorted = out_bam + ".tmp.sorted.bam"
        try:
            pysam.sort("-o", tmp_sorted, out_bam)
            os.replace(tmp_sorted, out_bam)
            pysam.index(out_bam)
            log("[sample] Sorting and indexing complete.")
        except Exception as exc:
            log(f"Warning: Failed to sort/index output BAM: {exc}")

    # Metrics
    if not no_metrics:
        log("[sample] Computing metrics...")
        try:
            output_depth = depth_from_bam(out_bam, region.contig, region.start, region.end)
            result = metrics_calculate(template_depth, output_depth)
            metrics_print(result, label_a="Template", label_b="Output")
        except Exception as exc:
            log(f"Warning: Could not compute metrics: {exc}")

    log(f"[sample] Done. Output written to: {out_bam}")
    return 0
