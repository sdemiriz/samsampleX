"""Depth-of-coverage plotting with adaptive downsampling for large regions."""

from __future__ import annotations

import sys
from typing import TextIO

import numpy as np

from .depth import DepthArray

# Downsampling thresholds: (max_positions, step_size)
# Regions smaller than the first threshold are plotted at full resolution.
# Larger regions use progressively coarser sampling.
DOWNSAMPLE_TIERS = [
    (10_000, 1),        # <= 10 kb: every base
    (100_000, 10),      # <= 100 kb: every 10th base
    (1_000_000, 100),   # <= 1 Mb: every 100th base
    (10_000_000, 1000), # <= 10 Mb: every 1000th base
]
DOWNSAMPLE_FALLBACK_STEP = 5000  # > 10 Mb


def _pick_step(n_positions: int) -> int:
    """Choose a sampling step size based on region length."""
    for threshold, step in DOWNSAMPLE_TIERS:
        if n_positions <= threshold:
            return step
    return DOWNSAMPLE_FALLBACK_STEP


def _downsample(positions: np.ndarray, *arrays: np.ndarray, step: int) -> tuple[np.ndarray, ...]:
    """Subsample arrays at *step* intervals, keeping the last point."""
    idx = np.arange(0, len(positions), step)
    if idx[-1] != len(positions) - 1:
        idx = np.append(idx, len(positions) - 1)
    return tuple(a[idx] for a in (positions, *arrays))


# ── TSV output ───────────────────────────────────────────────────────────────


def write_tsv(
    fp: TextIO,
    region_start: int,
    source: DepthArray,
    template: DepthArray,
    output: DepthArray,
    step: int,
) -> None:
    """Write depth data to TSV, downsampled at *step* intervals."""
    positions = np.arange(source.length)
    ds_pos, ds_src, ds_tpl, ds_out = _downsample(
        positions, source.depths, template.depths, output.depths, step=step,
    )

    fp.write("position\tsource_depth\ttemplate_depth\toutput_depth\n")
    for i in range(len(ds_pos)):
        pos_1based = region_start + 1 + int(ds_pos[i])
        fp.write(f"{pos_1based}\t{ds_src[i]}\t{ds_tpl[i]}\t{ds_out[i]}\n")


# ── PNG output ───────────────────────────────────────────────────────────────


def write_png(
    out_path: str,
    region_start: int,
    region_contig: str,
    source: DepthArray,
    template: DepthArray,
    output: DepthArray,
    step: int,
) -> None:
    """Render a depth comparison line plot to PNG using matplotlib."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    positions = np.arange(source.length)
    ds_pos, ds_src, ds_tpl, ds_out = _downsample(
        positions, source.depths, template.depths, output.depths, step=step,
    )

    genomic_pos = ds_pos + region_start

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(genomic_pos, ds_src, linewidth=0.6, alpha=0.8, color="#3366cc", label="Source")
    ax.plot(genomic_pos, ds_tpl, linewidth=0.6, alpha=0.8, color="#33aa55", label="Template")
    ax.plot(genomic_pos, ds_out, linewidth=0.6, alpha=0.8, color="#cc3333", label="Output")

    ax.set_xlabel(f"Position on {region_contig}")
    ax.set_ylabel("Depth")
    ax.set_title("Depth of Coverage Comparison")
    ax.legend(loc="upper right", fontsize="small")

    if step > 1:
        ax.annotate(
            f"Displayed every {step}th position",
            xy=(0.01, 0.01),
            xycoords="axes fraction",
            fontsize=7,
            color="grey",
        )

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ── Main entry point ─────────────────────────────────────────────────────────


def plot_run(
    source_bam: str,
    out_bam: str,
    region_str: str,
    template_bam: str | None = None,
    template_bed: str | None = None,
    out_png: str | None = None,
    out_tsv: str | None = None,
) -> int:
    """Run the plot subcommand. Returns 0 on success."""
    from .bed import bed_read_depths
    from .depth import depth_from_bam, region_parse, resolve_contig_name

    import pysam

    log = lambda msg: print(msg, file=sys.stderr)

    # Resolve region
    region = region_parse(region_str)

    with pysam.AlignmentFile(source_bam, "rb") as bam:
        resolved = resolve_contig_name(bam.header, region.contig)
        if resolved is None:
            log(f"Error: Contig '{region.contig}' not found in BAM")
            return 1
        region.contig = resolved
        if region.start < 0:
            region.start = 0
        if region.end < 0:
            region.end = bam.get_reference_length(resolved)

    log(f"[plot] Region: {region.contig}:{region.start + 1}-{region.end}")

    # Load depth arrays
    log(f"[plot] Loading source depths from: {source_bam}")
    source_depth = depth_from_bam(source_bam, region.contig, region.start, region.end)

    if template_bam:
        log(f"[plot] Loading template depths from BAM: {template_bam}")
        template_depth = depth_from_bam(template_bam, region.contig, region.start, region.end)
    else:
        log(f"[plot] Loading template depths from BED: {template_bed}")
        template_depth = bed_read_depths(template_bed, region.contig, region.start, region.end)

    log(f"[plot] Loading output depths from: {out_bam}")
    output_depth = depth_from_bam(out_bam, region.contig, region.start, region.end)

    if source_depth.length != template_depth.length or source_depth.length != output_depth.length:
        log("Error: Depth array length mismatch")
        return 1

    n = source_depth.length
    step = _pick_step(n)
    displayed = len(range(0, n, step))
    log(f"[plot] {n} positions, displaying ~{displayed} points (step={step})")

    if out_tsv:
        log(f"[plot] Writing TSV to: {out_tsv}")
        with open(out_tsv, "w") as fp:
            write_tsv(fp, region.start, source_depth, template_depth, output_depth, step)
    else:
        log(f"[plot] Writing PNG to: {out_png}")
        write_png(out_png, region.start, region.contig,
                  source_depth, template_depth, output_depth, step)

    log("[plot] Done.")
    return 0
