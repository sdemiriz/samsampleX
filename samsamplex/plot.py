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
    (1_000_000, 10),   # <= 1 Mb: every 10th base
    (10_000_000, 100), # <= 10 Mb: every 100th base
]
DOWNSAMPLE_FALLBACK_STEP = 500  # > 10 Mb


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
    import matplotlib.ticker as mticker

    positions = np.arange(source.length)
    ds_pos, ds_src, ds_tpl, ds_out = _downsample(
        positions, source.depths, template.depths, output.depths, step=step,
    )

    genomic_pos = ds_pos + region_start
    region_end = region_start + source.length

    # Match Plotter.py styling
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

    fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")

    colors = {"Source": "#66c2a5", "Template": "#fc8d62", "Output": "#8da0cb"}
    for label, y_data, color in [
        ("Source", ds_src, colors["Source"]),
        ("Template", ds_tpl, colors["Template"]),
        ("Output", ds_out, colors["Output"]),
    ]:
        ax.plot(genomic_pos, y_data, label=label, color=color, linewidth=1)
        ax.fill_between(genomic_pos, y_data, alpha=0.3, color=color)

    # Remove top, right, left spines (keep bottom)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # Horizontal grid lines
    ax.set_axisbelow(True)
    ax.grid(True, which="major", axis="y", color="gray", alpha=0.5, linestyle="-", linewidth=1)

    ax.tick_params(axis="y", length=0)
    ax.set_xlabel(f"Position on {region_contig}")
    ax.set_ylabel("Depth")
    ax.set_title("Depth of coverage")

    # Custom x-axis formatter for large values (e.g. millions)
    def _custom_formatter(x, pos):
        return f"{x/1e6:.3f}m" if x >= 1e6 else f"{x:.0f}"

    ax.xaxis.set_major_formatter(mticker.FuncFormatter(_custom_formatter))
    ax.xaxis.set_major_locator(mticker.LinearLocator(numticks=5))

    # Figure-level legend at bottom center
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(handles), frameon=False)

    # Region coordinates at bottom left
    fig.text(0, -0.05, f"{region_contig}:{region_start+1:,}-{region_end:,}", fontsize=10, va="bottom", ha="left")

    if step > 1:
        ax.annotate(
            f"Displayed every {step}th position",
            xy=(0.01, 0.01),
            xycoords="axes fraction",
            fontsize=7,
            color="grey",
        )

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
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
