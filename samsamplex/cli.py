"""Command-line interface for samsamplex."""

from __future__ import annotations

import argparse
import sys

from . import __version__


def _add_map_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "map",
        help="Extract depth of coverage from BAM to BED template",
    )
    p.add_argument("--template-bam", required=True, help="Input BAM file")
    p.add_argument("--region", required=True, help="Target region (samtools-style)")
    p.add_argument("--out-bed", default="out.bed", help="Output BED file [default: out.bed]")
    p.add_argument(
        "--collapse",
        type=int,
        default=0,
        help="Merge consecutive positions with depth diff <= INT [default: 0]",
    )


def _add_sample_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "sample",
        help="Sample reads from BAM to match template depth distribution",
    )
    p.add_argument("--source-bam", required=True, help="Input BAM file to sample from")
    p.add_argument(
        "--template-bed",
        required=True,
        nargs="+",
        help="Template BED file(s) with depth values",
    )
    p.add_argument("--region", required=True, help="Target region (samtools-style)")
    p.add_argument("--out-bam", default="out.bam", help="Output BAM file [default: out.bam]")
    p.add_argument(
        "--mode",
        default="random",
        choices=("min", "max", "mean", "random"),
        help="How to combine multiple templates [default: random]",
    )
    p.add_argument(
        "--stat",
        default="mean",
        choices=("mean", "min", "max", "median"),
        help="Statistic for summarising ratio over read span [default: mean]",
    )
    p.add_argument("--seed", type=int, default=42, help="Random seed [default: 42]")
    p.add_argument("--no-sort", action="store_true", help="Skip sorting/indexing output BAM")
    p.add_argument("--no-metrics", action="store_true", help="Skip metrics calculation")


def _add_plot_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "plot",
        help="Compare depth of coverage and output as PNG plot or TSV data",
    )
    p.add_argument("--source-bam", required=True, help="Source BAM file")
    p.add_argument("--out-bam", required=True, help="Output BAM file (from sampling)")
    p.add_argument("--region", required=True, help="Target region (samtools-style)")

    tpl = p.add_mutually_exclusive_group(required=True)
    tpl.add_argument("--template-bam", help="Template BAM file")
    tpl.add_argument("--template-bed", help="Template BED file")

    out = p.add_mutually_exclusive_group(required=True)
    out.add_argument("--out-png", help="Output PNG plot file")
    out.add_argument("--out-tsv", help="Output TSV data file")


def _add_stats_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "stats",
        help="Compare depth distributions between two BAM files",
    )
    p.add_argument("--bam-a", required=True, help="First BAM file (reference)")
    p.add_argument("--bam-b", required=True, help="Second BAM file (comparison)")
    p.add_argument("--region", required=True, help="Region to compare (samtools-style)")


# ── Subcommand handlers ─────────────────────────────────────────────────────


def _run_map(args: argparse.Namespace) -> int:
    from .bed import write_bed_output
    from .depth import depth_from_bam, get_contig_length, region_parse

    log = lambda msg: print(msg, file=sys.stderr)

    region = region_parse(args.region)

    log(f"[map] Template BAM: {args.template_bam}")
    log(f"[map] Region: {args.region}")
    log(f"[map] Collapse: {args.collapse}")
    log(f"[map] Output BED: {args.out_bed}")

    import pysam

    with pysam.AlignmentFile(args.template_bam, "rb") as bam:
        from .depth import resolve_contig_name

        resolved = resolve_contig_name(bam.header, region.contig)
        if resolved is None:
            log(f"Error: Contig '{region.contig}' not found in BAM")
            return 1
        region.contig = resolved
        contig_len = bam.get_reference_length(resolved)

    if region.start < 0:
        region.start = 0
    if region.end < 0:
        region.end = contig_len

    log(f"[map] Parsed region: {region.contig}:{region.start}-{region.end}")
    log("[map] Computing depth array (this may take a while)...")

    depth = depth_from_bam(args.template_bam, region.contig, region.start, region.end)
    log(f"[map] Computed depth for {depth.length} positions")

    suffix = " (collapsed)" if args.collapse > 0 else ""
    log(f"[map] Writing BED file{suffix}...")

    with open(args.out_bed, "w") as fp:
        write_bed_output(fp, depth, collapse=args.collapse)

    log(f"[map] Done. Output written to: {args.out_bed}")
    return 0


def _run_sample(args: argparse.Namespace) -> int:
    from .sample import sample_run

    return sample_run(
        source_bam=args.source_bam,
        template_beds=args.template_bed,
        region_str=args.region,
        out_bam=args.out_bam,
        mode=args.mode,
        stat=args.stat,
        seed=args.seed,
        no_sort=args.no_sort,
        no_metrics=args.no_metrics,
    )


def _run_plot(args: argparse.Namespace) -> int:
    from .plot import plot_run

    return plot_run(
        source_bam=args.source_bam,
        out_bam=args.out_bam,
        region_str=args.region,
        template_bam=args.template_bam,
        template_bed=args.template_bed,
        out_png=args.out_png,
        out_tsv=args.out_tsv,
    )


def _run_stats(args: argparse.Namespace) -> int:
    from .depth import depth_from_bam, region_parse
    from .metrics import metrics_calculate, metrics_print

    log = lambda msg: print(msg, file=sys.stderr)

    region = region_parse(args.region)

    import pysam

    with pysam.AlignmentFile(args.bam_a, "rb") as bam:
        from .depth import resolve_contig_name

        resolved = resolve_contig_name(bam.header, region.contig)
        if resolved is None:
            log(f"Error: Contig '{region.contig}' not found in BAM A")
            return 1
        region.contig = resolved
        if region.start < 0:
            region.start = 0
        if region.end < 0:
            region.end = bam.get_reference_length(resolved)

    log(f"Computing depth for BAM A: {args.bam_a}")
    log(f"Computing depth for BAM B: {args.bam_b}")
    log(f"Region: {region.contig}:{region.start + 1}-{region.end}")

    depth_a = depth_from_bam(args.bam_a, region.contig, region.start, region.end)
    depth_b = depth_from_bam(args.bam_b, region.contig, region.start, region.end)

    result = metrics_calculate(depth_a, depth_b)

    log("")
    log("========== Comparison Metrics ==========")
    log(f"BAM A (reference):       {args.bam_a}")
    log(f"BAM B (comparison):      {args.bam_b}")
    log(f"Region:                  {region.contig}:{region.start + 1}-{region.end}")
    log(f"Positions compared:      {depth_a.length}")
    log("-" * 41)
    metrics_print(result, label_a="BAM A", label_b="BAM B")

    return 0


# ── Entry point ──────────────────────────────────────────────────────────────


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="samsamplex",
        description=f"samsamplex v{__version__} — Depth-aware BAM file sampling",
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"samsamplex v{__version__}"
    )

    subparsers = parser.add_subparsers(dest="command")
    _add_map_parser(subparsers)
    _add_sample_parser(subparsers)
    _add_plot_parser(subparsers)
    _add_stats_parser(subparsers)

    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    dispatch = {
        "map": _run_map,
        "sample": _run_sample,
        "plot": _run_plot,
        "stats": _run_stats,
    }

    sys.exit(dispatch[args.command](args))


if __name__ == "__main__":
    main()
