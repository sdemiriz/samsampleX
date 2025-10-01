#!/usr/bin/env python

import os
import sys
import time
import argparse
from datetime import datetime as dt
from logging import INFO, basicConfig, getLogger

from subsample_reads.Loader import Loader
from subsample_reads.Mapper import Mapper
from subsample_reads.Plotter import Plotter
from subsample_reads.Comparator import Comparator

# Ensure log directory exists
os.makedirs("log", exist_ok=True)


def configure_logging(mode):
    """Configure logging with mode-specific suffix."""
    timestamp = dt.now().strftime("%d%m_%H%M%S")
    log_filename = f"log/{timestamp}_{mode}.txt"

    basicConfig(
        filename=log_filename,
        level=INFO,
        format="{asctime} [{levelname}]: {message}",
        style="{",
        datefmt="%H:%M:%S",
    )
    logger = getLogger("subsample_reads")

    # Log the command line invocation for reproducibility
    command_line = sys.argv.copy()
    command_line[0] = "subsample-reads"
    command_line = " ".join(command_line)
    logger.info(f"Command: {command_line}")
    logger.info("Begin log")

    return logger


def mapper_mode(args):
    """Produce a distribution of read counts from given BAM file."""
    logger = configure_logging("map")

    try:
        start_time = time.time()

        # Validate that if bed_files is provided, it matches the number of BAM files
        if args.bed and len(args.bed) != len(args.in_bam):
            raise ValueError(
                f"Number of BED files ({len(args.bed)}) must match number of BAM files ({len(args.in_bam)})"
            )

        Mapper(
            bam_paths=args.in_bam,
            contig=args.contig,
            start=args.start,
            end=args.end,
            interval_length=args.interval_length,
            interval_count=args.interval_count,
            bed_dir=args.bed_dir,
            bed=args.bed,
        )
        logger.info("End log\n")

    except Exception as e:
        logger.error(f"Mapper - Exception encountered. Details:\n{e}", exc_info=True)
        raise e

    finally:
        end_time = time.time()
        logger.info(f"Mapper - Time taken: {end_time - start_time} seconds")


def sample_mode(args):
    """Sample provided BAM file(s) based on regions in BED file."""
    logger = configure_logging("sample")

    try:
        start_time = time.time()

        in_bam = Loader(bam_path=args.in_bam)

        # Determine if we're in HLA-LA mode based on --prg flag
        if args.prg:
            # HLA-LA mode: use PRG-specific sampling with specified genome build
            in_bam.run_sampling(
                hlala_mode=True,
                bed_dir=args.bed_dir,
                bed_file=args.bed,
                main_seed=args.seed,
                out_bam=args.out_bam,
                genome_build=args.prg,
            )
        else:
            # Regular mode: use standard sampling
            in_bam.run_sampling(
                hlala_mode=False,
                bed_dir=args.bed_dir,
                bed_file=args.bed,
                main_seed=args.seed,
                out_bam=args.out_bam,
            )

        in_bam.close()
        logger.info("End log\n")
    except Exception as e:
        logger.error(f"Sample - Exception encountered. Details:\n{e}", exc_info=True)
        raise e

    finally:
        end_time = time.time()
        logger.info(f"Sample - Time taken: {end_time - start_time} seconds")


def compare_mode(args):
    """Compare two BAM files to compare read data before and after downsampling."""
    logger = configure_logging("compare")

    try:
        start_time = time.time()

        Comparator(
            bam_left_path=args.bam_left, bam_right_path=args.bam_right, out=args.out
        )
        logger.info("End log\n")
    except Exception as e:
        logger.error(f"Compare - Exception encountered. Details:\n{e}", exc_info=True)
        raise e

    finally:
        end_time = time.time()
        logger.info(f"Compare - Time taken: {end_time - start_time} seconds")


def plotter_mode(args):
    """Plot a distribution of depth of coverage and read count in each interval."""
    logger = configure_logging("plot")

    try:
        start_time = time.time()

        plotter = Plotter(
            in_bam=args.in_bam,
            map_bam=args.map_bam,
            out_bam=args.out_bam,
            bed=args.bed,
            out_plt=args.out_plt,
            no_det=args.no_det,
        )

        plotter.plot()
        logger.info("End log\n")
    except Exception as e:
        logger.error(f"Plotter - Exception encountered. Details:\n{e}", exc_info=True)
        raise e

    finally:
        end_time = time.time()
        logger.info(f"Plotter - Time taken: {end_time - start_time} seconds")


def main():
    """Main CLI entry point for subsample-reads toolkit."""
    parser = argparse.ArgumentParser(
        prog="python -m subsample-reads",
        description="Toolkit for mapping, sampling, and plotting BAM files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    subparsers = parser.add_subparsers(
        required=True,
        title="Functions",
        description="Use -h flag with any subcommand to learn usage.",
    )

    # Mapping
    mapper = subparsers.add_parser(
        "map",
        help="Generate read depth distribution(s) from supplied BAM file(s) and write to BED file(s).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    mapper.add_argument(
        "--in-bam",
        nargs="+",
        required=True,
        help="BAM files to sample read counts from.",
    )
    mapper.add_argument(
        "--contig",
        required=True,
        help="A valid contig name present in all provided BAM file(s).",
    )
    mapper.add_argument(
        "--start", required=True, help="Start coordinate of the region to map."
    )
    mapper.add_argument(
        "--end", required=True, help="End coordinate of the region to map."
    )
    bed_output = mapper.add_mutually_exclusive_group()
    bed_output.add_argument(
        "-d", "--bed-dir", default="bed/", help="Directory for output BED files."
    )
    bed_output.add_argument(
        "--bed",
        nargs="+",
        help="BED file names (must match number of input BAM files).",
    )
    intervals = mapper.add_mutually_exclusive_group(required=True)
    intervals.add_argument(
        "--interval-length", default=None, help="Length of intervals to generate."
    )
    intervals.add_argument(
        "--interval-count", default=None, help="Number of intervals to generate."
    )
    mapper.set_defaults(func=mapper_mode)

    # Sampling
    sample = subparsers.add_parser(
        "sample",
        help="Apply generated read depth distribution(s) from selected BED file(s) to supplied BAM file. Use --prg flag for HLA-LA PRG back-mapping.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sample.add_argument(
        "--in-bam", required=True, help="Target BAM file to downsample."
    )
    bed_selection = sample.add_mutually_exclusive_group()
    bed_selection.add_argument(
        "--bed-dir", default="bed/", help="Directory to fetch BED files from."
    )
    bed_selection.add_argument(
        "--bed", default=None, help="Specify one BED file to downsample from."
    )
    sample.add_argument(
        "--seed", default=42, type=int, help="Seed for random downsampling."
    )
    sample.add_argument("--out-bam", default="out.bam", help="Output BAM file.")
    sample.add_argument(
        "--prg",
        choices=["GRCh37", "GRCh38"],
        help="Enable HLA-LA PRG mode for downsampling with specified genome build (GRCh37 or GRCh38).",
    )
    sample.set_defaults(func=sample_mode)

    # Compare
    compare = subparsers.add_parser(
        "compare",
        help="Compare two BAM files to see how many reads overlap.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    compare.add_argument(
        "--bam-left", required=True, help="Path to first BAM for the comparison."
    )
    compare.add_argument(
        "--bam-right", required=True, help="Path to second BAM for the comparison."
    )
    compare.add_argument(
        "--out", default="out.tsv", help="Output file for cross-mapping info."
    )
    compare.set_defaults(func=compare_mode)

    # Plotting
    plotter = subparsers.add_parser(
        "plot",
        help="Plot BAM file(s) read depth together with supplied BED file(s).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    plotter.add_argument("--in-bam", help="Path to input/original BAM file.")
    plotter.add_argument("--map-bam", help="Path to mapping BAM file.")
    plotter.add_argument("--out-bam", help="Path to downsampled BAM file.")

    plotter.add_argument("--bed", default=None, help="Specific BED file to plot.")
    plotter.add_argument(
        "--out-plt", default="out.png", help="Path for the output plot."
    )
    plotter.add_argument(
        "--no-det",
        action="store_true",
        help="Set det parameter to False for the Plotter class.",
    )
    plotter.set_defaults(func=plotter_mode)

    args = parser.parse_args()

    # Execute the selected mode function
    args.func(args)


if __name__ == "__main__":
    main()
