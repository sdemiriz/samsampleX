"""Remap HLA*LA PRG-mapped reads back to canonical chr6 coordinates."""

from __future__ import annotations

import csv
import os
import sys
import tempfile
from pathlib import Path

import pysam

from .depth import region_parse, resolve_contig_name

# ── Short PRG haplotype names → GRCh38 alt contig names ─────────────────────

CONTIG_NAMES: dict[str, str] = {
    "apd": "chr6_GL000250v2_alt",
    "cox": "chr6_GL000251v2_alt",
    "dbb": "chr6_GL000252v2_alt",
    "mann": "chr6_GL000253v2_alt",
    "mcf": "chr6_GL000254v2_alt",
    "qbl": "chr6_GL000255v2_alt",
    "ssto": "chr6_GL000256v2_alt",
    "chr6": "chr6",
    "6": "chr6",
}

# ── HLA gene boundaries on chr6 (GENCODE v48, non-lncRNA) ───────────────────

GENE_MAPS_GRCH38: dict[str, tuple[str, int, int]] = {
    "A": ("chr6", 29941260, 29949572),
    "B": ("chr6", 31353872, 31367067),
    "C": ("chr6", 31268749, 31272130),
    "DMA": ("chr6", 32948613, 32969094),
    "DMB": ("chr6", 32934629, 32941028),
    "DOA": ("chr6", 33004182, 33009591),
    "DPA1": ("chr6", 33064569, 33080775),
    "DPB1": ("chr6", 33075936, 33089696),
    "DQA1": ("chr6", 32628179, 32647062),
    "DQB1": ("chr6", 32659467, 32668383),
    "DRA": ("chr6", 32439878, 32445046),
    "DRB1": ("chr6", 32577902, 32589848),
    "DRB3": ("chr6_GL000250v2_alt", 3824514, 3837642),
    "DRB4": ("chr6_GL000253v2_alt", 3840435, 3855431),
    "E": ("chr6", 30489509, 30494194),
    "F": ("chr6", 29722775, 29738528),
    "G": ("chr6", 29826967, 29831125),
    "H": ("chr6", 29887752, 29890482),
    "K": ("chr6", 29926459, 29929232),
    "L": ("chr6", 30259625, 30261703),
    "MICA": ("chr6", 31399784, 31415315),
    "MICB": ("chr6", 31494881, 31511124),
    "P": ("chr6", 29800415, 29802425),
    "TAP1": ("chr6", 32845209, 32853816),
    "TAP2": ("chr6", 32821833, 32838739),
    "V": ("chr6", 29792234, 29793136),
}

GENE_MAPS_GRCH37: dict[str, tuple[str, int, int]] = {
    "A": ("chr6", 29909037, 29917349),
    "B": ("chr6", 31321649, 31334844),
    "C": ("chr6", 31236526, 31239907),
    "DMA": ("chr6", 32916390, 32936871),
    "DMB": ("chr6", 32902406, 32908805),
    "DOA": ("chr6", 32971959, 32977368),
    "DPA1": ("chr6", 33032346, 33048552),
    "DPB1": ("chr6", 33043713, 33057473),
    "DQA1": ("chr6", 32595956, 32614839),
    "DQB1": ("chr6", 32627244, 32636160),
    "DRA": ("chr6", 32407655, 32412823),
    "DRB1": ("chr6", 32545679, 32557625),
    "DRB3": ("chr6_GL000250v2_alt", 3824514, 3837642),
    "DRB4": ("chr6_GL000253v2_alt", 3840435, 3855431),
    "E": ("chr6", 30457286, 30461971),
    "F": ("chr6", 29690552, 29706305),
    "G": ("chr6", 29794744, 29798902),
    "H": ("chr6", 29855529, 29858259),
    "K": ("chr6", 29894236, 29897009),
    "L": ("chr6", 30227402, 30229480),
    "MICA": ("chr6", 31367561, 31383092),
    "MICB": ("chr6", 31462658, 31478901),
    "P": ("chr6", 29768192, 29770202),
    "TAP1": ("chr6", 32812986, 32821593),
    "TAP2": ("chr6", 32789610, 32806516),
    "V": ("chr6", 29760011, 29760913),
}

# ── Alt contig → chr6 coordinate ranges (from UCSC Browser / liftover) ──────

ALT_CONTIG_MAPS_GRCH38: dict[str, tuple[str, int, int]] = {
    "chr6_GL000250v2_alt": ("chr6", 28734408, 33367716),
    "chr6_GL000251v2_alt": ("chr6", 28510120, 33383765),
    "chr6_GL000252v2_alt": ("chr6", 28734408, 33361299),
    "chr6_GL000253v2_alt": ("chr6", 28734408, 33258200),
    "chr6_GL000254v2_alt": ("chr6", 28734408, 33391865),
    "chr6_GL000255v2_alt": ("chr6", 28734408, 33411973),
    "chr6_GL000256v2_alt": ("chr6", 28691466, 33480577),
    "chr6": ("chr6", 1, 33480577),
}

ALT_CONTIG_MAPS_GRCH37: dict[str, tuple[str, int, int]] = {
    "chr6_GL000250v2_alt": ("chr6", 28702185, 33335493),
    "chr6_GL000251v2_alt": ("chr6", 28477897, 33351542),
    "chr6_GL000252v2_alt": ("chr6", 28702185, 33329076),
    "chr6_GL000253v2_alt": ("chr6", 28702185, 33225977),
    "chr6_GL000254v2_alt": ("chr6", 28702185, 33359642),
    "chr6_GL000255v2_alt": ("chr6", 28702185, 33379750),
    "chr6_GL000256v2_alt": ("chr6", 28659243, 33448354),
    "chr6": ("chr6", 60000, 33448354),
}

CHR6_LENGTHS: dict[str, int] = {
    "GRCh38": 170805979,
    "GRCh37": 171115067,
}


# ── sequences.txt loading ───────────────────────────────────────────────────


def load_sequences_txt(path: str) -> dict[str, str]:
    """Read HLA*LA sequences.txt and return a {FASTAID: Name} mapping.

    The file is tab-separated with columns: SequenceID, Name, FASTAID, Chr,
    Start_1based, Stop_1based.  Only the Name and FASTAID columns are used.
    """
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(
            f"sequences.txt not found at '{path}'. "
            "Specify the correct path with --prg-seq."
        )

    seq_map: dict[str, str] = {}
    with open(p, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            fasta_id = row["FASTAID"].strip()
            name = row["Name"].strip()
            if fasta_id:
                seq_map[fasta_id] = name

    return seq_map


def _resolve_prg_to_chr6_start(
    name: str,
    gene_maps: dict[str, tuple[str, int, int]],
    contig_names: dict[str, str],
    alt_contig_maps: dict[str, tuple[str, int, int]],
) -> int | None:
    """Resolve a sequences.txt Name to a chr6 start offset.

    Returns the chr6 start coordinate that PRG-local positions should be added
    to, or None if the name cannot be resolved.
    """
    gene_name = name.split("*")[0]

    if gene_name in gene_maps:
        _chrom, start, _end = gene_maps[gene_name]
        return start

    if name in contig_names:
        alt_contig = contig_names[name]
        if alt_contig in alt_contig_maps:
            _chrom, start, _end = alt_contig_maps[alt_contig]
            return start

    return None


def build_prg_offset_table(
    seq_map: dict[str, str],
    gene_maps: dict[str, tuple[str, int, int]],
    contig_names: dict[str, str],
    alt_contig_maps: dict[str, tuple[str, int, int]],
) -> dict[str, int]:
    """Pre-compute {PRG_N: chr6_start_offset} for every resolvable FASTAID.

    Called once at startup so that per-read coordinate conversion is a single
    dict lookup.
    """
    table: dict[str, int] = {}
    for fasta_id, name in seq_map.items():
        offset = _resolve_prg_to_chr6_start(name, gene_maps, contig_names, alt_contig_maps)
        if offset is not None:
            table[fasta_id] = offset
    return table


def modify_header(
    header: pysam.AlignmentHeader, genome_build: str,
) -> pysam.AlignmentHeader:
    """Return a new header with a chr6 SQ entry appended if not already present."""
    hd = header.to_dict()
    existing_contigs = {sq["SN"] for sq in hd.get("SQ", [])}

    if "chr6" not in existing_contigs:
        hd.setdefault("SQ", []).append(
            {"SN": "chr6", "LN": CHR6_LENGTHS[genome_build]}
        )

    return pysam.AlignmentHeader.from_dict(hd)


def remap_read(
    read: pysam.AlignedSegment,
    chr6_tid: int,
    offset_table: dict[str, int],
    header: pysam.AlignmentHeader,
) -> bool:
    """Remap a PRG-mapped read to chr6 coordinates in-place.

    Returns True if the read was successfully remapped, False if its primary
    contig could not be resolved (the read should be skipped).

    Mate contigs that cannot be resolved are set to unmapped rather than
    causing a hard failure.
    """
    ref_name = read.reference_name

    if ref_name in ("chr6", "6"):
        read.reference_id = chr6_tid
    elif ref_name is not None and ref_name.startswith("PRG"):
        if ref_name not in offset_table:
            return False
        read.reference_start = offset_table[ref_name] + read.reference_start
        read.reference_id = chr6_tid
    else:
        return False

    mate_name = read.next_reference_name
    if mate_name == ref_name or mate_name in ("chr6", "6"):
        read.next_reference_id = chr6_tid
    elif mate_name is not None and mate_name.startswith("PRG"):
        if mate_name in offset_table:
            read.next_reference_start = offset_table[mate_name] + read.next_reference_start
            read.next_reference_id = chr6_tid
        else:
            read.mate_is_unmapped = True
            read.next_reference_id = -1
            read.next_reference_start = 0
    else:
        read.mate_is_unmapped = True
        read.next_reference_id = -1
        read.next_reference_start = 0

    return True


def _get_prg_contigs(bam: pysam.AlignmentFile) -> list[str]:
    """Return contig names from the BAM header that are PRG, chr6, or bare 6."""
    return [
        c for c in bam.references
        if c.startswith("PRG") or c in ("chr6", "6")
    ]


# ── Main entry point ────────────────────────────────────────────────────────


def mapback_run(
    source_bam: str,
    region_str: str,
    out_bam: str,
    genome_build: str,
    prg_seq: str,
    no_sort: bool = False,
) -> int:
    """Remap PRG-mapped reads to chr6 and write a new BAM.  Returns 0 on success."""
    log = lambda msg: print(msg, file=sys.stderr)

    log(f"[mapback] Source BAM:    {source_bam}")
    log(f"[mapback] Region:        {region_str}")
    log(f"[mapback] Genome build:  {genome_build}")
    log(f"[mapback] sequences.txt: {prg_seq}")
    log(f"[mapback] Output BAM:    {out_bam}")

    # Load sequences.txt → {FASTAID: Name}
    try:
        seq_map = load_sequences_txt(prg_seq)
    except FileNotFoundError as exc:
        log(f"Error: {exc}")
        return 1
    log(f"[mapback] Loaded {len(seq_map)} entries from sequences.txt")

    # Select genome-build-specific maps
    if genome_build == "GRCh38":
        gene_maps = GENE_MAPS_GRCH38
        alt_contig_maps = ALT_CONTIG_MAPS_GRCH38
    else:
        gene_maps = GENE_MAPS_GRCH37
        alt_contig_maps = ALT_CONTIG_MAPS_GRCH37

    # Pre-compute PRG_N → chr6 offset (O(1) per-read lookup)
    offset_table = build_prg_offset_table(seq_map, gene_maps, CONTIG_NAMES, alt_contig_maps)
    log(f"[mapback] Built offset table with {len(offset_table)} resolvable PRG contigs")

    # Parse target region
    region = region_parse(region_str)

    with pysam.AlignmentFile(source_bam, "rb") as src:
        # Build output header with chr6
        out_header = modify_header(src.header, genome_build)
        chr6_tid = out_header.get_tid("chr6")

        # Resolve region bounds
        if region.start < 0:
            region.start = 0
        if region.end < 0:
            region.end = CHR6_LENGTHS[genome_build]

        log(f"[mapback] Target region: chr6:{region.start}-{region.end}")

        prg_contigs = _get_prg_contigs(src)
        log(f"[mapback] PRG/chr6 contigs in source: {len(prg_contigs)}")

        total = 0
        remapped = 0
        kept = 0
        skipped = 0

        with pysam.AlignmentFile(out_bam, "wb", header=out_header) as out:
            for contig in prg_contigs:
                for read in src.fetch(contig=contig):
                    if not read.is_mapped:
                        continue
                    total += 1

                    if not remap_read(read, chr6_tid, offset_table, out_header):
                        skipped += 1
                        continue
                    remapped += 1

                    r_start = read.reference_start
                    r_end = read.reference_end
                    if r_end is None:
                        continue

                    if max(r_start, region.start) < min(r_end, region.end):
                        out.write(read)
                        kept += 1

                    if total % 500_000 == 0:
                        log(f"[mapback]   Processed {total} reads...")

    log(f"[mapback] Processed {total} mapped reads")
    log(f"[mapback]   Remapped:  {remapped}")
    log(f"[mapback]   In region: {kept}")
    log(f"[mapback]   Skipped:   {skipped} (unresolvable contig)")

    if not no_sort:
        log("[mapback] Sorting output BAM...")
        tmp_sorted = out_bam + ".tmp.sorted.bam"
        try:
            pysam.sort("-o", tmp_sorted, out_bam)
            os.replace(tmp_sorted, out_bam)
            pysam.index(out_bam)
            log("[mapback] Sorting and indexing complete.")
        except Exception as exc:
            log(f"Warning: Failed to sort/index output BAM: {exc}")

    log(f"[mapback] Done. Output written to: {out_bam}")
    return 0
