"""Tests for mapback.py: sequences.txt loading, offset table, header modification."""

import textwrap
from pathlib import Path

import pytest

from samsamplex.mapback import (
    ALT_CONTIG_MAPS_GRCH38,
    CONTIG_NAMES,
    GENE_MAPS_GRCH38,
    CHR6_LENGTHS,
    _resolve_prg_to_chr6_start,
    build_prg_offset_table,
    load_sequences_txt,
    modify_header,
    remap_read,
)


# ── load_sequences_txt ──────────────────────────────────────────────────────


class TestLoadSequencesTxt:
    def test_basic_parsing(self, tmp_path: Path):
        p = tmp_path / "sequences.txt"
        p.write_text(textwrap.dedent("""\
            SequenceID\tName\tFASTAID\tChr\tStart_1based\tStop_1based
            1\tapd\tPRG_1\t\t\t
            10\tF*01:01:01:08\tPRG_10\t\t\t
        """))
        result = load_sequences_txt(str(p))
        assert result == {"PRG_1": "apd", "PRG_10": "F*01:01:01:08"}

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError, match="--prg-seq"):
            load_sequences_txt("/nonexistent/sequences.txt")

    def test_empty_fastaid_skipped(self, tmp_path: Path):
        p = tmp_path / "sequences.txt"
        p.write_text(textwrap.dedent("""\
            SequenceID\tName\tFASTAID\tChr\tStart_1based\tStop_1based
            1\tapd\tPRG_1\t\t\t
            2\tempty_row\t\t\t\t
        """))
        result = load_sequences_txt(str(p))
        assert "PRG_1" in result
        assert len(result) == 1


# ── _resolve_prg_to_chr6_start ──────────────────────────────────────────────


class TestResolvePrgToChr6Start:
    def test_hla_gene_allele(self):
        offset = _resolve_prg_to_chr6_start(
            "A*01:01:01:01", GENE_MAPS_GRCH38, CONTIG_NAMES, ALT_CONTIG_MAPS_GRCH38,
        )
        assert offset == 29941260

    def test_alt_contig_name(self):
        offset = _resolve_prg_to_chr6_start(
            "apd", GENE_MAPS_GRCH38, CONTIG_NAMES, ALT_CONTIG_MAPS_GRCH38,
        )
        assert offset == 28734408

    def test_chr6_itself(self):
        offset = _resolve_prg_to_chr6_start(
            "chr6", GENE_MAPS_GRCH38, CONTIG_NAMES, ALT_CONTIG_MAPS_GRCH38,
        )
        assert offset == 1

    def test_unresolvable_returns_none(self):
        offset = _resolve_prg_to_chr6_start(
            "UNKNOWN_GENE", GENE_MAPS_GRCH38, CONTIG_NAMES, ALT_CONTIG_MAPS_GRCH38,
        )
        assert offset is None


# ── build_prg_offset_table ──────────────────────────────────────────────────


class TestBuildPrgOffsetTable:
    def test_known_entries(self):
        seq_map = {
            "PRG_1": "apd",
            "PRG_10": "F*01:01:01:08",
            "PRG_999": "UNKNOWN",
        }
        table = build_prg_offset_table(
            seq_map, GENE_MAPS_GRCH38, CONTIG_NAMES, ALT_CONTIG_MAPS_GRCH38,
        )
        assert table["PRG_1"] == 28734408   # apd → chr6_GL000250v2_alt start
        assert table["PRG_10"] == 29722775  # F gene start
        assert "PRG_999" not in table       # unresolvable excluded

    def test_empty_seq_map(self):
        table = build_prg_offset_table(
            {}, GENE_MAPS_GRCH38, CONTIG_NAMES, ALT_CONTIG_MAPS_GRCH38,
        )
        assert table == {}


# ── modify_header ───────────────────────────────────────────────────────────


class TestModifyHeader:
    def _make_header(self, sq_entries: list[dict]) -> "pysam.AlignmentHeader":
        import pysam
        return pysam.AlignmentHeader.from_dict({"SQ": sq_entries})

    def test_adds_chr6_when_missing(self):
        header = self._make_header([{"SN": "PRG_1", "LN": 5000}])
        new_header = modify_header(header, "GRCh38")
        contigs = {sq["SN"] for sq in new_header.to_dict()["SQ"]}
        assert "chr6" in contigs
        chr6_entry = [sq for sq in new_header.to_dict()["SQ"] if sq["SN"] == "chr6"][0]
        assert chr6_entry["LN"] == CHR6_LENGTHS["GRCh38"]

    def test_no_duplicate_chr6(self):
        header = self._make_header([
            {"SN": "chr6", "LN": 170805979},
            {"SN": "PRG_1", "LN": 5000},
        ])
        new_header = modify_header(header, "GRCh38")
        chr6_count = sum(1 for sq in new_header.to_dict()["SQ"] if sq["SN"] == "chr6")
        assert chr6_count == 1

    def test_grch37_length(self):
        header = self._make_header([{"SN": "PRG_1", "LN": 5000}])
        new_header = modify_header(header, "GRCh37")
        chr6_entry = [sq for sq in new_header.to_dict()["SQ"] if sq["SN"] == "chr6"][0]
        assert chr6_entry["LN"] == CHR6_LENGTHS["GRCh37"]


# ── remap_read ──────────────────────────────────────────────────────────────


class TestRemapRead:
    def _make_bam_and_read(self, ref_name: str, ref_start: int, mate_name: str = None):
        """Create a minimal in-memory BAM and read for testing."""
        import pysam

        sq = [
            {"SN": "PRG_1", "LN": 50000},
            {"SN": "PRG_10", "LN": 50000},
            {"SN": "chr6", "LN": 170805979},
        ]
        header = pysam.AlignmentHeader.from_dict({"SQ": sq})
        chr6_tid = header.get_tid("chr6")

        read = pysam.AlignedSegment(header)
        read.query_name = "test_read"
        read.query_sequence = "ACGT" * 25
        read.flag = 0
        read.reference_id = header.get_tid(ref_name)
        read.reference_start = ref_start
        read.cigar = [(0, 100)]
        read.mapping_quality = 60

        if mate_name:
            read.next_reference_id = header.get_tid(mate_name)
            read.next_reference_start = ref_start + 200
        else:
            read.next_reference_id = read.reference_id
            read.next_reference_start = ref_start + 200

        return header, chr6_tid, read

    def test_prg_read_remapped(self):
        offset_table = {"PRG_1": 28734408}
        header, chr6_tid, read = self._make_bam_and_read("PRG_1", 1000, "PRG_1")

        result = remap_read(read, chr6_tid, offset_table, header)
        assert result is True
        assert read.reference_id == chr6_tid
        assert read.reference_start == 28734408 + 1000

    def test_chr6_read_tid_updated(self):
        offset_table = {}
        header, chr6_tid, read = self._make_bam_and_read("chr6", 30000000, "chr6")

        result = remap_read(read, chr6_tid, offset_table, header)
        assert result is True
        assert read.reference_id == chr6_tid
        assert read.reference_start == 30000000  # not shifted

    def test_unresolvable_contig_returns_false(self):
        import pysam

        sq = [{"SN": "chrX", "LN": 155270560}, {"SN": "chr6", "LN": 170805979}]
        header = pysam.AlignmentHeader.from_dict({"SQ": sq})
        chr6_tid = header.get_tid("chr6")

        read = pysam.AlignedSegment(header)
        read.query_name = "test_read"
        read.query_sequence = "ACGT" * 25
        read.flag = 0
        read.reference_id = header.get_tid("chrX")
        read.reference_start = 1000
        read.cigar = [(0, 100)]

        result = remap_read(read, chr6_tid, {}, header)
        assert result is False

    def test_unresolvable_mate_set_unmapped(self):
        import pysam

        sq = [
            {"SN": "PRG_1", "LN": 50000},
            {"SN": "PRG_99", "LN": 50000},
            {"SN": "chr6", "LN": 170805979},
        ]
        header = pysam.AlignmentHeader.from_dict({"SQ": sq})
        chr6_tid = header.get_tid("chr6")

        offset_table = {"PRG_1": 28734408}

        read = pysam.AlignedSegment(header)
        read.query_name = "test_read"
        read.query_sequence = "ACGT" * 25
        read.flag = 0
        read.reference_id = header.get_tid("PRG_1")
        read.reference_start = 1000
        read.cigar = [(0, 100)]
        read.next_reference_id = header.get_tid("PRG_99")
        read.next_reference_start = 500

        result = remap_read(read, chr6_tid, offset_table, header)
        assert result is True
        assert read.reference_id == chr6_tid
        assert read.mate_is_unmapped
