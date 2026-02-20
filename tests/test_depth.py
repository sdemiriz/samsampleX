"""Tests for depth.py: region_parse."""

import pytest

from samsamplex.depth import Region, region_parse


class TestRegionParse:
    def test_contig_only(self):
        r = region_parse("chr1")
        assert r.contig == "chr1"
        assert r.start == -1
        assert r.end == -1

    def test_contig_and_start(self):
        r = region_parse("chr1:1000")
        assert r.contig == "chr1"
        assert r.start == 999  # 1-based → 0-based
        assert r.end == -1

    def test_contig_start_end(self):
        r = region_parse("chr1:1000-2000")
        assert r.contig == "chr1"
        assert r.start == 999
        assert r.end == 2000  # end is exclusive, stays as-is

    def test_start_of_one(self):
        """Start=1 should map to 0-based 0."""
        r = region_parse("chrX:1-500")
        assert r.start == 0
        assert r.end == 500

    def test_numeric_contig(self):
        r = region_parse("21:100-200")
        assert r.contig == "21"
        assert r.start == 99
        assert r.end == 200

    def test_large_coordinates(self):
        r = region_parse("chr6:29941260-29949572")
        assert r.contig == "chr6"
        assert r.start == 29941259
        assert r.end == 29949572

    def test_invalid_empty(self):
        with pytest.raises(ValueError, match="Invalid region format"):
            region_parse("")

    def test_invalid_colon_only(self):
        with pytest.raises(ValueError, match="Invalid region format"):
            region_parse(":")

    def test_invalid_no_contig(self):
        with pytest.raises(ValueError, match="Invalid region format"):
            region_parse(":1000-2000")
