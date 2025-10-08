"""
Tests for FileHandler class.
"""

import pytest
from pathlib import Path
from subsample_reads.FileHandler import FileHandler


class TestFileHandler:
    """Test cases for FileHandler class."""

    def test_check_file_exists_valid_file(self, tmp_path):
        """Test check_file_exists with a valid file."""

        test_file = tmp_path / "test_file.txt"
        test_file.write_text("test content")

        # Should not raise any exception
        FileHandler.check_file_exists(test_file)

    def test_check_file_exists_nonexistent_file(self):
        """Test check_file_exists with a nonexistent file."""
        with pytest.raises(FileNotFoundError):
            FileHandler.check_file_exists("DOES_NOT_EXIST.txt")

    def test_check_file_exists_directory(self, tmp_path):
        """Test check_file_exists with a directory instead of file."""
        # tmp_path itself is a directory, not a file
        with pytest.raises(FileNotFoundError) as exc_info:
            FileHandler.check_file_exists(tmp_path)

        # Verify the error message indicates it's not a file
        assert "not a file" in str(exc_info.value)
