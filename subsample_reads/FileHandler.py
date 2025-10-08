from pathlib import Path
import logging
from typing import Union

logger = logging.getLogger(__name__)


class FileHandler:
    """Handle file presence or absence"""

    @staticmethod
    def check_file_exists(path: Union[str, Path]) -> None:
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"File not found: {p.absolute()}")
        if not p.is_file():
            raise FileNotFoundError(f"Path exists but is not a file: {p.absolute()}")
