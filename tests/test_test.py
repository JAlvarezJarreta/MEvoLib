from contextlib import nullcontext as does_not_raise
from pathlib import Path
from typing import ContextManager

import pytest

from mevolib import _utils


@pytest.mark.parametrize(
    "filename, output, expected",
    [
        ("/home/user/workbench", "/home/user/workbench", does_not_raise()),
        ("workbench", str(Path("workbench").absolute()), does_not_raise()),
        ("docs", str(Path("docs").absolute()), pytest.raises(ValueError)),
    ],
)
def test_get_abspath(filename: str, output: str, expected: ContextManager) -> None:
    """TODO"""
    with expected:
        assert _utils.get_abspath(filename) == output


@pytest.mark.parametrize(
    "filename, content",
    [
        ("myfile.txt", "content")
    ],
)
def test_write_content(tmp_path: Path, filename: str, content: str) -> None:
    """TODO"""
    filepath = tmp_path / filename
    _utils.write_content(filepath, content)
    with filepath.open("r") as fin:
        text = fin.read()
    assert text == content
