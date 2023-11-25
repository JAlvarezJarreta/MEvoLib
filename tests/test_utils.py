# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""`mevolib._utils` module tests."""

import pathlib

import pytest

from mevolib import _utils


# def test_get_abspath_1():
#     """Unit test of `mevolib._utils.get_abspath()` for an absolute path."""
#     assert _utils.get_abspath("/home/user/workbench")  == "/home/user/workbench"


# def test_get_abspath_2():
#     """Unit test of `mevolib._utils.get_abspath()` for a relative path."""
#     assert _utils.get_abspath("utils.py") == "utils.py"


def test_create_file(tmp_dir: pathlib.Path):
    d = tmp_dir / "sub"
    d.mkdir()
    p = d / "hello.txt"
    p.write_text("content", encoding="utf-8")
    assert p.read_text(encoding="utf-8") == "content"
    assert len(list(tmp_dir.iterdir())) == 1
