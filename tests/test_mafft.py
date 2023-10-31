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
"""`mevolib.align` module tests."""

import filecmp
import pytest

from mevolib import align

class MockStdOut:
    
    def __init__(self, filename):
        with filename.open("r") as fin:
            self.stdout = fin.read()


def test_get_alignment(tmp_path, mocker) -> None:
    """TODO"""
    expected_aln = pytest.flatfiles_dir / "f001.mafft_default.aln"
    mocked_subprocess = mocker.patch("subprocess.run")
    mocked_subprocess.return_value = MockStdOut(expected_aln)
    output_path = tmp_path / "output.aln"
    align.get_alignment(
        "mafft", str(pytest.flatfiles_dir / "f001.fasta"), "fasta", outfile=str(output_path)
    )
    assert filecmp.cmp(output_path, pytest.flatfiles_dir / "f001.mafft_default.aln", shallow=False)
