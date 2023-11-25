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

from contextlib import nullcontext as does_not_raise
import pytest

from mevolib import align


@pytest.mark.parametrize(
    "input,output,expectation",
    [
        (0, ["clustalo"], does_not_raise()),
        (3, ["clustalo", "mafft", "muscle"], does_not_raise()),
        (4, [], pytest.raises(ValueError)),
    ],
)
def test_get_tools(input: int, output: list, expectation: Exception):
    """Unit test of the method `mevolib.align.get_tools()`."""
    with expectation:
        assert align.get_tools(input) == output


def test_get_keywords():
    """Unit test of the method `mevolib.align.get_keywords()`."""
    assert align.get_keywords() == {}
