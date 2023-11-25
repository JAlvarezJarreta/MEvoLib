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
"""Functions aimed to support different needs on the configuration and execution processes of MEvoLib."""

import os
import multiprocessing
from pathlib import Path
import tempfile


NUMCORES = multiprocessing.cpu_count()


def get_abspath(filename: str) -> str:
    """Returns the absolute path of `filename`.

    If `filename` is already an absolute path, it is returned without changes. Otherwise, the current
    working directory is used as the file's absolute path.

    Args:
        filename: File name to get the absolute path from.

    """
    if os.path.isdir(filename):
        raise ValueError("expected file, not dir")
    if not os.path.isabs(filename):
        return os.path.join(os.getcwd(), filename)
    else:
        return filename


def write_content(filepath, content) -> None:
    """TODO"""
    with Path(filepath).open("w") as fout:
        fout.write(content)


def get_tempfile_path() -> str:
    """Returns the path of a new temporary file name without creating it."""
    return os.path.join(
        tempfile.gettempdir(),
        tempfile.gettempprefix() + next(tempfile._get_candidate_names())
    )
