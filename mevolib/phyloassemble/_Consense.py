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
"""MEvoLib's variables library functions to ease the usage of Consense from PHYLIP."""

import os
import tempfile
from typing import List

from Bio import Phylo

from mevolib._utils import get_abspath


SPRT_INFILE_FORMATS = ['newick']
KEYWORDS = {'default': ['R', '2', 'Y']} # majority rule consensus


def gen_args(args: str, infile_path: str, outfile: str) -> List:
    """Returns the list of arguments generated from `args` and the infile path requested to call Consense.

    Args:
        args: Keyword or arguments to use in the call of Consense, excluding infile and outfile arguments.
        infile_path: Input alignment file path.
        outfile: Consensus tree output file.

    """
    if outfile:
        outfile_path = get_abspath(outfile)
    else:
        # Output files will be saved in temporary files to retrieve the consensus tree
        outfile_path = os.path.join(
            tempfile.gettempdir(),
            tempfile.gettempprefix() + next(tempfile._get_candidate_names())
        )
    # Create full command line list
    argument_list = [infile_path, outfile_path]
    return argument_list

def gen_stdin_content(args: str) -> str:
    """Returns the standard input content generated from `args`.

    Args:
        args: Keyword or arguments to use in the call of Consense, excluding infile and outfile arguments.
            If `args` is not a keyword, the second character will be used as separator of the different
            arguments.

    """
    if args in KEYWORDS:
        options = KEYWORDS[args]
    else:
        options = [opt for opt in args.split(args[1])]
    # If the output file already exists, overwrite it
    options.append('R')
    stdin_content = '\n'.join(options) + '\n'
    return stdin_content


def get_results(command: List) -> Bio.Phylo.BaseTree:
    """Returns consensus tree from the files generated during the execution of `command`.

    Args:
        command: Consense's command line executed.

    Raises:
        IOError: If the consensus tool did not generate a consensus tree (indicated by user's options).

    """
    outfile_path = command[2]
    try:
        consensus_tree = Phylo.read(outfile_path, 'newick')
    except IOError as e:
        cleanup(command)
        raise IOError('The given arguments did not generate a consensus tree file')
    else:
        return consensus_tree


def cleanup(command: List) -> None:
    """Removes the temporary file created (if any) in `gen_args()` function.

    Args:
        command: Consense's command line executed.

    """
    logfile_path = get_abspath(command[2])
    if os.path.dirname(logfile_path) == tempfile.gettempdir():
        os.remove(logfile_path)
