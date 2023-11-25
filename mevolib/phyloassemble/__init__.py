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
"""APIs to handle different phylogenetic assembly tools (supertrees, consensus trees, ...)."""

import os
import tempfile
import subprocess

from Bio import Phylo

from mevolib._utils import get_abspath
from mevolib.phyloassemble import Consense


_STREE_TOOL_TO_LIB = { }
_CONS_TOOL_TO_LIB = { 'consense': Consense }


def get_tools() -> dict:
    """Returns a dict of supertree and consensus tree tools included in the current version of MEvoLib."""
    return dict([
               ('supertree', list(_STREE_TOOL_TO_LIB.keys())),
               ('consensus', list(_CONS_TOOL_TO_LIB.keys()))
           ])


def get_keywords(tool: str) -> dict:
    """Returns a dictionary containing the keywords and their corresponding arguments.

    Args:
        tool: Name of the supertree or consensus tool.

    Raises:
        ValueError: If the tool introduced is not included in ``MEvoLib.PhyloAssemble``.
    """
    tool = tool.lower()
    tool_lib_keys = _STREE_TOOL_TO_LIB.keys() | _CONS_TOOL_TO_LIB.keys()
    if tool not in tool_lib_keys:
        raise ValueError(f'The tool "{tool}" is not included in MEvoLib.PhyloAssemble')
    keyword_dict = {}
    if tool in _STREE_TOOL_TO_LIB:
        tool_lib_dict = _STREE_TOOL_TO_LIB
    else:
        tool_lib_dict = _CONS_TOOL_TO_LIB
    for key, value in iter(tool_lib_dict[tool].KEYWORDS.items()) :
        keyword_dict[key] = ' '.join(value)
    return keyword_dict


def get_consensus_tree(binary: str, infile: str, infile_format: str, args: str = 'default',
                       outfile: str = None, outfile_format: str = 'newick' ) -> Bio.Phylo.BaseTree:
    """
    Calculate the consensus tree of the input trees file with the given
    arguments. The resultant consensus tree is returned as a Bio.Phylo.BaseTree
    object and saved in the output file (if provided). If 'infile' or 'outfile'
    contain a relative path, the current working directory will be used to get
    the absolute path. If the output file already exists, the old file will be
    overwritten without any warning.

    Args:
        binary: Name or path of the consensus tool.
        infile: Input phylogenetic trees file.
        infile_format: Input file format.
        args: Keyword or arguments to use in the call of the consensus tool,
            excluding infile and outfile arguments. By default, 'default'
            arguments are used.
            * For Consense, the second character will be used as separator of
            the different arguments. 
        outfile: Consensus tree output file.
        outfile_format: Output file format. By default, NEWICK format.

    Returns:
        Bio.Phylo.BaseTree
            Resultant consensus tree.

    Raises:
        ValueError: If the tool introduced isn't included in MEvoLib.
        IOError: If the input path or the input file provided doesn't exist.
        IOError: If the consensus tool did not generate a consensus tree.

    * The input file format must be supported by Bio.Phylo.
    * The output file format must be supported by Bio.Phylo.
    """
    # Get the variables associated with the given consensus tool
    bin_path, bin_name = os.path.split(binary)
    bin_name = bin_name.lower()
    if bin_name not in _CONS_TOOL_TO_LIB:
        raise ValueError(f'The consensus tool "{bin_name}" is not included in MEvoLib.PhyloAssemble')
    tool_lib = _CONS_TOOL_TO_LIB[bin_name]
    sprt_infile_formats = tool_lib.SPRT_INFILE_FORMATS
    gen_args = tool_lib.gen_args
    gen_stdin_content = tool_lib.gen_stdin_content
    get_results = tool_lib.get_results
    cleanup = tool_lib.cleanup
    # Get the command line to run in order to get the consensus tree
    infile_path = get_abspath(infile)
    # If the input file format is not supported by the consensus tool, convert
    # it to a temporary supported file
    if infile_format.lower() not in sprt_infile_formats:
        tmpfile = tempfile.NamedTemporaryFile()
        Phylo.convert(infile_path, infile_format, tmpfile.name, sprt_infile_formats[0])
        infile_path = tmpfile.name
    # Create full command line list
    command = [binary] + gen_args(args, infile_path, outfile)
    # Generate the standard input file content
    stdin_content = gen_stdin_content(args)
    # Create the input file with the given options
    with tempfile.NamedTemporaryFile(mode='w+') as stdin_file:
        stdin_file.write(stdin_content)
        stdin_file.seek(0)
        # Run the consensus process
        subprocess.run(command, stdin=stdin_file, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       universal_newlines=True, check=True)
        consensus_tree = get_results(command)
        cleanup(command)
        # Return the resultant consensus tree as a ``Bio.Phylo.BaseTree`` object
        return consensus_tree
