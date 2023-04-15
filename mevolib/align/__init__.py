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
"""Functions aimed to provide an easy interface to handle different alignment tools."""

import argparse
from io import StringIO
import os
import tempfile
from typing import Optional
import subprocess

from Bio import SeqIO, AlignIO

from mevolib.align import _ClustalOmega, _Mafft, _Muscle
from mevolib._utils import get_abspath

import Bio.Align

_TOOL_TO_LIB = {
    'clustalo': _ClustalOmega,
    'mafft': _Mafft,
    'muscle': _Muscle,
}


def get_tools() -> list:
    """Returns a list of alignment software tools included in the current version of MEvoLib."""
    return list(_TOOL_TO_LIB.keys())


def get_keywords(tool: str) -> dict:
    """Returns a dictionary containing the keywords and their corresponding arguments.

    Args:
        tool: Name of the alignment tool.

    Raises:
        ValueError: If the tool introduced is not included in `mevolib.align`.

    """
    tool = tool.lower()
    if tool not in _TOOL_TO_LIB:
        raise ValueError(f'The alignment tool "{tool}" is not included in MEvoLib.Align')
    keyword_dict = {}
    for key, value in _TOOL_TO_LIB[tool].KEYWORDS.items():
        keyword_dict[key] = ' '.join(value)
    return keyword_dict


def get_alignment(binary: str, infile: str, infile_format: str, args: str = 'default',
                  outfile: Optional[str] = None, outfile_format: str = 'fasta',
                  **kwargs: dict) -> Bio.Align.MultipleSeqAlignment:
    """Returns the result of aligning the sequences using the given alignment tool and arguments.

    The resultant alignment is returned as a Bio.Align.MultipleSeqAlign object and saved in the output
    file (if provided). If `infile` or `outfile` contain a relative path, the current working directory
    will be used to get the absolute path. If the output file already exists, the old file will be
    overwritten without any warning.

    The alignment tool might not be included in MEvoLib, but it can still be used passing in `**kwargs`
    the keys "informats" and "incmd" with the list of supported infile formats and the infile argument
    (empty string if it does not need any), respectively.

    Args:
        binary: Name or path of the alignment tool.
        infile: Unaligned input sequence file.
        infile_format: Input file format.
        args: Keyword or arguments to use in the call of the alignment tool, excluding infile and
            outfile arguments. By default, "default" arguments are used.
        outfile: Alignment output file.
        outfile_format: Output file format. By default, FASTA format.

    Raises:
        IOError: If the input path or the input file provided does not exist.
        RuntimeError: If the call to the alignment tool command raises an exception.

    * The input file format must be supported by Bio.SeqIO.
    * The output file format must be supported by Bio.AlignIO.
    """
    # Get the variables associated with the given alignment tool, or get those
    # values from **kwargs
    bin_path, bin_name = os.path.split(binary)
    bin_name = bin_name.lower()
    if bin_name in _TOOL_TO_LIB:
        tool_lib = _TOOL_TO_LIB[bin_name]
        sprt_infile_formats = tool_lib.SPRT_INFILE_FORMATS
        infile_cmd = tool_lib.INFILE_CMD
        keywords = tool_lib.KEYWORDS
    else:
        # Include the required variables through **kwargs dictionary
        sprt_infile_formats = kwargs['informats']
        infile_cmd = kwargs['incmd']
        keywords = dict()
    # Get the command line to run in order to get the resultant alignment
    infile_path = get_abspath(infile)
    # If the input file format is not supported by the alignment tool, convert it to a temporary
    # supported file
    if infile_format.lower() not in sprt_infile_formats:
        tmpfile = tempfile.NamedTemporaryFile()
        SeqIO.convert(infile_path, infile_format, tmpfile.name, sprt_infile_formats[0])
        infile_path = tmpfile.name
    # Get argument list from keyword dictionary or 'args' string
    if args in keywords:
        arg_list = keywords[args]
    else:
        # Remove possible empty strings in the given arguments
        arg_list = [arg for arg in args.split(' ')]
    # Create full command line list (removing empty elements)
    command = [x for x in [binary] + arg_list + [infile_cmd, infile_path] if x]
    # Run the alignment process handling any Runtime exception
    try:
        output = subprocess.run(command, stderr=subprocess.DEVNULL, universal_newlines=True, check=True,
                                stdout=subprocess.PIPE).stdout
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f'Running "{" ".join(e.cmd)}" raised an exception')
    else:
        alignment = AlignIO.read(StringIO(output), 'fasta')
        if outfile:
            # Save the resultant alignment in the given outfile and format
            outfile_path = get_abspath(outfile)
            AlignIO.write(alignment, outfile_path, outfile_format)
        # Return the resultant alignment as a Bio.Align.MultipleSeqAligment object
        return alignment

def main():
    """Default call for Align module."""     
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tool", required=True, help="Alignment tool")
    parser.add_argument("-i", "--input", required=True, help="FASTA file of unaligned sequences")
    parser.add_argument("-o", "--output", required=True, help="Output file name")
    args = parser.parse_args()
    alignment = get_alignment(args.tool, args.input, 'fasta', 'default', args.output)  
