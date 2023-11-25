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
"""Functions aimed to provide an easy interface to handle alignment tools."""

# __all__ = ["get_tools", "get_keywords", "get_alignment"]

from collections import namedtuple
from io import StringIO
from pathlib import Path
import tempfile
from typing import Optional
import subprocess

import Bio
from Bio import SeqIO, AlignIO

from mevolib._utils import NUMCORES


_ALIGN_TOOL = namedtuple("AlignTool", ["input_formats", "input_arg", "keywords"])
_CONFIGURED_TOOLS = {
    "clustalo": _ALIGN_TOOL(
        ["fasta", "clustal", "msf", "phylip", "selec", "stockholm"],
        "-i",
        {
            "default": ["--auto", "--output-order=input-order"],
        },
    ),
    "mafft": _ALIGN_TOOL(
        ["fasta"],
        "",
        {
            "default": ["--auto", "--quiet", "--thread", str(NUMCORES)],
            "linsi": ["--localpair", "--maxiterate", "1000", "--quiet", "--thread", str(NUMCORES)],
        },
    ),
    "muscle": _ALIGN_TOOL(
        ["fasta"],
        "-in",
        {
            "default": ["-quiet"],
            "fastdna": ["-maxiters", "1", "-diags", "-quiet"],
            "fastprot": ["-maxiters", "1", "-diags", "-sv", "-distance1", "kbit20_3", "-quiet"],
            "largein": ["-maxiters", "2", "-quiet"],
        },
    ),
}


def get_tools() -> list:
    """Returns the list of alignment software tools pre-configured in MEvoLib."""
    return list(_CONFIGURED_TOOLS.keys())


def get_keywords(tool: str) -> dict:
    """Returns the available keywords and corresponding arguments for the given align tool.

    If `tool` has not been pre-configured, an empty dictionary will be returned.

    Args:
        tool: Name of the alignment tool.

    """
    try:
        keyword_dict = {key: " ".join(val) for key, val in _CONFIGURED_TOOLS[tool.lower()].KEYWORDS.items()}
    except KeyError:
        keyword_dict = {}
    return keyword_dict


def get_alignment(
    binary: Path,
    infile: Path,
    input_format: str,
    args: str = "default",
    outfile: Optional[Path] = None,
    outfile_format: str = "fasta",
    **kwargs: dict,
) -> Bio.Align.MultipleSeqAlignment:
    """Returns the result of aligning the sequences using the given alignment tool and arguments.

    The resultant alignment is returned as a ``Bio.Align.MultipleSeqAlign`` object and saved in the output
    file (if provided). If `infile` or `outfile` contains a relative path, the current working directory
    will be used to get the absolute path. If the output file already exists, the old file will be
    overwritten without any warning.

    The alignment tool might not be included in MEvoLib, but it can still be used passing in `**kwargs`
    the keys "input_formats" and "input_arg" with the list of supported input formats and the input argument
    (empty string if it does not need any), respectively.

    Args:
        binary: Name or path of the alignment tool.
        infile: Unaligned input sequence file.
        input_format: Input file format.
        args: Keyword or arguments to use in the call of the alignment tool, excluding infile and
            outfile arguments. By default, "default" arguments are used.
        outfile: Alignment output file.
        outfile_format: Output file format. By default, FASTA format.

    """
    # Get the variables associated with the given alignment tool, or get those values from **kwargs
    bin_name = Path(binary).name.lower()
    bin_config = _CONFIGURED_TOOLS.get(
        bin_name,
        _ALIGN_TOOL(kwargs["input_formats"], kwargs["input_arg"], {})
    )
    # Prepare command line arguments
    infile_path = infile.absolute().resolve()
    if input_format.lower() not in bin_config.input_formats:
        # If the input format is not supported by the alignment tool, convert it to a supported one
        tmpfile = tempfile.NamedTemporaryFile()
        SeqIO.convert(infile_path, input_format, tmpfile.name, bin_config.input_formats[0])
        infile_path = tmpfile.name
    arg_list = bin_config.keywords.get(args, list(args.split(" ")))
    # Build command list, removing any empty elements
    command = [x for x in [binary] + arg_list + [bin_config.input_arg, infile_path] if x]
    output = subprocess.run(
        command, stderr=subprocess.DEVNULL, universal_newlines=True, check=True, stdout=subprocess.PIPE
    ).stdout
    alignment = AlignIO.read(StringIO(output), "fasta")
    if outfile:
        # Save the resultant alignment in the given outfile and format
        outfile_path = outfile.absolute().resolve()
        AlignIO.write(alignment, outfile_path, outfile_format)
    return alignment
