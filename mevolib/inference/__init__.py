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
"""Functions aimed to provide an easy interface to handle different phylogenetic inference and bootstrapping tools."""

import argparse
import tempfile
import subprocess
from typing import Optional
from pathlib import Path

import Bio.Phylo.BaseTree
from Bio import AlignIO, Phylo

from . import _RAxML
from mevolib.inference import _FastTree, _RAxML

_PHYLO_TOOL_TO_LIB = {"fasttree": _FastTree, "raxml": _RAxML}

_BOOTS_TOOL_TO_LIB = {}


def get_tools() -> dict:
    """
    Returns :
        dict
            Dictionary of phylogenetic inference and bootstrapping software
            tools included in the current version of MEvoLib.
    """
    return dict(
        [("inference", list(_PHYLO_TOOL_TO_LIB.keys())), ("bootstrap", list(_BOOTS_TOOL_TO_LIB.keys()))]
    )


def get_keywords(tool: str) -> dict:
    """
    Arguments :
        tool: Name of the phylogenetic inference or bootstrapping tool.

    Returns :
        dict: Dictionary containing the keywords and their corresponding
            arguments.

    Raises :
        ValueError: If the tool introduced isn't included in MEvoLib.Inference.
    """
    tool = tool.lower()
    tool_lib_keys = _PHYLO_TOOL_TO_LIB.keys() | _BOOTS_TOOL_TO_LIB.keys()
    if tool not in tool_lib_keys:
        raise ValueError(f"The tool '{tool}' isn't included in 'MEvoLib.Inference'")
    # else : # tool in tool_lib_keys
    keyword_dict = {}
    if tool in _PHYLO_TOOL_TO_LIB:
        tool_lib_dict = _PHYLO_TOOL_TO_LIB
    else:
        tool_lib_dict = _BOOTS_TOOL_TO_LIB
    for key, value in tool_lib_dict.items():
        keyword_dict[key] = value
    return keyword_dict


def get_phylogeny(
    binary: str,
    infile: str,
    infile_format: str,
    args: Optional[str] = "default",
    outfile: Optional[str] = None,
    outfile_format: Optional[str] = "newick",
    bootstraps: Optional[int] = 0,
    tmp_file: Optional[str] = None,
    parent_file: Optional[str] = None,
    seed: Optional[int] = None,
) -> (Bio.Phylo.BaseTree, float):
    """
    Infers the phylogeny from the input alignment using the phylogenetic
    inference tool and arguments given. The resultant phylogeny is returned as a
    Bio.Phylo.BaseTree object and saved in the ouput file (if provided). If
    'infile' or 'outfile' contain a relative path, the current working directory
    will be used to get the absolute path. If the output file already exists,
    the old file will be overwritten without any warning.

    Arguments :
        binary: Name or path of the phylogenetic inference tool.
        infile: Sequence alignment file.
        infile_format: Input file format.
        args: Keyword or arguments to use in the call of the phylogenetic
            inference tool, excluding infile and outfile arguments. By default,
            'default' arguments are used.
        outfile: Phylogenetic tree output file.
        outfile_format: Output file format. By default, NEWICK format.
        bootstraps: Number of bootstraps to generate. By default, 0 (only use the
          input alignment).
        tmp_file: Path of the directory we want to save the selected tool's output data into (Only required while testing for
            reproducibility purposes or just if the user wants a specific folder to allocate the result of the execution
            into. Otherwise, a random one will be provided).
        parent_file: Path of the directory we want to save all the tmp_file's data into (Only required while testing for
            reproducibility purposes or just if the user wants a specific folder to allocate the result of the tmp_files
            into.).
        seed: Specify a random number seed for the parsimony inferences (Only required for while testing for
            reproducibility purposes. Otherwise, a random one will be provided).

    Returns :
        Bio.Phylo.BaseTree: Resultant phylogenetic tree.
        float: Log-likelihood score of the phylogeny.

    Raises :
        ValueError: If the tool introduced isn't included in MEvoLib.
        IOError: If the input path or the input file provided doesn't exist.
        RuntimeError: If the call to the phylogenetic inference tool command raises an
            exception.

    * The input file format must be supported by Bio.AlignIO.
    * The output file format must be supported by Bio.Phylo.
    """
    # Get the variables associated with the given phylogenetic inference tool
    bin_name = Path(binary).name
    bin_name = bin_name.lower()
    if bin_name in _PHYLO_TOOL_TO_LIB:
        tool_lib = _PHYLO_TOOL_TO_LIB[bin_name]
        sprt_infile_formats = tool_lib.SPRT_INFILE_FORMATS
        gen_args = tool_lib.gen_args
        get_results = tool_lib.get_results
        cleanup = tool_lib.cleanup
    else:
        raise ValueError(
            f'The phylogenetic inference tool "{bin_name}" isn\'t included in ' "MEvoLib.Inference"
        )
    # Get the command line to run in order to get the resultant phylogeny
    infile_path = Path(infile).absolute()
    # If the input file format is not supported by the phylogenetic inference
    # tool, convert it to a temporary supported file
    if infile_format.lower() not in sprt_infile_formats:
        if tmp_file is None:
            tmp_file = tempfile.NamedTemporaryFile()
            infile_path = tmp_file.name

        out_file = tmp_file.name

        AlignIO.convert(infile_path, infile_format, out_file, sprt_infile_formats[0])
        infile_path = Path(out_file).absolute()

    if parent_file is None:
        parent_file = tmp_file
    # Create full command line list
    if bin_name == "raxml":
        command = ["raxmlHPC"] + gen_args(args, infile_path, bootstraps, tmp_file, seed)
    else:
        command = [bin_name] + gen_args(args, infile_path, bootstraps, tmp_file)
    # Run the phylogenetic inference process handling any Runtime exception
    try:
        output = subprocess.run(
            command, stderr=subprocess.DEVNULL, universal_newlines=True, check=True, stdout=subprocess.PIPE
        ).stdout
    except subprocess.CalledProcessError as e:
        cleanup(command, parent_file)
        raise RuntimeError(f'Running "{e.cmd}" raised an exception')
    else:
        phylogeny, score = get_results(command, output)
        if outfile:
            # Save the resultant phylogeny in the given outfile and format
            outfile_path = Path(outfile).absolute()
            Phylo.write(phylogeny, outfile_path, outfile_format)
        cleanup(command, parent_file)
        # Return the resultant phylogeny as a Bio.Phylo.BaseTree object and its
        # log-likelihood score
        return phylogeny, score

def main():
    """Default call for Inference module."""
    parser = argparse.ArgumentParser(
        description="Infers the phylogeny of a given set of aligned genes using the given inference"
        + " tool and arguments"
    )
    parser.add_argument(
        "-t", "--tool", required=True, help="Name or path of the phylogenetic inference tool."
    )
    parser.add_argument("-i", "--input", required=True, help="Aligned sequences input file")
    parser.add_argument("-if", "--informat", required=False, help="Input file format (fasta, by default)")
    parser.add_argument(
        "-a",
        "--args",
        required=False,
        help="Keyword or arguments to use in the call of the phylogenetic"
        + " inference tool, excluding infile and outfile arguments. By default, 'default ' arguments are used.",
    )
    parser.add_argument("-o", "--output", required=True, help="Output file name (without extension)")
    parser.add_argument("-of", "--outformat", required=False, help="Output file format. By default, 'NEWICK' format.")
    parser.add_argument(
        "-b",
        "--bootstraps",
        required=False,
        help="Number of bootstraps to generate. By default, 0 (only uses the input alignment).",
    )
    args = parser.parse_args()

    # We split the input file to obtain its name
    filename = Path(args.input).stem
    out_format = "newick" if args.outformat is None else args.outformat
    get_phylogeny(
        binary=args.tool,
        infile=args.input,
        infile_format="fasta" if args.informat is None else args.informat,
        args="default" if args.args is None else args.args,
        outfile=f"{filename}_inference.{out_format}",
        outfile_format=out_format,
    )
