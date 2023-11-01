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
"""MEvoLib's variables and library functions to ease the usage of FastTree."""

import os
import tempfile
from io import StringIO
from pathlib import Path
from typing import Optional

import Bio.Phylo.BaseTree
from Bio import Phylo

SPRT_INFILE_FORMATS = ["fasta", "phylip"]

KEYWORDS = {
    "default": ["-gtr", "-nt", "-nopr", "-quiet"],
    # Nucleotide evolution models
    "GTR+CAT": ["-gtr", "-nt", "-nopr", "-quiet"],
    "GTR+G": ["-gtr", "-nt", "-nocat", "-gamma", "-nopr", "-quiet"],
    "JC+CAT": ["-nt", "-nopr", "-quiet"],
    "JC+G": ["-nt", "-nocat", "-gamma", "-nopr", "-quiet"],
    # Protein evolution models
    "JTT+CAT": ["-nopr", "-quiet"],
    "WAG+CAT": ["-wag", "-nopr", "-quiet"],
}


def gen_args(args: str, infile_path: str, bootstraps: int, log_tmpfile: Optional[str] = None) -> list:
    """
    Return the argument list generated from 'args', the infile path and the
    bootstraps requested.

    Arguments :
        args: Keyword or arguments to use in the call of FastTree, excluding
            infile and outfile arguments.
        infile_path: Input alignment file path.
        bootstraps: Number of bootstraps to generate.
        log_tmpfile: Path of the directory we want to save the FastTree output data into (Only required while testing for
            reproducibility purposes or just if the user wants a specific folder to allocate the result of the execution
            into. Otherwise, a random one will be provided).

    Returns :
        list: List of arguments (excluding binary file) to call FastTree.
    """
    if args in KEYWORDS:
        argument_list = list(KEYWORDS[args])
    else:
        argument_list = args.split(" ")
    # Add the log file generation to get the log-likelihood score of the
    # resultant phylogeny
    if "-log" not in argument_list:
        if log_tmpfile is None:
            log_tmpfile = tempfile.NamedTemporaryFile(delete=False).name
        argument_list += ["-log", log_tmpfile]
    # Add the bootstrapping generation option if 'boostraps' is greater than 0
    if bootstraps > 0:
        argument_list += ["-boot", str(bootstraps)]
    # Add the input file option
    argument_list += [infile_path]
    return argument_list


def get_results(command: list, output: str) -> tuple[Bio.Phylo.BaseTree, float]:
    """
    Extract resultant phylogeny and its log-likelihood score from 'output' and
    files generated during the execution of 'command'.

    Arguments :
        command: FastTree's command line executed.
        output: Output from 'command' execution.

    Returns :
        Bio.Phylo.BaseTree: Resultant phylogenetic tree.
        float: Log-likelihood score of the phylogeny.
    """

    phylogeny = Phylo.read(StringIO(output), "newick")
    # Read the log file to get the log-likelihood score of the final phylogeny
    index = command.index("-log") + 1
    logfile_path = Path(command[index])
    with logfile_path.open("r") as logfile:
        # It is located at the last line that matches "TreeLogLk.*" pattern
        for line in logfile.readlines():
            if "TreeLogLk" in line:
                score = float(line.split("\t")[2])
    return (phylogeny, score)


def cleanup(command: list, log_tmpfile: Optional[str] = None) -> None:
    """
    Remove the temporary files and directories created (if any) in gen_args()
    function.

    Arguments :
        command: FastTree's command line executed.
        log_tmpfile: Path of the directory we may have saved the FastTree output data into, got from gen_args function.
            (Only required while testing for reproducibility purposes or just if the user wants a specific folder to
            allocate the result of the execution into. Otherwise, a random one will be provided).
    """
    index = command.index("-log") + 1
    logfile_path = Path(command[index]).absolute()
    if log_tmpfile is None:
        dir_name = tempfile.gettempdir()
    else:
        dir_name = log_tmpfile
    if (Path(logfile_path).parent == dir_name) and Path.exists(logfile_path):
        os.remove(logfile_path)
