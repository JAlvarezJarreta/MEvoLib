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
"""MEvoLib's variables and library functions to ease the usage of RAxML."""

import os
import random
import tempfile
import shutil
from io import StringIO
from pathlib import Path
from typing import Optional

import Bio.Phylo.BaseTree
from Bio import Phylo

from mevolib._utils import NUMCORES

SPRT_INFILE_FORMATS = ["fasta", "phylip"]

KEYWORDS = {
    "default": ["-m", "GTRCAT", "--silent"],
    # Nucleotide evolution models
    "JC+CAT": ["-m", "GTRCAT", "--JC69", "--silent"],
    "JC+G": ["-m", "GTRGAMMA", "--JC69", "--silent"],
    "K80+CAT": ["-m", "GTRCAT", "--K80", "--silent"],
    "K80+G": ["-m", "GTRGAMMA", "--K80", "--silent"],
    "HKY+CAT": ["-m", "GTRCAT", "--HKY85", "--silent"],
    "HKY+G": ["-m", "GTRGAMMA", "--HKY85", "--silent"],
    "GTR+CAT": ["-m", "GTRCAT", "--silent"],
    "GTR+G": ["-m", "GTRGAMMA", "--silent"],
    # Protein evolution models
    "JTT+CAT": ["-m", "PROTCATJTT", "--silent"],
    "JTT+G": ["-m", "PROTGAMMAJTT", "--silent"],
    "WAG+CAT": ["-m", "PROTCATWAG", "--silent"],
    "WAG+G": ["-m", "PROTGAMMAWAG", "--silent"],
}


def gen_args(
    args: str,
    infile_path: str,
    bootstraps: int,
    tmpdir_path: Optional[str] = None,
    seed: Optional[int] = None,
) -> list:
    """
    Return the argument list generated from 'args', the infile path and the
    bootstraps requested.

    Arguments :
        args: Keyword or arguments to use in the call of RAxML, excluding infile and outfile arguments ("-s" and "-n").
        infile_path: Input alignment file path.
        bootstraps: Number of bootstraps to generate.
        tmpdir_path: Path of the directory we want to save the RAxML output data into (Only required while testing for
            reproducibility purposes or just if the user wants a specific folder to allocate the result of the execution
            into. Otherwise, a random one will be provided).
        seed: Specify a random number seed for the parsimony inferences (Only required for while testing for
            reproducibility purposes. Otherwise, a random one will be provided).

    Returns :
        list
            List of arguments (excluding binary file) to call RAxML.
    """
    argument_list = []
    if args in KEYWORDS:
        if seed is None:
            seed = random.randint(1, 999999)
        pid = os.getpid()
        if tmpdir_path is None:
            tmpdir_path = tempfile.mkdtemp()
        argument_list = ["-p", str(seed), "-n", str(pid), "-T", str(NUMCORES), "-w", tmpdir_path] + KEYWORDS[
            args
        ]
    else:  # args not in KEYWORDS
        argument_list = [arg for arg in args.split(" ")]
    # Add the bootstrapping generation option if 'boostraps' is greater than 0
    if bootstraps > 0:
        argument_list += ["-N", str(bootstraps)]
    # Add the input file option
    argument_list += ["-s", infile_path]
    return argument_list


def get_results(command: list, output: str) -> tuple[Bio.Phylo.BaseTree.Tree, float]:
    """
    Extract resultant phylogeny and its log-likelihood score from 'output' and
    files generated during the execution of 'command'.

    Arguments :
        command: RAxML's command line executed.
        output: Output from 'command' execution.

    Returns :
        Bio.Phylo.BaseTree.Tree: Resultant phylogenetic tree.
        float: Log-likelihood score of the phylogeny.
    """
    index = command.index("-n") + 1
    outfiles_id = command[index]
    # Get the output directory from 'command' or use current working directory
    if "-w" in command:
        index = command.index("-w") + 1
        outdir_path = command[index]
    else:  # '-w' not in command
        outdir_path = os.getcwd()
    # Read the final phylogeny from bestTree file
    treefile_path = Path.joinpath(outdir_path, "RAxML_bestTree." + outfiles_id)
    phylogeny = Phylo.read(treefile_path, "newick")
    # Read the info file to get the log-likelihood score of the final phylogeny
    infofile_path = Path.joinpath(outdir_path, "RAxML_info." + outfiles_id)
    with open(infofile_path, "r") as infofile:
        # It is located at the last line that matches "TreeLogLk.*" pattern
        for line in infofile.readlines():
            if "GAMMA-based Score" in line:
                score = float(line.split(" ")[-1])
    return (phylogeny, score)


def cleanup(command, tmpdir_path: Optional[str] = None):
    """
    Remove the temporary files and directories created (if any) in gen_args()
    function.

    Arguments :
        command: RAxML's command line executed.
        tmpdir_path: Path of the directory we may have saved the RAxMl output data into, got from gen_args function.
            (Only required while testing for reproducibility purposes or just if the user wants a specific folder to
            allocate the result of the execution into. Otherwise, a random one will be provided).
    """
    if "-w" in command:
        index = command.index("-w") + 1
        outdir_path = command[index]
        if tmpdir_path is None:
            dir_name = tempfile.gettempdir()
        else:
            dir_name = tmpdir_path
        if (Path(outdir_path).parent == dir_name) and Path.exists(outdir_path):
            os.remove(outdir_path)
