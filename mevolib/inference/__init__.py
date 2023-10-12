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


from __future__ import absolute_import

import os
import tempfile
import subprocess
from typing import Optional
from io import StringIO
from pathlib import Path

import Bio.Phylo.BaseTree
from Bio import AlignIO, Phylo

from . import _RAxML
from mevolib.inference import _FastTree, _RAxML
from mevolib._utils import get_abspath

_PHYLO_TOOL_TO_LIB = { 'fasttree': _FastTree,
                       'raxml': _RAxML }

_BOOTS_TOOL_TO_LIB = { }


def get_tools() -> dict:
    """
    Returns :
        dict
            Dictionary of phylogenetic inference and bootstrapping software
            tools included in the current version of MEvoLib.
    """
    return ( dict([('inference', list(_PHYLO_TOOL_TO_LIB.keys())),
                   ('bootstrap', list(_BOOTS_TOOL_TO_LIB.keys()))]))



def get_keywords (tool: str) -> dict:
    """
    Arguments :
        tool  ( string )
            Name of the phylogenetic inference or bootstrapping tool.

    Returns :
        dict
            Dictionary containing the keywords and their corresponding
            arguments.

    Raises :
        ValueError
            If the tool introduced isn't included in MEvoLib.Inference.
    """
    tool = tool.lower()
    tool_lib_keys = _PHYLO_TOOL_TO_LIB.keys() | _BOOTS_TOOL_TO_LIB.keys()
    if ( tool not in tool_lib_keys ) :
        raise ValueError('The tool "{}" isn\'t included in ' \
                         'MEvoLib.Inference'.format(tool))
    # else : # tool in tool_lib_keys
    keyword_dict = {}
    if ( tool in _PHYLO_TOOL_TO_LIB ) :
        tool_lib_dict = _PHYLO_TOOL_TO_LIB
    else : # tool in _BOOTS_TOOL_TO_LIB
        tool_lib_dict = _BOOTS_TOOL_TO_LIB
    for key,value in tool_lib_dict.items() :
        keyword_dict[key] = ' '.join(value)
    return ( keyword_dict ) 



def get_phylogeny(binary: str, infile: str, infile_format: str, args: Optional[str] = 'default',
                    outfile : Optional[str] = None , outfile_format : Optional[str] = 'newick',
                    bootstraps: Optional[int] = 0) -> (Bio.Phylo.BaseTree, float):
    """
    Infer the phylogeny from the input alignment using the phylogenetic
    inference tool and arguments given. The resultant phylogeny is returned as a
    Bio.Phylo.BaseTree object and saved in the ouput file (if provided). If
    'infile' or 'outfile' contain a relative path, the current working directory
    will be used to get the absolute path. If the output file already exists,
    the old file will be overwritten without any warning.

    Arguments :
        binary  ( string )
            Name or path of the phylogenetic inference tool.
        infile  ( string )
            Sequence alignment file.
        infile_format  ( string )
            Input file format.
        args  ( Optional[string] )
            Keyword or arguments to use in the call of the phylogenetic
            inference tool, excluding infile and outfile arguments. By default,
            'default' arguments are used.
        outfile  ( Optional[string] )
            Phylogenetic tree output file.
        outfile_format  ( Optional[string] )
            Output file format. By default, NEWICK format.
        bootstraps  ( Optional[int] )
            Number of bootstraps to generate. By default, 0 (only use the input
            alignment).

    Returns :
        Bio.Phylo.BaseTree
            Resultant phylogenetic tree.
        float
            Log-likelihood score of the phylogeny.

    Raises :
        ValueError
            If the tool introduced isn't included in MEvoLib.
        IOError
            If the input path or the input file provided doesn't exist.
        RuntimeError
            If the call to the phylogenetic inference tool command raises an
            exception.

    * The input file format must be supported by Bio.AlignIO.
    * The output file format must be supported by Bio.Phylo.
    """
    # Get the variables associated with the given phylogenetic inference tool
    bin_name = Path.split(binary)
    bin_name = bin_name.lower()
    if ( bin_name in _PHYLO_TOOL_TO_LIB ) :
        tool_lib = _PHYLO_TOOL_TO_LIB[bin_name]
        sprt_infile_formats = tool_lib.SPRT_INFILE_FORMATS
        gen_args = tool_lib.gen_args
        get_results = tool_lib.get_results
        cleanup = tool_lib.cleanup
    else : # bin_name not in _PHYLO_TOOL_TO_LIB
        message = 'The phylogenetic inference tool "{}" isn\'t included in ' \
                  'MEvoLib.Inference'.format(bin_name)
        raise ValueError(message)
    # Get the command line to run in order to get the resultant phylogeny
    infile_path = get_abspath(infile)
    # If the input file format is not supported by the phylogenetic inference
    # tool, convert it to a temporary supported file
    if ( infile_format.lower() not in sprt_infile_formats ) :
        tmpfile = tempfile.NamedTemporaryFile()
        AlignIO.convert(infile_path, infile_format, tmpfile.name,
                        sprt_infile_formats[0])
        infile_path = tmpfile.name
    # Create full command line list
    command = [binary] + gen_args(args, infile_path, bootstraps)
    # Run the phylogenetic inference process handling any Runtime exception
    try :
        output = subprocess.run(command, stderr=subprocess.DEVNULL,
                                         universal_newlines=True)
    except subprocess.CalledProcessError as e :
        cleanup(command)
        message = 'Running "{}" raised an exception'.format(' '.join(e.cmd))
        raise RuntimeError(message)
    else :
        phylogeny, score = tool_lib.get_results(command, output)
        if outfile:
            # Save the resultant phylogeny in the given outfile and format
            outfile_path = get_abspath(outfile)
            Phylo.write(phylogeny, outfile_path, outfile_format)
        cleanup(command)
        # Return the resultant phylogeny as a Bio.Phylo.BaseTree object and its
        # log-likelihood score
        return phylogeny, score

