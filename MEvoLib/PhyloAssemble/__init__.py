#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  __init__.py
# Last version :  v1.0 ( 12/Feb/2016 )
# Description :  Functions aimed to provide an easy interface to handle
#       different phylogenetic assembly tools (supertrees, consensus trees,
#       ...).
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  12/Feb/2016
#   VERSION :  v1.0
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

from __future__ import absolute_import

import os
import tempfile
import subprocess

from . import _Consense

from MEvoLib._utils import get_abspath
from MEvoLib._py3k import viewkeys, viewitems, DEVNULL, StringIO


#-------------------------------------------------------------------------------

#_STREE_TOOL_TO_LIB = { 'superfine': _SuperFine }

_CONS_TOOL_TO_LIB = { 'consense': _Consense }


#-------------------------------------------------------------------------------

def get_keywords ( tool ) :
    """
    Arguments :
        tool  ( string )
            Name of the supertree or consensus tool.

    Returns :
        dict
            Dictionary containing the keywords and their corresponding
            arguments.

    Raises :
        ValueError
            If the tool introduced isn't included in MEvoLib.PhyloAssemble.
    """
    tool = tool.lower()
    #tool_lib_keys = viewkeys(_STREE_TOOL_TO_LIB) | viewkeys(_CONS_TOOL_TO_LIB)
    tool_lib_keys = viewkeys(_CONS_TOOL_TO_LIB)
    if ( tool not in tool_lib_keys ) :
        message = 'The tool "{}" isn\'t included in ' \
                  'MEvoLib.PhyloAssemble'.format(tool)
        raise ValueError(message)
    # else : # tool in tool_lib_keys
    keyword_dict = {}
#    if ( tool in _STREE_TOOL_TO_LIB ) :
#        tool_lib_dict = _STREE_TOOL_TO_LIB
#    else : # tool in _CONS_TOOL_TO_LIB
    if ( tool in _CONS_TOOL_TO_LIB ) :
        tool_lib_dict = _CONS_TOOL_TO_LIB
    for key, value in iter(viewitems(tool_lib_dict[tool].KEYWORDS)) :
        keyword_dict[key] = ' '.join(value)
    return ( keyword_dict )




def get_consensus_tree ( binary, infile, infile_format, args = 'default',
                         outfile = None, outfile_format = 'newick' ) :
    """
    Calculate the consensus tree of the input trees file with the given
    arguments. The resultant consensus tree is returned as a Bio.Phylo.BaseTree
    object and saved in the output file (if provided). If 'infile' or 'outfile'
    contain a relative path, the current working directory will be used to get
    the absolute path. If the output file already exists, the old file will be
    overwritten without any warning.

    Arguments :
        binary  ( string )
            Name or path of the consensus tool.
        infile  ( string )
            Input phylogenetic trees file.
        infile_format  ( string )
            Input file format.
        args  ( Optional[string] )
            Keyword or arguments to use in the call of the consensus tool,
            excluding infile and outfile arguments. By default, 'default'
            arguments are used.
            * For Consense, the second character will be used as separator of
            the different arguments. 
        outfile  ( Optional[string] )
            Consensus tree output file.
        outfile_format  ( Optional[string] )
            Output file format. By default, NEWICK format.

    Returns :
        Bio.Phylo.BaseTree
            Resultant consensus tree.

    Raises :
        ValueError
            If the tool introduced isn't included in MEvoLib.
        IOError
            If the input path or the input file provided doesn't exist.
        RuntimeError
            If the call to the phylogenetic inference tool command raises an
            exception.
        IOError
            If the consensus tool didn't generate a consensus tree (indicated by
            user's options/arguments).

    * The input file format must be supported by Bio.Phylo.
    * The output file format must be supported by Bio.Phylo.
    """
    # Get the variables associated with the given consensus tool
    bin_path, bin_name = os.path.split(binary)
    bin_name = bin_name.lower()
    if ( bin_name in _CONS_TOOL_TO_LIB ) :
        tool_lib = _CONS_TOOL_TO_LIB[bin_name]
        sprt_infile_formats = tool_lib.SPRT_INFILE_FORMATS
        gen_args = tool_lib.gen_args
        gen_stdin_content = tool_lib.gen_stdin_content
        get_results = tool_lib.get_results
        cleanup = tool_lib.cleanup
    else : # bin_name not in _CONS_TOOL_TO_LIB
        message = 'The consensus tool "{}" isn\'t included in ' \
                  'MEvoLib.PhyloAssemble'.format(bin_name)
        raise ValueError(message)
    # Get the command line to run in order to get the consensus tree
    infile_path = get_abspath(infile)
    # If the input file format is not supported by the consensus tool, convert
    # it to a temporary supported file
    if ( infile_format.lower() not in sprt_infile_formats ) :
        tmpfile = tempfile.NamedTemporaryFile()
        Phylo.convert(infile_path, infile_format, tmpfile.name,
                      sprt_infile_formats[0])
        infile_path = tmpfile.name
    # Create full command line list
    command = [binary] + gen_args(args, infile_path, outfile)
    # Generate the standard input file content
    stdin_content = gen_stdin_content(args)
    # Create the input file with the given options
    with tempfile.NamedTemporaryFile(mode='w+') as stdin_file :
        stdin_file.write(stdin_content)
        stdin_file.seek(0)
        # Run the consensus process handling any Runtime exception
        try :
            subprocess.check_call(command, stdin=stdin_file, stdout=DEVNULL,
                                  stderr=DEVNULL, universal_newlines=True)
        except subprocess.CalledProcessError as e :
            cleanup(command)
            message = 'Running "{}" raised an exception'.format(' '.join(e.cmd))
            raise RuntimeError(message)
        else :
            consensus_tree = get_results(command)
            cleanup(command)
            # Return the resultant consensus tree as a Bio.Phylo.BaseTree object
            return ( consensus_tree )


#-------------------------------------------------------------------------------
