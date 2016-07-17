#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  _Consense.py
# Last version :  v1.00 ( 02/Feb/2016 )
# Description :  MEvoLib's variables library functions to ease the usage of
#       Consense from PHYLIP
#       (<http://evolution.genetics.washington.edu/phylip.html>).
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  02/Feb/2016
#   VERSION :  v1.00
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

from __future__ import absolute_import

import os
import tempfile

from Bio import Phylo

from MEvoLib._utils import get_abspath


#-------------------------------------------------------------------------------

SPRT_INFILE_FORMATS = ['newick']

KEYWORDS = {'default': ['R', '2', 'Y']} # majority rule consensus


#-------------------------------------------------------------------------------

def gen_args ( args, infile_path, outfile ) :
    """
    Return the argument list generated from 'args' and the infile path
    requested.

    Arguments :
        args  ( string )
            Keyword or arguments to use in the call of Consense, excluding
            infile and outfile arguments.
        infile_path  ( string )
            Input alignment file path.
        outfile  ( string )
            Consensus tree output file.

    Returns :
        list
            List of arguments (excluding binary file) to call Consense.
    """
    if ( outfile ) :
        outfile_path = get_abspath(outfile)
    else :
        # Output files will be saved in temporary files to retrieve the
        # consensus tree
        outfile_path = os.path.join(tempfile.gettempdir(),
                                    tempfile.gettempprefix() + \
                                        next(tempfile._get_candidate_names()))
    # Create full command line list
    argument_list = [infile_path, outfile_path]
    return ( argument_list )



def gen_stdin_content ( args ) :
    """
    Arguments :
        args  ( string )
            Keyword or arguments to use in the call of Consense, excluding
            infile and outfile arguments. If 'args' is not a keyword, the second
            character will be used as separator of the different arguments.

    Returns :
        string
            Standard input content generated from 'args'.
    """
    if ( args in KEYWORDS ) :
        options = KEYWORDS[args]
    else : # args not in KEYWORDS
        options = [opt  for opt in args.split(args[1])]
    # If the output file already exists, overwritte it
    options.append('R')
    stdin_content = '\n'.join(options) + '\n'
    return ( stdin_content )



def get_results ( command ) :
    """
    Extract resultant consensus tree from the files generated during the
    execution of 'command'.

    Arguments :
        command  ( list )
            Consense's command line executed.

    Returns :
        Bio.Phylo.BaseTree
            Resultant consensus tree.

    Raises :
        IOError
            If the consensus tool didn't generate a consensus tree (indicated by
            user's options/arguments).
    """
    outfile_path = command[2]
    try :
        consensus_tree = Phylo.read(outfile_path, 'newick')
    except IOError as e :
        cleanup(command)
        message = 'The given arguments don\'t generate a consensus tree file'
        raise IOError(message)
    else :
        return ( consensus_tree )



def cleanup ( command ) :
    """
    Remove the temporary file created (if any) in gen_args() function.

    Arguments :
        command  ( list )
            Consense's command line executed.
    """
    logfile_path = get_abspath(command[2])
    if ( os.path.dirname(logfile_path) == tempfile.gettempdir() ) :
        os.remove(logfile_path)


#-------------------------------------------------------------------------------
