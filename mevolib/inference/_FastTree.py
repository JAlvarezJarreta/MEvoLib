#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  _FastTree.py
# Last version :  v2.0 ( 25/Sep/2023 )
# Description :  MEvoLib's variables and library functions to ease the usage of
#       FastTree.
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  25/Sep/2023
#   VERSION :  v2.0
#   AUTHOR(s) :  S. Moragon Jimenez
#   CHANGES :  * Added StringIO and Bio.Phylo.BaseTree imports as part of the 
#                   module version update to Python 3.10.6.
#              * Added input variables and return type in methods, so that
#                   it is easier to visualize the whole IO parameters of the
#                   functions.
#
#   DATE :  26/Jan/2016
#   VERSION :  v1.01
#   AUTHOR(s) :  J. Alvarez-Jarreta
#   CHANGES :  * Fixed a bug in gen_args() method that produced a modification
#                  of the KEYWORDS dictionary after the first usage of the
#                  function that lasted until the reload of the module.
#              * Fixed a bug in cleanup() method that raised an error if the
#                  temporal log file was already deleted.
#
#   DATE :  26/Jan/2016
#   VERSION :  v1.00
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

from __future__ import absolute_import

import os
import tempfile
import Bio.Phylo.BaseTree

from Bio import Phylo


from mevolib._utils import get_abspath
from io import StringIO


#-------------------------------------------------------------------------------

SPRT_INFILE_FORMATS = ['fasta', 'phylip']

KEYWORDS = {'default': ['-gtr', '-nt', '-nopr', '-quiet'],
            # Nucleotide evolution models
            'GTR+CAT': ['-gtr', '-nt', '-nopr', '-quiet'],
            'GTR+G':   ['-gtr', '-nt', '-nocat', '-gamma', '-nopr', '-quiet'],
            'JC+CAT':  ['-nt', '-nopr', '-quiet'],
            'JC+G':    ['-nt', '-nocat', '-gamma', '-nopr', '-quiet'],
            # Protein evolution models
            'JTT+CAT': ['-nopr', '-quiet'],
            'WAG+CAT': ['-wag', '-nopr', '-quiet']}


#-------------------------------------------------------------------------------

def gen_args (args: str, infile_path: str, bootstraps: int) -> list:
    """
    Return the argument list generated from 'args', the infile path and the
    bootstraps requested.

    Arguments :
        args  ( string )
            Keyword or arguments to use in the call of FastTree, excluding
            infile and outfile arguments.
        infile_path  ( string )
            Input alignment file path.
        bootstraps  ( int )
            Number of bootstraps to generate.

    Returns :
        list
            List of arguments (excluding binary file) to call FastTree.
    """
    if ( args in KEYWORDS ) :
        argument_list = list(KEYWORDS[args])
    else : # args not in KEYWORDS
        argument_list = [arg  for arg in args.split(' ')]
    # Add the log file generation to get the log-likelihood score of the
    # resultant phylogeny
    if ( '-log' not in argument_list ) :
        log_tmpfile = tempfile.NamedTemporaryFile(delete=False)
        argument_list += ['-log', log_tmpfile.name]
    # Add the bootstrapping generation option if 'boostraps' is greater than 0
    if ( bootstraps > 0 ) :
        argument_list += ['-boot', str(bootstraps)]
    # Add the input file option
    argument_list += [infile_path]
    return ( argument_list )



def get_results (command: list, output: str) -> Bio.Phylo.BaseTree:
    """
    Extract resultant phylogeny and its log-likelihood score from 'output' and
    files generated during the execution of 'command'.

    Arguments :
        command  ( list )
            FastTree's command line executed.
        output  ( string )
            Output from 'command' execution.

    Returns :
        Bio.Phylo.BaseTree
            Resultant phylogenetic tree.
        float
            Log-likelihood score of the phylogeny.
    """
    phylogeny = Phylo.read(StringIO(output), 'newick')
    # Read the log file to get the log-likelihood score of the final phylogeny
    index = command.index('-log') + 1
    logfile_path = get_abspath(command[index])
    with open(logfile_path, 'r') as logfile :
        # It is located at the last line that matches "TreeLogLk.*" pattern
        for line in logfile.readlines() :
            if ( 'TreeLogLk' in line ) :
                score = float(line.split('\t')[2])
    return ( phylogeny, score )



def cleanup (command: list) -> None:
    """
    Remove the temporary files and directories created (if any) in gen_args()
    function.

    Arguments :
        command  ( list )
            FastTree's command line executed.
    """
    index = command.index('-log') + 1
    logfile_path = get_abspath(command[index])
    if ( (os.path.dirname(logfile_path) == tempfile.gettempdir()) and
         os.path.lexists(logfile_path) ) :
        os.remove(logfile_path)


#-------------------------------------------------------------------------------
