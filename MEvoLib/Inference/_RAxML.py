#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  _RAxML.py
# Last version :  v1.0 ( 26/Jan/2016 )
# Description :  MEvoLib's variables and library functions to ease the usage of
#       RAxML.
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  26/Jan/2016
#   VERSION :  v1.0
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

from __future__ import absolute_import

import os
import random
import tempfile
import shutil

from Bio import Phylo

from MEvoLib._utils import NUMCORES
from MEvoLib._py3k import StringIO


#-------------------------------------------------------------------------------

SPRT_INFILE_FORMATS = ['fasta', 'phylip']

KEYWORDS = {'default': ['-m', 'GTRCAT', '--silent'],
            # Nucleotide evolution models
            'JC+CAT':  ['-m', 'GTRCAT', '--JC69', '--silent'],
            'JC+G':    ['-m', 'GTRGAMMA', '--JC69', '--silent'],
            'K80+CAT': ['-m', 'GTRCAT', '--K80', '--silent'],
            'K80+G':   ['-m', 'GTRGAMMA', '--K80', '--silent'],
            'HKY+CAT': ['-m', 'GTRCAT', '--HKY85', '--silent'],
            'HKY+G':   ['-m', 'GTRGAMMA', '--HKY85', '--silent'],
            'GTR+CAT': ['-m', 'GTRCAT', '--silent'],
            'GTR+G':   ['-m', 'GTRGAMMA', '--silent'],
            # Protein evolution models
            'JTT+CAT': ['-m', 'PROTCATJTT', '--silent'],
            'JTT+G':   ['-m', 'PROTGAMMAJTT', '--silent'],
            'WAG+CAT': ['-m', 'PROTCATWAG', '--silent'],
            'WAG+G':   ['-m', 'PROTGAMMAWAG', '--silent']}


#-------------------------------------------------------------------------------

def gen_args ( args, infile_path, bootstraps ) :
    """
    Return the argument list generated from 'args', the infile path and the
    bootstraps requested.

    Arguments :
        args  ( string )
            Keyword or arguments to use in the call of RAxML, excluding infile
            and outfile arguments ("-s" and "-n").
        infile_path  ( string )
            Input alignment file path.
        bootstraps  ( int )
            Number of bootstraps to generate.

    Returns :
        list
            List of arguments (excluding binary file) to call RAxML.
    """
    argument_list = []
    if ( args in KEYWORDS ) :
        seed = random.randint(1, 999999)
        pid = os.getpid()
        tmpdir_path = tempfile.mkdtemp()
        argument_list = ['-p', str(seed), '-n', str(pid), '-T', str(NUMCORES),
                         '-w', tmpdir_path] + KEYWORDS[args]
    else : # args not in KEYWORDS
        argument_list = [arg  for arg in args.split(' ')]
    # Add the bootstrapping generation option if 'boostraps' is greater than 0
    if ( bootstraps > 0 ) :
        argument_list += ['-N', str(bootstraps)]
    # Add the input file option
    argument_list += ['-s', infile_path]
    return ( argument_list )



def get_results ( command, output ) :
    """
    Extract resultant phylogeny and its log-likelihood score from 'output' and
    files generated during the execution of 'command'.

    Arguments :
        command  ( list )
            RAxML's command line executed.
        output  ( string )
            Output from 'command' execution.

    Returns :
        Bio.Phylo.BaseTree
            Resultant phylogenetic tree.
        float
            Log-likelihood score of the phylogeny.
    """
    index = command.index('-n') + 1
    outfiles_id = command[index]
    # Get the output directory from 'command' or use current working directory
    if ( '-w' in command ) :
        index = command.index('-w') + 1
        outdir_path = command[index]
    else : # '-w' not in command
        outdir_path = os.getcwd()
    # Read the final phylogeny from bestTree file
    treefile_path = os.path.join(outdir_path, 'RAxML_bestTree.' + outfiles_id)
    phylogeny = Phylo.read(treefile_path, 'newick')
    # Read the info file to get the log-likelihood score of the final phylogeny
    infofile_path = os.path.join(outdir_path, 'RAxML_info.' + outfiles_id)
    with open(infofile_path, 'r') as infofile :
        # It is located at the last line that matches "TreeLogLk.*" pattern
        for line in infofile.readlines() :
            if ( 'GAMMA-based Score' in line ) :
                score = float(line.split(' ')[-1])
    return ( phylogeny, score )



def cleanup ( command ) :
    """
    Remove the temporary files and directories created (if any) in gen_args()
    function.

    Arguments :
        command  ( list )
            RAxML's command line executed.
    """
    if ( '-w' in command ) :
        index = command.index('-w') + 1
        outdir_path = command[index]
        print(os.path.dirname(outdir_path))
        if ( os.path.dirname(outdir_path) == tempfile.gettempdir() ) :
            shutil.rmtree(outdir_path)


#-------------------------------------------------------------------------------
