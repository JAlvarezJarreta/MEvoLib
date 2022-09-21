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
# Last version :  v1.10 ( 16/Jul/2016 )
# Description :  Functions aimed to provide an easy interface to handle
#       different clustering methods.
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  16/Jul/2016
#   VERSION :  v1.10
#   AUTHOR(s) :  J. Alvarez-Jarreta
#   CHANGES :  * Added get_tools() method to have an easy access to all the
#                  available clustering methods and tools.
#              * Fixed a problem with get_subsets() method that forced to pass
#                  all the arguments (including None values) to every method.
#                  Now it accepts keyworded arguments.
#
#   DATE :  02/Feb/2016
#   VERSION :  v1.00
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

from __future__ import absolute_import, division

import multiprocessing
import math

from Bio import SeqIO

from . import Genes
from . import NaiveRows
from . import NaiveCols
from . import PRD

from mevolib._utils import get_abspath, NUMCORES
from mevolib._py3k import viewkeys


#-------------------------------------------------------------------------------

_METHOD_TO_FUNC = { 'genes': Genes.map_seqs,
                    'rows': NaiveRows.map_seqs,
                    'cols': NaiveCols.map_seqs,
                    'prd': PRD.map_seqs }


#-------------------------------------------------------------------------------

def get_tools ( ) :
    """
    Returns :
        list
            List of clustering methods and software tools included in the
            current version of MEvoLib.
    """
    return ( list(viewkeys(_METHOD_TO_FUNC)) )



def get_subsets ( method, seqfile, fileformat = 'genbank', *args, **kwargs ) :
    """
    Division of all the sequences stored in the sequence input file into subsets
    applying the 'method' function. If 'seqfile' contains a relative path, the
    current working directory will be used to get the absolute path.

    Arguments :
        method  ( string )
            Desired partition method (case-insensitive): genes, naive rows or
            cols, padded-Recursive-DCM3.
        seqfile  ( string )
            Input sequences file.
        fileformat  ( string ) 
            Input file format.
        args & kwargs
            Non-keyworded and keyworded arguments passed to the selected method.
            

    Returns :
        dict
            Dictionary with the set identifiers as keys and the corresponding
            sequences as values in lists of SeqRecord objects.

    Raises :
        ValueError
            If there is no corresponding method to 'method' value.
        IOError
            If the path or the file provided doesn't exist.
        IOError
            If the file format provided doesn't correspond to the actual one.

    * The file format must be supported by Bio.SeqIO.
    * For "rows" method, if the number of input sequences is lower than the
    number of sets multiplied by the number of cores, the resulting sets might
    be fewer than the number requested.
    """
    method_key = method.lower()
    if ( method_key not in _METHOD_TO_FUNC ) :
        message = 'The method "{}" isn\'t included in ' \
                  'MEvoLib.Cluster'.format(method)
        raise ValueError(message)
    # else : # method_key in _METHOD_TO_FUNC
    # Get the mapping function and the sequence file path
    mapseqs_func = _METHOD_TO_FUNC[method_key]
    filepath = get_abspath(seqfile)
    if ( method_key in ['prd', 'genes'] ) :
        # Non data-driven (throught input slicing) parallelizable methods
        seq_list = (x  for x in SeqIO.parse(filepath, fileformat))
        set_dict = mapseqs_func(seq_list, *args, **kwargs)
    else :
        # Data-driven (throught input slicing) parallelizable methods
        manager = multiprocessing.Manager()
        seq_list = manager.list([x  for x in SeqIO.parse(filepath, fileformat)])
        num_seqs = len(seq_list)
        # Launch one process per available CPU core
        slice_size = int(math.ceil(num_seqs / NUMCORES))
        pool = multiprocessing.Pool(processes=NUMCORES)
        results = [pool.apply_async(mapseqs_func,
                                    args=(seq_list[start:start+slice_size],) + \
                                          args)
                       for start in range(0, num_seqs, slice_size)]
        # Build the final sets dictionary merging the results of all executed
        # processes
        output = [p.get() for p in iter(results)]
        set_dict = output[0]
        for key in iter(set_dict) :
            for result in output[1:] :
                set_dict[key].extend(result[key])
    return ( set_dict )


#-------------------------------------------------------------------------------
