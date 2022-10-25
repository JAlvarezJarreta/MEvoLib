#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  PRD.py
# Last version :  v1.01 ( 19/Jul/2016 )
# Description :  Clustering with padded-Recursive-DMC3 decomposition (PRD) from
#       DACTAL system (<http://www.cs.utexas.edu/~phylo/software/dactal/>).
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  19/Jul/2016
#   VERSION :  v1.01
#   AUTHOR(s) :  J. Alvarez-Jarreta
#   CHANGES :  * Minor error in one exception output corrected.
#
#   DATE :  20/Jan/2016
#   VERSION :  v1.00
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

from __future__ import absolute_import
from __future__ import annotations

import os
import tempfile
import multiprocessing
import subprocess
import dendropy
 
from Bio import Phylo
from pathlib import Path

from mevolib._utils import NUMCORES
from mevolib._py3k import DEVNULL


#-------------------------------------------------------------------------------

def _prd_decomposition ( tree_file, subset_size, overlapping, binary = 'dcm' ) :
    """
    Apply the padded-Recursive-DMC3 decomposition (PRD) from DACTAL system to
    the tree file with the given overlapping. The sets that are over the subset
    size are marked to extract their corresponding subtree and process it
    recursively until all the decomposition sets are below the size threshold.

    Arguments :
        tree_file  ( string )
            Input tree file.
        subset_size  ( int )
            Maximum subset size.
        overlapping  ( string )
            Number of overlapping sequences between any two resultant subsets.
        binary  ( Optional[string] )
            Name or path of the DCM binary file.

    Returns :
        list
            List of lists, one per set, composed by each set's corresponding
            sequence ids.

    Raises :
        RuntimeError
            If the call to the dcm command raises an exception.
        IOError
            If the dcm tool can't generate a decomposition for the 'subset_size'
            and 'overlapping' values given.

    * The tree file must be in NEWICK format.
    """
    # Build the command line and launch it handling any Runtime exception
    command = [binary, '3', tree_file, overlapping]
    try :
        output = subprocess.check_output(command, stderr=DEVNULL,
                                         universal_newlines=True)
    except subprocess.CalledProcessError as e :
        message = 'Running "{}" raised an exception'.format(' '.join(e.cmd))
        raise RuntimeError(message)
    else :
        # Split the output in lines, removing the first and last lines which are
        # not part of the solution
        sets_decomp = output.split('\n')[1:-1]
        if ( len(sets_decomp) == 1 ) :
            raise IOError('DCM can\'t obtain a decomposition with the given ' \
                          '"subset_size" and "overlapping" values')
        # else : # len(output_subsets) > 1
        set_list = []
        further_decomp = []
        for set_ids_list in sets_decomp :
            # Get the sequence identifiers, removing the text format elements
            seq_id_list = set_ids_list.split(' ')[1:-1]
            if ( len(seq_id_list) <= subset_size ) :
                # The set is below the threshold, so there is no need for
                # further decomposition
                set_list.append(seq_id_list)
            else : # len(subset) > self._subset
                # Extract the corresponding subtree of this set and save it in a
                # temporal file for a further recursive decomposition
                subtree = dendropy.Tree.get(path=tree_file, schema='newick')
                subtree.retain_taxa_with_labels(seq_id_list)
                tmpfile = tempfile.NamedTemporaryFile(delete=False)
                subtree.write(path=tmpfile.name, schema='newick')
                further_decomp.append(tmpfile.name)
        return ( set_list, further_decomp )


#-------------------------------------------------------------------------------

def map_seqs ( record_list, tree_file, file_format, subset_size, overlapping,
               binary = 'dcm' ) :
    """
    Generate a map of the sequences in sets, of at most 'subset_size', with the
    specified overlapping using the padded-Recursive-DMC3 decomposition (PRD)
    from DACTAL system. If 'tree_file' contains a relative path, the current
    working directory will be used to get the absolute path.

    Arguments :
        record_list  ( list )
            List of SeqRecord objects (from Biopython).
        tree_file  ( string )
            Input tree file.
        file_format  ( string )
            Tree file format.
        subset_size  ( int )
            Maximum subset size.
        overlapping  ( int )
            Number of overlapping sequences between any two resultant subsets.
        binary  ( Optional[string] )
            Name or path of the DCM binary file.

    Returns :
        dict
            Dictionary with the set identifiers as keys and the corresponding
            sequences as values in lists of SeqRecord objects.

    Raises :
        ValueError
            When 'subset_size' < (4 * 'overlapping').
        RuntimeError
            If the call to the dcm command raises an exception.
        IOError
            If the dcm tool can't generate a decomposition for the 'subset_size'
            and 'overlapping' values given.

    * The tree file format must be supported by Bio.Phylo.
    """
    if ( subset_size < (4 * overlapping) ) :
        raise ValueError('The maximum subset size must be greater than or ' \
                         'equal to 4 times the overlapping value')
    # else : # subset_size >= (4 * overlapping)
    # If the input file format is not supported by the PRD process, convert it
    # to a temporary supported file
    infile_path = Path(tree_file).resolve()
    if ( file_format.lower() != 'newick' ) :
        tmpfile = tempfile.NamedTemporaryFile()
        Phylo.convert(infile_path, file_format, tmpfile.name, 'newick')
        infile_path = tmpfile.name
    # The first decomposition process will be always executed, so there is no
    # need to overload this stage with the multiprocess generation
    set_list, further_decomp = _prd_decomposition(infile_path, subset_size,
                                                  str(overlapping), binary)
    # Parallelization of the recursive decomposition of the different subtrees.
    # All new subtrees are attached to 'further_decomp' file list so we can
    # launch at most one process per core, speeding up the whole process
    start = 0
    to_process = len(further_decomp[start:])
    pool = multiprocessing.Pool(processes=NUMCORES)
    while ( to_process > 0 ) :
        end = start + min(to_process, NUMCORES)
        results = [pool.apply_async(_prd_decomposition,
                                    args=(further_decomp[i], subset_size,
                                          str(overlapping), binary,))
                           for i in range(start, end)]
        # Collect the results of all the processes launched
        for pool_result in results :
            output = pool_result.get()
            set_list += output[0]
            further_decomp += output[1]
        start = end
        to_process = len(further_decomp[start:])
    # Remove all the temporal files created for the multirpocessing stage
    for file_path in further_decomp :
        os.remove(file_path)
    record_dict = {record.id: record  for record in record_list}
    # Map all the resultant sets with an unique set id and replace the sequence
    # ids by their corresponding Bio.SeqRecord object
    set_dict = {}
    num_zeros = len(str(len(set_list)))
    for index, seq_id_list in enumerate(set_list, 1) :
        set_id = 'prdset{}'.format(str(index).zfill(num_zeros))
        set_dict[set_id] = []
        for seq_id in seq_id_list :
            set_dict[set_id].append(record_dict[seq_id])
    return ( set_dict )


#-------------------------------------------------------------------------------
