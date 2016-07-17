#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  NaiveCols.py
# Last version :  v1.00 ( 01/Dec/2015 )
# Description :  Clustering where each resulting set is composed by a specific
#       range slice of all the input sequences (column/fragment division).
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  01/Dec/2015
#   VERSION :  v1.00
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

from __future__ import absolute_import, division

import math

from Bio.SeqRecord import SeqRecord


#-------------------------------------------------------------------------------

def map_seqs ( record_list, num_sets ) :
    """
    Naive splicing in 'num_sets' sets of the sequences at 'record_list'. The
    maximum fragment length per set is calculated as follows:

        frag_length_per_set = ceiling( maximum_sequence_length / 'num_sets' )

    Arguments :
        record_list  ( list )
            List of SeqRecord objects (from Biopython).
        seq_length  ( int )
            Number of sets of sequence slices.

    Returns :
        dict
            Dictionary with the set identifiers as keys and the corresponding
            sequence fragments as values in lists of SeqRecord objects.
    """
    # Get maximum length to determine fragment size for the given number of sets
    max_len = len(max(record_list, key=lambda x: len(x.seq)))
    max_frag_len = int(math.ceil(max_len / num_sets))
    # Minimum string length for the given number of sets (with zero-filling)
    num_zeros = len(str(num_sets))
    # Generate all the set ids and their corresponding starting site
    frag_list = [('cset{}'.format(str(i).zfill(num_zeros)), value)
                    for i, value in enumerate(range(0, max_len, max_frag_len),
                                              1)]
    set_dict = {}
    for record in iter(record_list) :
        for set_id, start in frag_list :
            if ( len(record) < start ) :
                # The current sequence can't be divided into more sets
                break
            else : # len(record) >= start
                end = start + max_frag_len
                frag_record = SeqRecord(record.seq[start:end], id=record.id,
                                        name=record.name, description=set_id)
                set_dict.setdefault(set_id, []).append(frag_record)
    return ( set_dict )


#-------------------------------------------------------------------------------
