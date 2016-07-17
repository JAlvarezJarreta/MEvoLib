#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  NaiveRows.py
# Last version :  v1.00 ( 01/Dec/2015 )
# Description :  Clustering where each resulting set is composed by a specific
#       range of input sequences (row/sequence division).
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
import random


#-------------------------------------------------------------------------------

def map_seqs ( record_list, num_sets ) :
    """
    Naive distribution in 'num_sets' sets of the sequences at 'record_list'. The
    maximum number of sequences per set is calculated as follows:

        sequences_per_set = ceiling( total_number_of_sequences / 'num_sets' )

    Arguments :
        record_list  ( list )
            List of SeqRecord objects (from Biopython).
        num_sets  ( int )
            Number of sets.

    Returns :
        dict
            Dictionary with the set identifiers as keys and the corresponding
            sequences as values in lists of SeqRecord objects.
    """
    num_seqs = len(record_list)
    # Determine the minimum and maximum number of sequences per set
    min_seqs_set = int(math.floor(num_seqs / num_sets))
    max_seqs_set = min_seqs_set + 1
    # Get the number of sets that will have the maximum number of sequences to
    # balance the distribution of sequences per set
    big_sets = num_seqs % num_sets
    small_sets = num_sets - big_sets
    # Get a random distribution of the set size list
    set_size_list = [min_seqs_set] * small_sets + [max_seqs_set] * big_sets
    random.shuffle(set_size_list)
    # Minimum string length for the given number of sets (with zero-filling)
    num_zeros = len(str(num_sets))
    set_dict = {}
    start = 0
    for index, set_size in enumerate(set_size_list, 1) :
        set_id = 'rset{}'.format(str(index).zfill(num_zeros))
        end = start + set_size
        set_dict[set_id] = record_list[start:end]
        start = end
    return ( set_dict )
    

#-------------------------------------------------------------------------------
