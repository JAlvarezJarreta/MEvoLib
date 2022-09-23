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
"""Clustering where each resulting set is composed by a specific range slice of all the input sequences 
(column/fragment division)."""

import math

from typing import List

from Bio.SeqRecord import SeqRecord


def map_seqs (record_list: List, num_sets: int) -> dict:
    """Naive splicing in 'num_sets' sets of the sequences at 'record_list'. The maximum fragment length per 
    set is calculated as follows:

        frag_length_per_set = ceiling( maximum_sequence_length / 'num_sets' )

    Args:
        record_list: List of SeqRecord objects (from Biopython).
        seq_length: Number of sets of sequence slices.

    Returns:
        dict: Dictionary with the set identifiers as keys and the corresponding sequence fragments as values 
        in lists of SeqRecord objects.
    """
    # Get maximum and minimum length to determine fragment size for the given number of sets
    max_len = len(max(record_list, key = lambda x: len(x.seq)))
    max_frag_len = int(math.ceil(float(max_len) / num_sets))
    min_frag_len = max_frag_len - 1
    # Get the number of sets that will have the maximum length sequences to balance the distribution per set
    big_sets = max_len % num_sets
    small_sets = num_sets - big_sets	
    # Minimum string length for the given number of sets (with zero-filling)
    num_zeros = len(str(num_sets))
    # Generate all the set ids and their corresponding starting site
    first_range = big_sets * max_frag_len
    frag_list = [('cset{}'.format(str(i).zfill(num_zeros)), value)
                    for i, value in enumerate(range(0, first_range,
					                                max_frag_len), 1)]
    frag_list += [('cset{}'.format(str(i).zfill(num_zeros)), value)
                    for i, value in enumerate(range(first_range, max_len,
					                                min_frag_len), big_sets+1)]
    set_dict = {}
    for record in iter(record_list):
        for set_id, start in frag_list:
            if len(record) < start:
                # The current sequence can't be divided into more sets
                break
            else: # len(record) >= start
                end = start + max_frag_len
                frag_record = SeqRecord(record.seq[start:end], id=record.id,
                                        name=record.name, description=set_id)
                set_dict.setdefault(set_id, []).append(frag_record)
    return set_dict