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
"""Clustering where each resulting set is composed by a specific range of input sequences 
(row/sequence division)."""

import math
import random



def map_seqs (record_list: list, num_sets: int) -> dict:
    """Naive distribution in `num_sets` sets of the sequences at `record_list`.
    
    The maximum number of sequences per set is calculated as follows:

        sequences_per_set = ceiling( total_number_of_sequences / 'num_sets' )

    Args:
        record_list: List of SeqRecord objects (from Biopython).
        num_sets: Number of sets.

    Returns:
        dict: Dictionary with the set identifiers as keys and the corresponding sequences as values in lists 
        of SeqRecord objects.
    """
    num_seqs = len(record_list)
    # Determine the minimum and maximum number of sequences per set
    min_seqs_set = int(math.floor(float(num_seqs) / num_sets))
    max_seqs_set = min_seqs_set + 1
    # Get the number of sets that will have the maximum number of sequences to balance the distribution of 
    # sequences per set
    big_sets = num_seqs % num_sets
    small_sets = num_sets - big_sets
    # Get a random distribution of the set size list
    set_size_list = [min_seqs_set] * small_sets + [max_seqs_set] * big_sets
    random.shuffle(set_size_list)
    # Minimum string length for the given number of sets (with zero-filling)
    num_zeros = len(str(num_sets))
    set_dict = {}
    start = 0
    for index, set_size in enumerate(set_size_list, 1):
        set_id = 'rset{}'.format(str(index).zfill(num_zeros))
        end = start + set_size
        set_dict[set_id] = record_list[start:end]
        start = end
    return set_dict