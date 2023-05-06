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
"""Functions aimed to provide an easy interface to handle"""

from __future__ import annotations

import argparse
import multiprocessing
import math

from Bio import SeqIO

from mevolib.cluster import Genes, NaiveRows, NaiveCols, PRD

from mevolib._utils import NUMCORES

from pathlib import Path

from mevolib._utils import get_abspath
from typing import Tuple

from time import sleep


_METHOD_TO_FUNC = { 'genes': Genes.map_seqs,
                    'rows': NaiveRows.map_seqs,
                    'cols': NaiveCols.map_seqs,
                    'prd': PRD.map_seqs }


def get_tools() -> list:
    """Returns a list of clustering methods and software tools included in the current version of MEvoLib."""
    return list(_METHOD_TO_FUNC.keys())


def get_subsets(method: str, seqfile: str, fileformat: str = 'genbank', *args: tuple, **kwargs: dict) -> dict:
    """Division of all the sequences stored in the input file into subsets applying the `method` function.
    
    If `seqfile` contains a relative path, the current working directory will be used to get the absolute path.

    Args:
        method: Desired partition method (case-insensitive): genes, naive rows or cols, padded-Recursive-DCM3.
        seqfile: Input sequences file.
        fileformat: Input file format.
        args & kwargs: Non-keyworded and keyworded arguments passed to the selected method.
          
    Raises:
        ValueError: If there is no corresponding method to `method` value.
        IOError: If the path or the file provided doesn't exist.
        IOError: If the file format provided doesn't correspond to the actual one.
    
    Returns: 
        A dict with the set identifiers as keys and the corresponding sequences as values in lists of 
            SeqRecord objects.


    * The file format must be supported by Bio.SeqIO.
    * For "rows" method, if the number of input sequences is lower than the number of sets multiplied by the 
    number of cores, the resulting sets might be fewer than the number requested.
    """
    method_key = method.lower()
    if method_key not in _METHOD_TO_FUNC:
        raise ValueError(f'The method "{method}" is not included in mevolib.cluster')
    # else: # method_key in _METHOD_TO_FUNC
    # Get the mapping function and the sequence file path
    mapseqs_func = _METHOD_TO_FUNC[method_key]
    filepath = Path(seqfile).resolve()
    if method_key in ['prd', 'genes']:
        # Non data-driven (throught input slicing) parallelizable methods
        seq_list = (x  for x in SeqIO.parse(filepath, fileformat))
        set_dict = mapseqs_func(seq_list, *args, **kwargs)
    else:
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
        # Build the final sets dictionary merging the results of all executed processes
        output = [p.get() for p in iter(results)]
        set_dict = output[0]
        for key in iter(set_dict):
            for result in output[1:]:
                set_dict[key].extend(result[key])
    return set_dict

def parallel_seqio_write(args: Tuple) -> None:
    # Unpack DNA sequence and destination FASTA filename
    seq, filename = args
    SeqIO.write(seq, filename, 'fasta')


def main():
    """Default call for Genes module."""     
    parser = argparse.ArgumentParser(
        description="Performs the division of every sequence stored in the input file into subsets by genes"
    )
    parser.add_argument("-i", "--input", required=True, help="Input sequences file")
    parser.add_argument("-f", "--format", default="genbank", help="Input file format")
    parser.add_argument("-o", "--output", required=True, help="Output file name (without extension)")
    args = parser.parse_args()
    gene_dict = get_subsets('genes', args.input, args.format, log_file=args.output + '.log') 
    
    #  We dump the split sequences stored in the dictionary into a fasta file.
    args_iterator = ((value, f"{key}.fasta") for key, value in gene_dict.items() if len(value)>0)
    pool = multiprocessing.Pool(processes=NUMCORES)
    pool.imap_unordered(parallel_seqio_write, args_iterator)
    sleep(5)


