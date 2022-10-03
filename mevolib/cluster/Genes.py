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
"""Clustering where each resulting set is composed by a gene of all the input sequences 
(from GenBank data or reference sequence)."""

import tempfile
import itertools
import math
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from mevolib import align
from mevolib.data import rCRS


_REF_SEQ_DICT = {'rCRS': rCRS}

# Dictionary with the main qualifiers associated with all possible features of each sequence record at GenBank
_FEAT_QUAL_DICT = { 'assembly_gap': ['gap_type'],
                    'C_region': ['gene', 'standard_name'],
                    'CDS': ['gene', 'product', 'standard_name'],
                    'centromere': ['standard_name'],
                    'D-loop': ['gene'],
                    'D_segment': ['gene', 'product', 'standard_name'],
                    'exon': ['gene', 'product', 'standard_name'],
                    'gap': ['note'],
                    'gene': ['gene', 'product', 'standard_name'],
                    'iDNA': ['gene', 'standard_name'],
                    'intron': ['gene', 'standard_name'],
                    'J_segment': ['gene', 'product', 'standard_name'],
                    'LTR': ['gene', 'standard_name'],
                    'mat_peptide': ['gene', 'product', 'standard_name'],
                    'misc_binding': ['gene'],
                    'misc_difference': ['gene', 'standard_name'],
                    'misc_feature': ['gene', 'product', 'standard_name'],
                    'misc_recomb': ['gene', 'standard_name'],
                    'misc_RNA': ['gene', 'product', 'standard_name'],
                    'misc_structure': ['gene', 'standard_name'],
                    'mobile_element': ['gene', 'standard_name'],
                    'modified_base': ['gene'],
                    'mRNA': ['product', 'gene', 'standard_name'],
                    'ncRNA': ['product', 'gene', 'standard_name'],
                    'N_region': ['gene', 'product', 'standard_name'],
                    'old_sequence': ['gene'],
                    'operon': ['standard_name'],
                    'oriT': ['gene', 'standard_name'],
                    'polyA_site': ['gene'],
                    'precursor_RNA': ['gene', 'product', 'standard_name'],
                    'prim_transcript': ['gene', 'standard_name'],
                    'primer_bind': ['gene', 'standard_name'],
                    'protein_bind': ['gene', 'standard_name'],
                    'regulatory': ['gene', 'standard_name'],
                    'repeat_region': ['gene', 'standard_name'],
                    'rep_origin': ['gene', 'standard_name'],
                    'rRNA': ['product', 'gene', 'standard_name'],
                    'S_region': ['gene', 'product', 'standard_name'],
                    'sig_peptide': ['gene', 'product', 'standard_name'],
                    #'source': ['organism', 'mol_type'],
                    'stem_loop': ['gene', 'standard_name'],
                    'STS': ['gene', 'standard_name'],
                    'telomere': ['standard_name'],
                    'tmRNA': ['product', 'gene', 'standard_name'],
                    'transit_peptide': ['gene', 'product', 'standard_name'],
                    'tRNA': ['product', 'gene', 'standard_name'],
                    'unsure': ['gene'],
                    'V_region': ['gene', 'product', 'standard_name'],
                    'V_segment': ['gene', 'product', 'standard_name'],
                    'variation': ['gene', 'product', 'standard_name'],
                    '3\'UTR': ['gene', 'standard_name'],
                    '5\'UTR': ['gene', 'standard_name'] }


def _normalization(record: Bio.SeqRecord, refseq_record: Bio.SeqRecord, 
                   alignment_bin: str) -> tuple(Bio.Seq, list):
    """Normalizes the input sequence with the reference sequence.
    
    The normalization consists on aligning both sequences and removing those sites where a gap has
    been introduced in the reference sequence, making the features of the reference sequence
    applicable to the new sequence. It returns the new sequence and the reference sequence's features.

    Args:
        record: Sequence to normalize.
        refseq_record: Reference sequence of the same type as `record`.
        alignment_bin: Binary path to the alignment tool that will be used in the normalization process.

    Returns:
        Bio.Seq: Normalized sequence of `record`.
        list: List of SeqFeature objects (from Biopython) from reference sequence.

    Raises:
        RuntimeError: If the call to the alignment tool command raises an exception.
    """
    tmpfile = tempfile.NamedTemporaryFile()
    SeqIO.write([refseq_record, record], tmpfile.name, 'fasta')
    alignment = align.get_alignment(alignment_bin, tmpfile.name, 'fasta')
    # Get the normalized sequence by removing the sites that correspond to gaps introduced in the reference 
    # sequence during the alignment process
    record_seq = ''.join((x  for i, x in enumerate(alignment[1])
                                 if alignment[0][i] != '-'))
    return (Seq(record_seq, refseq_record.seq.alphabet), refseq_record.features)


def _string_filter(input_str: str) -> str:
    """Returns a lower case string of `input_str` with known GenBank's misspellings corrected.

    Args:
        input_str: Input string to filter.

    """
    # Ignore case
    final_str = input_str.lower()
    # Correct misspellings
    final_str = final_str.replace('oxydase', 'oxidase')
    return final_str


def get_features() -> List:
    """Returns a list of all possible feature keywords that can be found in any GenBank's sequence record."""
    return list(_FEAT_QUAL_DICT.keys())


def map_seqs (record_list: list, feature_filter: Optional[list] = None, ref_seq: Optional[str] = None, 
              alignment_bin: Optional[str] = None, log_file: Optional[str] = None) -> dict:
    """Gene splicing of the sequences at `record_list`.

    By default, the gene location is extracted from the feature list of each sequence. If there is
    no list, that sequence is classified as "unprocessable" or, if a reference sequence is given,
    the reference features are used to extract the different genes (through a normalization process
    using an alignment tool). All the features are returned unless a list of feature keywords are
    passed through `feature_filter`. If a log file path is given and any file exists with that name,
    the file will be overwritten without any warning.
    with that name, the file will be overwritten without any warning.

    Args:
        record_list: List of SeqRecord objects (from Biopython).
        feature_filter: List of feature keywords the user wants to be returned (from all the possible ones).
        ref_seq: Keyword (from MEvoLib.Data) or file path (GENBANK format) of the reference sequence.
        alignment_bin: Binary path of the alignment tool (only required if a reference sequence is passed).
        log_file: Absolute path for the log file.

    Returns:
        dict: Dictionary with the set identifiers as keys and the corresponding sequence fragments as values 
            in lists of SeqRecord objects.

    Raises:
        IOError: If the reference sequence's file path doesn't exist.
        RuntimeError: If the call to the alignment tool command raises an exception.

    * Reference sequence's file must be in GENBANK format.
    """
    # Load the desired feature keywords as keys of the gene dictionary and a term dictionary with a list of 
    # sequences for each qualifier of any selected feature
    if feature_filter:
        gene_dict = dict((key, {}) for key in feature_filter)
        term_dict = dict((key, {}) for key in feature_filter)
    else:
        gene_dict = dict((key, {}) for key in _FEAT_QUAL_DICT.keys())
        term_dict = dict((key, {}) for key in _FEAT_QUAL_DICT.keys())
    # Get the reference sequence's SeqRecord object or create an unprocessable list for those sequences 
    # without gene information
    if ref_seq in _REF_SEQ_DICT:
        refseq_record = _REF_SEQ_DICT[ref_seq].RECORD
    elif ref_seq:
        refseq_record = SeqIO.read(ref_seq, 'gb')
    else:
        unprocessable = []
    num_seqs = 0
    # Iterate over all the records to get their gene division
    for record in record_list:
        num_seqs += 1
        if len(record.features) <= 1:
            # GenBank's "source" feature key is mandatory
            if ref_seq:
                record.seq, record.features = _normalization(record,
                                                             refseq_record,
                                                             alignment_bin)
            else:
                unprocessable.append(record)
                continue
        # else: # len(record.features) > 1
        record_features = (feat  for feat in record.features[1:]
                                     if feat.type in gene_dict)
        for feature in record_features:
            # Create a set of qualifiers of the record from the main fields of GenBank 
            # (pre-saved in _FEAT_QUAL_DICT)
            record_qualifiers = set()
            for qualifier_key in iter(_FEAT_QUAL_DICT[feature.type]):
                if qualifier_key in feature.qualifiers:
                    record_qualifiers.update((_string_filter(x)  for x in
                                             feature.qualifiers[qualifier_key]))
            if not record_qualifiers:
                # 'record_qualifiers' is empty
                record_qualifiers.add(feature.type)
            # Generate a string of the qualifiers' set to store it as a description of the gene SeqRecord 
            # object
            qualifier_id = ':'.join(sorted(record_qualifiers,
                                           key=lambda item: (len(item), item)))
            feature_record = SeqRecord(feature.extract(record).seq,
                                       id=record.id, name=record.id,
                                       description=qualifier_id)
            # Add new terms to the corresponding entry of the dictionary for the given feature, or add the 
            # sequence record id to the existing entry
            for pair in itertools.combinations(qualifier_id.split(':'), 2):
                if pair not in term_dict[feature.type]:
                    term_dict[feature.type][pair] = set([record.id])
                else: # pair in term_dict[feature.type]
                    term_dict[feature.type][pair].add(record.id)
            # Merge possible matching qualifiers for the same type of feature
            qualifiers_to_merge = []
            for key in gene_dict[feature.type].keys():
                key_set = set(key.split(':'))
                if not record_qualifiers.isdisjoint(key_set):
                    if record_qualifiers <= key_set:
                        record_qualifiers.update(key_set)
                    elif record_qualifiers > key_set:
                        qualifiers_to_merge.append(key)
                    else:
                        # 'record_qualifiers' and 'key_set' differ but their intersection is not empty
                        record_qualifiers.update(key_set)
                        qualifiers_to_merge.append(key)
            # Generate new qualifier string
            qualifier_id = ':'.join(sorted(record_qualifiers,
                                           key=lambda item: (len(item), item)))
            # Add the new gene SeqRecord object to the dictionary
            if qualifier_id not in gene_dict[feature.type]:
                gene_dict[feature.type][qualifier_id] = [feature_record]
            else: # qualifier_id in gene_dict[feature.type]
                gene_dict[feature.type][qualifier_id].append(feature_record)
            # Merge those qualifiers that belong to the same gene
            for qualifier_key in qualifiers_to_merge:
                if qualifier_key != qualifier_id:
                    gene_dict[feature.type][qualifier_id].extend(
                        gene_dict[feature.type][qualifier_key])
                    del gene_dict[feature.type][qualifier_key]
    # The error calculation has been extracted from the following sampling statistics equation:
    #
    #                   N * Z^2 * p * (1-p)
    #         n = -------------------------------
    #              (N-1) * e^2 + Z^2 * p * (1-p)
    #
    # where N is the number of sequences, n is the minimum sampling size (threshold), e is the error fixed to
    # 0,01, Z is fixed to get a 0,99 confidence interval and p is assumed to be 0,5.
    e_value = 0.01
    z_value = 2.58
    p_value = 0.5
    coef = math.pow(z_value, 2) * p_value * (1.0 - p_value)
    threshold = math.ceil((num_seqs * coef) / \
                          ((num_seqs - 1.0) * math.pow(e_value, 2) + coef))
    # If no log file path is provided, save log content in a named temporary
    # file that won't be deleted after the function ends
    if not log_file:
        log_file = (tempfile.NamedTemporaryFile(delete=False)).name
    # Clean-up empty features and merge qualifiers dict keys with features dict keys to get a {str: list} 
    # dict for all the genes
    set_dict = {}
    with open(log_file, 'w') as log:
        for feat_key, feat_value in gene_dict.items():
            if feat_value:
                log.write('> {}\n'.format(feat_key))
                for qual_key, qual_value in iter(feat_value.items()):
                    # Generation of the content of the set dictionary that will
                    # be returned
                    new_key = '{}.{}'.format(feat_key, qual_key.split(':')[0])
                    set_dict.setdefault(new_key, []).extend(qual_value)
                    # For every existing pair of qualifiers, if the number of 
                    # records that hold both is below the calculated threshold,
                    # it might be the result of a typo in those records'
                    # information (further review of the log file is advisable)
                    for pair in itertools.combinations(qual_key.split(':'), 2):
                        if pair in term_dict[feat_key]:
                            sampling_size = len(term_dict[feat_key][pair])
                            if sampling_size < threshold:
                                seq_list = list(term_dict[feat_key][pair])
                                text = '\t{}\n'.format(' || '.join(pair))
                                for i in range(0, sampling_size // 6 + 1):
                                    text += '\t\t{}\n'.format(
                                                ' '.join(seq_list[i*6:(i+1)*6]))
                                if (sampling_size % 6) != 0:
                                    text += '\n'
                                log.write(text)
                log.write('\n')
    # If no reference sequence has been introduced, include in the gene dict those sequences that couldn't be 
    # processed due to lack of information
    if not ref_seq:
        set_dict['unprocessable'] = unprocessable
    return set_dict
