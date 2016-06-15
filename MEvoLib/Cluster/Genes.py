#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  Genes.py
# Last version :  v1.2 ( 15/Jun/2016 )
# Description :  Clustering where each resulting set is composed by a gene of
#       all the input sequences (from GenBank data or reference sequence).
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  15/Jun/2016
#   VERSION :  v1.2
#   AUTHOR(s) :  J. Alvarez-Jarreta
#   CHANGES :  * Fixed two bugs in map_seqs() method: i) an empty dictionary was
#                  returned if the input file was not in GENBANK format;
#                  ii) the error calculation was raising an exception if the
#                  number of sequences was lower than the sampling size.
#
#   DATE :  11/May/2016
#   VERSION :  v1.1
#   AUTHOR(s) :  J. Alvarez-Jarreta
#   CHANGES :  * The map_seqs() method now generates a report of possible wrong
#                  pairs of terms (maybe a typo in GenBank information?) based
#                  on sampling statistics applied to the input sequences.
#              * Fixed a bug in map_seqs() method with a possible random
#                  ordering of terms for the same gene making it not unique in
#                  the generation of the qualifier id.
#
#   DATE :  07/Feb/2016
#   VERSION :  v1.0
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

from __future__ import absolute_import

import tempfile
import itertools
import math
import warnings

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from MEvoLib import Align
from MEvoLib.Data import rCRS
from MEvoLib._py3k import viewkeys, viewitems


#-------------------------------------------------------------------------------

_REF_SEQ_DICT = { 'rCRS': rCRS }

# Dictionary with the main qualifiers associated with all possible features of
# each sequence record at GenBank
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


#-------------------------------------------------------------------------------

def _normalization ( record, refseq_record, alignment_bin ) :
    """
    Normalization of the input sequence with the reference sequence. The
    normalization consists on aligning both sequences and removing those sites
    where a gap has been introduced in the reference sequence, making the
    features of the reference sequence applicable to the new sequence. It
    returns the new sequence and the reference sequence's features.
    Arguments :
        record  ( Bio.SeqRecord )
            Sequence to normalize.
        refseq_record  ( Bio.SeqRecord )
            Reference sequence of the same type as 'record'.
        alignment_bin  ( string )
            Binary path to the alignment tool that will be used in the
            normalization process.
    Returns :
        Bio.Seq
            Normalized sequence of 'record'.
        list
            List of SeqFeature objects (from Biopython) from reference sequence.
    Raises :
        RuntimeError
            If the call to the alignment tool command raises an exception.
    """
    tmpfile = tempfile.NamedTemporaryFile()
    SeqIO.write([refseq_record, record], tmpfile.name, 'fasta')
    alignment = Align.get_alignment(alignment_bin, tmpfile.name, 'fasta')
    # Get the normalized sequence by removing the sites that correspond to gaps
    # introduced in the reference sequence during the alignment process
    record_seq = ''.join((x  for i, x in enumerate(alignment[1])
                                 if alignment[0][i] != '-'))
    return ( Seq(record_seq, refseq_record.seq.alphabet),
             refseq_record.features )



def _string_filter ( input_str ) :
    """
    Arguments :
        input_str  ( string )
            Input string to filter.
    Returns :
        string
            Lower case of 'input_str' with known GenBank's misspellings
            corrected.
    """
    # Ignore case
    final_str = input_str.lower()
    # Correct misspellings
    final_str = final_str.replace('oxydase', 'oxidase')
    return ( final_str )


#-------------------------------------------------------------------------------

def get_features ( ) :
    """
    Returns :
        list
            List of all possible feature keywords that can be found in any
            GenBank's sequence record.
    """
    return ( [x  for x in iter(viewkeys(_FEAT_QUAL_DICT))] )



def map_seqs ( record_list, feature_filter = None, ref_seq = None,
               alignment_bin = None, log_file = None ) :
    """
    Gene splicing of the sequences at 'record_list'. By default, the gene
    location is extracted from the feature list of each sequence. If there is no
    list, that sequence is classified as "unprocessable" or, if a reference
    sequence is given, the reference features are used to extract the different
    genes (through a normalization process using an alignment tool). All the
    features are returned unless a list of feature keywords are passed through
    'feature_filter' parameter. If a log file path is given and any file exists
    with that name, the file will be overwritten without any warning.
    Arguments :
        record_list  ( list )
            List of SeqRecord objects (from Biopython).
        feature_filter  ( Optional[list] )
            List of feature keywords the user wants to be returned (from all the
            possible ones).
        ref_seq  ( Optional[string] )
            Keyword (from MEvoLib.Data) or file path (GENBANK format) of the
            reference sequence.
        alignment_bin  ( Optional[string] )
            Binary path of the alignment tool (only required if a reference
            sequence is passed).
        log_file  ( Optional[string] )
            Absolute path for the log file.
    Returns :
        dict
            Dictionary with the set identifiers as keys and the corresponding
            sequence fragments as values in lists of SeqRecord objects.
    Raises :
        IOError
            If the reference sequence's file path doesn't exist.
        RuntimeError
            If the call to the alignment tool command raises an exception.
    * Reference sequence's file must be in GENBANK format.
    """
    # Load the desired feature keywords as keys of the gene dictionary and a
    # term dictionary with a list of sequences for each qualifier of any
    # selected feature
    if ( feature_filter ) :
        gene_dict = dict((key, {})  for key in iter(feature_filter))
        term_dict = dict((key, {})  for key in iter(feature_filter))
    else : # feature_filter is None
        gene_dict = dict((key, {})  for key in iter(viewkeys(_FEAT_QUAL_DICT)))
        term_dict = dict((key, {})  for key in iter(viewkeys(_FEAT_QUAL_DICT)))
    # Get the reference sequence's SeqRecord object or create an unprocessable
    # list for those sequences without gene information
    if ( ref_seq in _REF_SEQ_DICT ) :
        refseq_record = _REF_SEQ_DICT[ref_seq].RECORD
    elif ( ref_seq ) : # ref_seq != None
        refseq_record = SeqIO.read(ref_seq, 'gb')
    else : # ref_seq is None
        unprocessable = []
    num_seqs = 0
    # Iterate over all the records to get their gene division
    for record in iter(record_list) :
        num_seqs += 1
        if ( len(record.features) <= 1 ) :
            # GenBank's "source" feature key is mandatory
            if ( ref_seq ) :
                record.seq, record.features = _normalization(record,
                                                             refseq_record,
                                                             alignment_bin)
            else : # ref_seq is None
                unprocessable.append(record)
                continue
        # else : # len(record.features) > 1
        record_features = (feat  for feat in record.features[1:]
                                     if feat.type in gene_dict)
        for feature in record_features :
            # Create a set of qualifiers of the record from the main fields of
            # GenBank (pre-saved in _FEAT_QUAL_DICT)
            record_qualifiers = set()
            for qualifier_key in iter(_FEAT_QUAL_DICT[feature.type]) :
                if ( qualifier_key in feature.qualifiers ) :
                    record_qualifiers.update((_string_filter(x)  for x in
                                             feature.qualifiers[qualifier_key]))
            if ( not record_qualifiers ) :
                # 'record_qualifiers' is empty
                record_qualifiers.add(feature.type)
            # Generate a string of the qualifiers' set to store it as a
            # description of the gene SeqRecord object
            qualifier_id = ':'.join(sorted(record_qualifiers,
                                           key=lambda item: (len(item), item)))
            feature_record = SeqRecord(feature.extract(record).seq,
                                       id=record.id, name=record.id,
                                       description=qualifier_id)
            # Add new terms to the corresponding entry of the dictionary for
            # the given feature, or add the sequence record id to the existing
            # entry
            for pair in itertools.combinations(qualifier_id.split(':'), 2) :
                if ( pair not in term_dict[feature.type] ) :
                    term_dict[feature.type][pair] = set([record.id])
                else : # pair in term_dict[feature.type]
                    term_dict[feature.type][pair].add(record.id)
            # Merge possible matching qualifiers for the same type of feature
            qualifiers_to_merge = []
            for key in iter(viewkeys(gene_dict[feature.type])) :
                key_set = set(key.split(':'))
                if ( not record_qualifiers.isdisjoint(key_set) ) :
                    if ( record_qualifiers <= key_set ) :
                        record_qualifiers.update(key_set)
                    elif ( record_qualifiers > key_set ) :
                        qualifiers_to_merge.append(key)
                    else :
                        # 'record_qualifiers' and 'key_set' differ but their
                        # intersection is not empty
                        record_qualifiers.update(key_set)
                        qualifiers_to_merge.append(key)
            # Generate new qualifier string
            qualifier_id = ':'.join(sorted(record_qualifiers,
                                           key=lambda item: (len(item), item)))
            # Add the new gene SeqRecord object to the dictionary
            if ( qualifier_id not in gene_dict[feature.type] ) :
                gene_dict[feature.type][qualifier_id] = [feature_record]
            else : # qualifier_id in gene_dict[feature.type]
                gene_dict[feature.type][qualifier_id].append(feature_record)
            # Merge those qualifiers that belong to the same gene
            for qualifier_key in iter(qualifiers_to_merge) :
                gene_dict[feature.type][qualifier_id].extend(
                                         gene_dict[feature.type][qualifier_key])
                del gene_dict[feature.type][qualifier_key]
    # The error calculation has been extracted from the following sampling
    # statistics equation:
    #
    #                   N * Z^2 * p * (1-p)
    #         n = -------------------------------
    #              (N-1) * e^2 + Z^2 * p * (1-p)
    #
    # where N is the number of sequences, n is the sampling size, e is the
    # error, Z is fixed to get a 0,95 confidence interval and p is assumed
    # to be 0,5. The pre-calculable part has been saved at 'coef' variable.
    z_value = 1.96
    p_value = 0.5
    coef = math.sqrt((math.pow(z_value, 2) * p_value * (1.0 - p_value)) / \
                     (num_seqs - 1.0))
    # If no log file path is provided, save log content in a named temporary
    # file that won't be deleted after the function ends
    if ( not log_file ) :
        log_file = (tempfile.NamedTemporaryFile(delete=False)).name
    # Clean-up empty features and merge qualifiers dict keys with features dict
    # keys to get a {str: list} dict for all the genes
    set_dict = {}
    with open(log_file, 'w') as log :
        for feat_key, feat_value in iter(viewitems(gene_dict)) :
            if ( feat_value ) :
                log.write('> {}\n'.format(feat_key))
                for qual_key, qual_value in iter(viewitems(feat_value)) :
                    # Generation of the content of the set dictionary that will
                    # be returned
                    new_key = '{}.{}'.format(feat_key, qual_key.split(':')[0])
                    set_dict.setdefault(new_key, []).extend(qual_value)
                    # Calculate the error from the established values for a
                    # correct sampling statistics application on the dataset
                    sampling_size = float(len(qual_value))
                    # If the number of sequences is lower than the sampling
                    # size, we cannot calculate the square root of a negative
                    # number so we skip these qualifier with a warning
                    if ( num_seqs < sampling_size ) :
                        message = 'The number of sequencess is lower than the' \
                                  ' sampling size for the feature ' \
                                  '"{}"'.format(feat_key)
                        warnings.warn(message, RuntimeWarning, stacklevel=2)
                        continue
                    error = math.sqrt((num_seqs - sampling_size) / \
                                      sampling_size) * coef
                    threshold = math.ceil(sampling_size * error)
                    # For every existing pair of qualifiers, if the number of 
                    # records that hold both is below the calculated threshold,
                    # it might be the result of a typo in those records'
                    # information (further review of the log file is advisable)
                    for pair in itertools.combinations(qual_key.split(':'), 2) :
                        if ( (pair in term_dict[feat_key]) and
                             (len(term_dict[feat_key][pair]) <= threshold) ) :
                            log.write('\t{}\n\t\t{}\n'.format(' || '.join(pair),
                                          ', '.join(term_dict[feat_key][pair])))
                log.write('\n')
    # If no reference sequence has been introduced, include in the gene dict
    # those sequences that couldn't be processed due to lack of information
    if ( not ref_seq ) :
        set_dict['unprocessable'] = unprocessable
    return ( set_dict )


#-------------------------------------------------------------------------------