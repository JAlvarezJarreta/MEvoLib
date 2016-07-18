#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  BioSeqs.py
# Last version :  v1.11 ( 17/Jul/2016 )
# Description :  Definition and implementation of the class 'BioSeqs'.
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  17/Jul/2016
#   VERSION :  v1.11
#   AUTHOR(s) :  J. Alvarez-Jarreta
#   CHANGES :  * Minor bugs fixed.
#
#   DATE :  15/Jul/2016
#   VERSION :  v1.10
#   AUTHOR(s) :  J. Alvarez-Jarreta
#   CHANGES :  * Fixed an important bug in the BioSeqs.from_seqfile() method:
#                  the default alphabet assigned from Bio.SeqIO.read() and
#                  Bio.SeqIO.parse() methods when reading certain file formats
#                  (e.g. FASTA) raises an exception when writing those sequences
#                  with the Bio.SeqIO.write() method in a GENBANK file. Now a
#                  more specific alphabet is assigned to the Bio.Seq objects
#                  created.
#              * Fixed a bug in the BioSeqs.write() method: the unexpected
#                  exceptions were not raising after handling them.
#              * Fixed some documentation errors.
#
#   DATE :  25/Nov/2015
#   VERSION :  v1.00
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

from __future__ import absolute_import

import os
from datetime import datetime
import warnings
import sys
import itertools
from operator import itemgetter
import errno
import copy

import numpy
from Bio import SeqIO, Entrez, Alphabet

from MEvoLib._utils import get_abspath
from MEvoLib._py3k import StringIO, viewkeys, viewvalues, viewitems


#-------------------------------------------------------------------------------

# Dictionary with the corresponding retrieval type of each entrez database
# supported by BioSeqs class
entrez_db_dict = {'nuccore': 'gbwithparts', 'nucest': 'gb', 'nucgss': 'gb',
                  'popset': 'gb', 'protein': 'gp'}


#-------------------------------------------------------------------------------

def _get_entrez_db_rettype ( entrez_db ) :
    """
    Check if the given entrez database is supported by the BioSeqs class.
    Supported entrez databases are stored in 'entrez_db_dict'.

    Arguments :
        entrez_db  ( string )
            Entrez database name.

    Returns :
        string
            Corresponding retrieval type for the given entrez database.

    Raises :
        ValueError
            If the introduced database is not supported.
    """
    if ( entrez_db in entrez_db_dict ) :
        return ( entrez_db_dict[entrez_db] )
    else :
        message = '"{}" isn\'t a supported NCBI\'s Entrez DB'.format(entrez_db)
        raise ValueError(message)



def _estimate_batch_size ( record ) :
    """
    Arguments :
        record  ( string )
            Record of the sequence in text format.

    Returns :
        int
            Advisable batch size given the recommendations of Entrez and
            the record size.

    * It is calculated considering a 1Mbps download speed and a desirable
    maximum fetch time of 5 minutes per batch.
    """
    record_size = sys.getsizeof(record)
    # Entrez advises a ceiling of 200 identifiers per request, and return at
    # least 1 if the record is very big
    return ( int(max(min(200, numpy.floor(3.75e7/record_size)), 1)) )


#-------------------------------------------------------------------------------

class BioSeqs :
    """
    A BioSeqs object holds SeqRecord objects (Bio.SeqRecord from Biopython) and
    information about their source.
    
    Attributes :
        data  ( dict )
            Dictionary where the SeqRecord objects are stored.
        _report  ( Private[list] )
            List of tuples with the information of the source or sources of the
            SeqRecord objects stored at 'data'. This information is handled
            internally and shouldn't be modified externally in any way.

    Examples :

        Even though the user can create a BioSeqs object by hand, there are
        three main methods provided to make things easier:
        
        >>> from MEvoLib.BioSeqs import BioSeqs
        >>> seqsdb = BioSeqs.from_seqfile('seqs.fasta', 'fasta')
        >>> print(seqsdb)
        DB: {NC_012920, JT02934, HT23292, HGR8364, ASN29384,...}
        Num. sequences: 286
        History:
        2015/11/25 19:22:35    local    /home/usr1/seqs.fasta    fasta

        >>> bioseqsdb = BioSeqs.from_bioseqs('bioseqsfile.gb')
        >>> print(bioseqsdb)
        DB: {NC_012920, HT23292, HGR8364, ASN29384, MN73652,...}
        Num. sequences: 231
        History:
        2015/11/25 19:26:03    local    /home/usr1/old_seqs.fasta    fasta

        >>> hmtDNA = BioSeqs.from_entrez(entrez_db='nuccore',
        ... email='eg@test.com', max_fetch=10,
        ... query='"Homo sapiens"[Organism] AND mitochondrion[All Fields]')
        >>> print(hmtDNA)
        DB: {KR902533.1, KT002469.1, KR026958.1, KP300793.1, KR902534.1,...}
        Num. sequences: 10
        History:
        2015/11/25 19:34:47    entrez    nuccore    "Homo sapiens"[Organism] AND
        mitochondrion[All Fields]
        
        The user can also get several statistics (number of sequences, length
        mean, standard deviation of the length, minimum length, maximum length)
        about the stored data:
        
        >>> len(hmtDNA)
        10
        >>> hmtDNA.statistics()
        (10, 14926.0, 4924.3348992528927, 153, 16571)
        
        The database can be updated by joining two BioSeqs objects:
        
        >>> len(seqsdb)
        286
        >>> seqsdb.join(hmtDNA)
        >>> len(seqsdb)
        296
        
        Or including a new sequence file:

        >>> len(bioseqsdb)
        231
        >>> bioseqsdb.include('seqs.fasta', 'fasta')
        >>> len(bioseqsdb)
        286

        Or from the last "entrez" entry in the '_report' list:

        >>> hmtDNA.update()
        >>> len(hmtDNA)
        36633
        
        Finally, all the information can be saved in two files, one with the
        SeqRecord objects (in GENBANK format) and the other one with the
        information in the '_report' list:
        
        >>> seqsdb.write('new_seqs.gb')
        
    """

    def __init__ ( self, seq_dict, report ) :
        """
        Create a BioSeqs object with 'seq_dict' and 'report' information.

        Arguments :
            seq_dict  ( dict )
                Dictionary of sequences.
            report  ( list )
                Historical report of the sequence information.
        """
        self.data = seq_dict
        self._report = report



    @classmethod
    def from_bioseqs ( cls, bioseqs_file ) :
        """
        Create a BioSeqs object retrieving all the information from previously
        saved BioSeqs sequence and report files. If 'bioseqs_file' contains a
        relative path, the current working directory will be used to get the
        absolute path.
        
        Arguments :
            bioseqs_file  ( string )
                Sequence file generated by BioSeqs.write().

        Raises :
            ValueError
                If the number of sequences read doesn't match the number stored
                in the report document.
        """
        data_filepath = get_abspath(bioseqs_file)
        report_filepath = os.path.splitext(data_filepath)[0] + '.rep'
        # Load all the contents into a new BioSeqs object
        seq_dict = SeqIO.to_dict(SeqIO.parse(data_filepath, 'genbank'))
        report = []
        with open(report_filepath, 'r') as report_file :
            str_num_seqs = report_file.readline()
            num_seqs = int(str_num_seqs.split(':')[-1])
            if ( len(seq_dict) != num_seqs ) :
                message = 'The number of sequences at report file doesn\'t ' \
                          'match the number of sequences loaded'
                raise ValueError(message)
            # Ignore "History:" line
            report_file.readline()
            for line in report_file.readlines() :
                date_time, src_type, src, details = line.strip().split('    ')
                report.append((date_time, src_type, src, details))
        return ( cls(seq_dict, report) )



    @classmethod
    def from_seqfile ( cls, seqfile, fileformat ) :
        """
        Create a BioSeqs object retrieving all the information stored at the
        sequence file provided. If 'seqfile' contains a relative path, the
        current working directory will be used to get the absolute path.
        
        Arguments :
            seqfile  ( string )
                Input sequences file.
            fileformat  ( string ) 
                Input file format.

        Raises :
            IOError
                If the path or the file provided doesn't exist.

        * The file format must be supported by Bio.SeqIO.
        * If the file format provided doesn't correspond to the actual file
        format, an empty sequence dictionary will be created.
        """
        filepath = get_abspath(seqfile)
        # Read the sequence file and create a new BioSeqs object, generating a
        # new report list
        seq_dict = {}
        for record in SeqIO.parse(filepath, fileformat) :
            # When reading or parsing from certain sequence file format
            # (e.g. FASTA), Bio.SeqIO gives a default alphabet to the Seq object
            # created that will raise an error when writing it in a GENBANK
            # file. Thus, we change that alphabet to a more specific one,
            # checking if it is a DNA or a protein sequence
            if ( isinstance(record.seq.alphabet,
                            Alphabet.SingleLetterAlphabet) ) :
                record.seq.alphabet = Alphabet.IUPAC.ExtendedIUPACDNA()
                if ( not Alphabet._verify_alphabet(record.seq) ) :
                    record.seq.alphabet = Alphabet.IUPAC.ExtendedIUPACProtein()
            seq_dict[record.id] = record
        date_time = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
        report = [(date_time, 'local', filepath, fileformat)]
        return ( cls(seq_dict, report) )



    @classmethod
    def from_entrez ( cls, email, entrez_db, query, max_fetch = sys.maxsize ) :
        """
        Create a BioSeqs object fetching the sequences that matches the query
        at the provided database from NCBI's Entrez.
        
        Arguments :
            email  ( string )
                E-mail required by Bio.Entrez.
            entrez_db  ( string )
                NCBI's Entrez database.
            query  ( string )
                Query to fetch the desired information.
            max_fetch  ( Optional[int] )
                Maximum number of sequences to fetch.

        * The e-mail information is considered sensible information and it won't
        be saved in any public or private variable of the object.
        """
        Entrez.email = email
        db_rettype = _get_entrez_db_rettype(entrez_db)
        seq_dict = {}
        # Execute Entrez.esearch() to get the total number of sequences that
        # matches the query in the Entrez database
        handle = Entrez.esearch(db=entrez_db, term=query, rettype='count')
        available_seqs = int(Entrez.read(handle)['Count'])
        handle.close()
        num_seqs = min(max_fetch, available_seqs)
        if ( available_seqs < 1 ) :
            warnings.warn('The query provided didn\'t return any sequence')
        elif ( max_fetch < 1 ) :
            warnings.warn('"max_fetch" forces to create an empty dictionary')
        else : # num_seqs >= 1
            # Execute again Entrez.esearch() giving the total number of
            # sequences to get the complete list of Entrez database's sequence
            # identifiers
            sequence_ids = []
            for index in range(0, num_seqs, 100000) :
                handle = Entrez.esearch(db=entrez_db, term=query, restart=index,
                                        retmax=num_seqs)
                record = Entrez.read(handle)
                handle.close()
                sequence_ids += record['IdList']
            # Fetch the first sequence and estimate the batch size
            fetch_handle = Entrez.efetch(db=entrez_db, id=sequence_ids[0],
                                         retmode='text', rettype=db_rettype)
            record_str = fetch_handle.read()
            fetch_handle.close()
            record = SeqIO.read(StringIO(record_str), 'genbank')
            seq_dict[record.id] = record
            batch_size = _estimate_batch_size(record_str)
            # In batches of 'batch_size', fetch the Entrez database information
            # of each sequence in text format through Entrez.efetch()
            for start in range(1, num_seqs, batch_size) :
                end = start + batch_size
                fetch_handle = Entrez.efetch(db=entrez_db,
                                             id=sequence_ids[start:end],
                                             retmode='text', rettype=db_rettype)
                for record in SeqIO.parse(fetch_handle, 'genbank') :
                    seq_dict[record.id] = record
                fetch_handle.close()
        # Generate the corresponding report tuple
        date_time = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
        report = [(date_time, 'entrez', entrez_db, query)]
        return ( cls(seq_dict, report) )



    def __len__ ( self ) :
        """
        Returns :
            int
                Total number of sequences.
        """
        return ( len(self.data) )



    def __str__ ( self ) :
        """
        Returns :
            string
                Fancy print of the data stored in the BioSeqs object.
        """
        output = 'DB: {'
        end = min(5, len(self))
        for index, key in enumerate(itertools.islice(self.data.keys(), end)) :
            output += key
            if ( index < (end - 1) ) :
                output += ', '
        if ( len(self) != end ) :
            output += ',...'
        output += '}}\nNum. sequences: {:d}\nHistory:\n'.format(len(self))
        output += '\n'.join(['    '.join(x)  for x in self._report])
        return ( output )



    def include ( self, seqfile, fileformat ) :
        """
        Add the information of the sequence file to the BioSeqs object. For any
        matching sequence, the new information (from the external sequence file)
        will replace the existing one.
        
        Arguments :
            seqfile  ( string )
                Input sequences file.
            fileformat  ( string ) 
                Input file format.

        Raises :
            IOError
                If the path or the file provided doesn't exist.

        * The file format must be supported by Bio.SeqIO.
        * If the file format provided doesn't correspond to the actual file
        format, nothing will be added.
        """
        new_bioseqs = BioSeqs.from_seqfile(seqfile, fileformat)
        self.data.update(new_bioseqs.data)
        if ( new_bioseqs.data ) :
            self._report.extend(new_bioseqs._report)



    def join ( self, new_bioseqs ) :
        """
        Join the information of the BioSeqs object with another BioSeqs object.
        For any matching sequence, the new information (from the external
        BioSeqs object) will replace the existing one.
        
        Arguments :
            new_bioseqs  ( BioSeqs )
                The external BioSeqs object.
        """
        self.data.update(new_bioseqs.data)
        self._report.extend(new_bioseqs._report)
        self._report.sort(key=itemgetter(0))



    def update ( self, email ) :
        """
        Update the BioSeqs object from the last NCBI's Entrez database and query
        values stored in the report list. All the sequences stored must have
        their genbank identifier information in the annotations property. The
        deleted sequences from the database will be deleted in the object and
        the new sequences will be fetched and stored.
        
        Arguments :
            email  ( string )
                E-mail required by Bio.Entrez.

        Raises :
            ValueError
                If there is no entrez entry in the report list.
            ValueError
                If any sequence hasn't its GenBank identifier information in
                the annotations property.

        * The e-mail information is considered sensible information and it won't
        be saved in any public or private variable of the object.
        """
        # Get the last entrez entry of the report
        for record in reversed(self._report) :
            if ( record[1] == 'entrez' ) :
                date_time, src_type, entrez_db, query = record
                break
        else :
            message = 'No entrez entry found in object\'s report'
            raise ValueError(message)
        # Perform the update process in a copy of the dictionary to avoid
        # incomplete updates due to unexpected HTTP exceptions
        seq_dict = copy.copy(self.data)
        Entrez.email = email
        db_rettype = _get_entrez_db_rettype(entrez_db)
        # Execute Entrez.esearch() to get the total number of sequences that
        # matches the query in the Entrez database
        handle = Entrez.esearch(db=entrez_db, term=query, rettype='count')
        num_seqs = int(Entrez.read(handle)['Count'])
        handle.close()
        if ( num_seqs == 0 ) :
            warnings.warn('The query stored didn\'t return any sequence')
        else :
            # Execute again Entrez.esearch() giving the total number of
            # sequences to get the complete list of Entrez database's sequence
            # identifiers
            updated_seq_ids = set()
            for index in range(0, num_seqs, 100000) :
                handle = Entrez.esearch(db=entrez_db, term=query, restart=index,
                                        retmax=num_seqs)
                record = Entrez.read(handle)
                handle.close()
                updated_seq_ids.update(record['IdList'])
            # Get an "entrez identifier: accession" dictionary of the stored
            # sequences
            gi_acc_dict = {}
            try :
                for key, value in viewitems(seq_dict) :
                    gi_acc_dict[value.annotations['gi']] = key
            except KeyError :
                message = 'Missing genbank identifier'
                raise ValueError(message)
            else :
                # Use that dictionary to check which of the stored identifiers
                # have been removed from the Entrez database
                deprecated_seq_ids = viewkeys(gi_acc_dict)
                ids_to_remove = deprecated_seq_ids - updated_seq_ids
                # Remove all the deprecated sequences
                for gi_value in ids_to_remove :
                    accession = gi_acc_dict[gi_value]
                    del seq_dict[accession]
                    del gi_acc_dict[gi_value]
                # Finally, get the list of new identifiers to fetch
                ids_to_fetch = list(updated_seq_ids.difference(
                                                            deprecated_seq_ids))
                num_to_fetch = len(ids_to_fetch)
                if ( num_to_fetch > 0 ) :
                    # Fetch the first sequence and estimate the batch size
                    fetch_handle = Entrez.efetch(db=entrez_db,
                                                 id=ids_to_fetch[0],
                                                 retmode='text',
                                                 rettype=db_rettype)
                    record_str = fetch_handle.read()
                    fetch_handle.close()
                    record = SeqIO.read(StringIO(record_str), 'genbank')
                    seq_dict[record.id] = record
                    batch_size = _estimate_batch_size(record_str)
                    # In batches of 'batch_size', fetch the Entrez database
                    # information of each new sequence in text format through
                    # Entrez.efetch()
                    for start in range(1, num_to_fetch, batch_size) :
                        end = start + batch_size
                        fetch_handle = Entrez.efetch(db=entrez_db,
                                                     id=ids_to_fetch[start:end],
                                                     retmode='text',
                                                     rettype=db_rettype)
                        for record in SeqIO.parse(fetch_handle, 'genbank') :
                            seq_dict[record.id] = record
                        fetch_handle.close()
                    # The process has ended correctly so we can replace the old
                    # dictionary with the new one
                    self.data = seq_dict
                    # Generate the corresponding report tuple
                    date_time = datetime.now().strftime('%Y/%m/%d %H:%M:%S')
                    self._report.append((date_time, 'entrez', entrez_db, query))



    def write ( self, bioseqs_file ) :
        """
        Save all sequences stored at the BioSeqs object in the 'bioseqs_file'
        (in GENBANK format). A file with a detailed report of the sequences will
        be created replacing the extension of 'bioseqs_file' by ".rep". If
        'bioseqs_file' contains a relative path, the current working directory
        will be used to get the absolute path. If any file already exists, it
        will be overwritten without warning.
        
        Arguments :
            bioseqs_file  ( string )
                New BioSeqs sequence file.

        Raises :
            IOError
                If the path provided doesn't exist.
        """
        data_filepath = get_abspath(bioseqs_file)
        report_filepath = os.path.splitext(data_filepath)[0] + '.rep'
        # Generate a single string with all the report content
        str_report = '\n'.join(['    '.join(x)  for x in self._report])
        # Write all the information in the BioSeqs files
        try :
            SeqIO.write(viewvalues(self.data), data_filepath, 'genbank')
            with open(report_filepath, 'w') as report_file :
                report_file.write('Num. sequences: {:d}\nHistory:\n' \
                                  '{:s}'.format(len(self), str_report))
        except IOError :
            raise
        except :
            if ( os.path.lexists(data_filepath) ) :
                os.remove(data_filepath)
            if ( os.path.lexists(report_filepath) ) :
                os.remove(report_filepath)
            raise



    def statistics ( self ) :
        """
        Calculate the total number of sequences stored and the mean, standard
        deviation, minimum and maximum values of their length.
        
        Returns :
            int
                Total number of sequences.
            numpy.float64
                Mean length of the sequences.
            numpy.float64
                Standard deviation of the length of the sequences.
            numpy.int64
                Minimum length of the sequences.
            numpy.int64
                Maximum length of the sequences.
        """
        seqs_length = [len(record)  for record in self.data.values()]
        mean_value = numpy.mean(seqs_length)
        std_value = numpy.std(seqs_length)
        min_value = numpy.amin(seqs_length)
        max_value = numpy.amax(seqs_length)
        return ( len(self), mean_value, std_value, min_value, max_value )


#-------------------------------------------------------------------------------
