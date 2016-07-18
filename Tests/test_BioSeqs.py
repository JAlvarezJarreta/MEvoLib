#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  test_BioSeqs.py
# Last version :  v1.00 ( 17/Jul/2016 )
# Description :  Test suite for MEvoLib.Fetch.BioSeqs module.
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  17/Jul/2016
#   VERSION :  v1.00
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

import os
import unittest
import re

from Bio import SeqIO, Entrez

from MEvoLib.Fetch.BioSeqs import BioSeqs
from MEvoLib._py3k import urlopen, URLError, viewitems


#-------------------------------------------------------------------------------

OFFLINE_MODE = False
# Check internet avaibility connecting to "google.com" through one of its
# IP addresses to avoid DNS lookup
try:
    response = urlopen('http://216.58.214.174', timeout=5)
except URLError :
    OFFLINE_MODE = True


#-------------------------------------------------------------------------------

class BioSeqsCreation ( unittest.TestCase ) :
    """
    Test BioSeqs reading and writing.
    """

    def setUp ( self ) :
        self.files_to_clean = set()


    def tearDown ( self ) :
        for filename in self.files_to_clean :
            if ( os.path.isfile(filename) ) :
                os.remove(filename)


    def test_seqfile_source ( self ) :
        """
        Test BioSeqs.from_seqfile() and BioSeqs.write() methods.
        """
        infile = 'Fasta/f001.fasta'
        self.assertTrue(os.path.isfile(infile))
        seq_db = BioSeqs.from_seqfile(infile, 'fasta')
        outfile = 'tmp_test.gb'
        outrepfile = 'tmp_test.rep'
        self.files_to_clean.add(outfile)
        self.files_to_clean.add(outrepfile)
        seq_db.write(outfile)
        self.assertTrue(os.path.isfile(outfile))
        # Check the content of both sequence files
        indict = SeqIO.to_dict(SeqIO.parse(infile, 'fasta'))
        outdict = SeqIO.to_dict(SeqIO.parse(outfile, 'gb'))
        self.assertEqual(len(indict), len(outdict))
        for key, value in viewitems(indict) :
            self.assertEqual(str(value.seq), str(outdict[key].seq))
        # Check the content of the report file
        with open(outrepfile, 'r') as repfile :
            for line in repfile.readlines() :
                self.assertTrue(('Num. sequences: 50' in line) or
                    ('History:' in line) or
                    (bool(re.match(r"""\d\d\d\d/\d\d/\d\d\ \d\d:\d\d:\d\d[ ]+
                                       local[ ]+.*Tests/Fasta/f001\.fasta
                                       [ ]+fasta""", line, re.VERBOSE))))


    def test_bioseqs_source ( self ) :
        """
        Test BioSeqs.from_bioseqs() method and len() property.
        """
        infile = 'BioSeqs/f001.gb'
        inrepfile = 'BioSeqs/f001.rep'
        self.assertTrue(os.path.isfile(infile))
        self.assertTrue(os.path.isfile(inrepfile))
        seq_db = BioSeqs.from_bioseqs(infile)
        # Check the content of the BioSeqs' object
        indict = SeqIO.to_dict(SeqIO.parse(infile, 'gb'))
        self.assertEqual(len(seq_db), len(indict))
        for key, value in viewitems(indict) :
            self.assertEqual(str(value.seq), str(seq_db.data[key].seq))
        # Check the content of the BioSeqs' report
        with open(inrepfile, 'r') as repfile :
            line = repfile.readline().strip() # Num. sequences: 50
            self.assertEqual(len(seq_db), int(line[-2:]))
            line = repfile.readline() # History:
            line = repfile.readline().strip() # [First source information]
            source_info = line.split('    ')
            self.assertEqual(seq_db._report[0], tuple(source_info))


    @unittest.skipIf(OFFLINE_MODE, 'Offline mode activated')
    def test_entrez_source ( self ) :
        """
        Test BioSeqs.from_entrez() method and len() property.
        """
        query = '"homo sapiens"[porgn] AND mitochondrion[Filter] AND ' \
                'mRNA[Filter]'
        seq_db = BioSeqs.from_entrez(entrez_db='nuccore', query=query,
            email='mevolib.test@zaramit.org', max_fetch=10)
        # Check the number of sequences fetched
        self.assertEqual(len(seq_db), 10)
        # Check report information
        self.assertIn('entrez', seq_db._report[0][1])
        self.assertIn('nuccore', seq_db._report[0][2])
        self.assertIn(query, seq_db._report[0][3])



class BioSeqsMethods ( unittest.TestCase ) :
    """
    Test BioSeqs methods.
    """

    def test_str ( self ) :
        """
        Test str() property.
        """
        infile = 'Fasta/f001.fasta'
        self.assertTrue(os.path.isfile(infile))
        seq_db = BioSeqs.from_seqfile(infile, 'fasta')
        outstr = str(seq_db)
        # Check the output string
        lines = outstr.split('\n')
        self.assertTrue(bool(re.match(r"""DB:\ \{clpA_\d+,\ clpA_\d+,\ clpA_\d+,
            \ clpA_\d+,\ clpA_\d+,\.\.\.\}""", lines[0], re.VERBOSE)))
        self.assertEqual(lines[1], 'Num. sequences: 50')
        self.assertEqual(lines[2], 'History:')
        self.assertTrue(bool(re.match(r"""\d\d\d\d/\d\d/\d\d\ \d\d:\d\d:\d\d[ ]+
            local[ ]+.*Tests/Fasta/f001\.fasta[ ]+fasta""", lines[3],
            re.VERBOSE)))


    def test_include ( self ) :
        """
        Test BioSeqs.include() method.
        """
        infile1 = 'Fasta/f001.fasta'
        infile2 = 'Phylip/f003.phylip'
        self.assertTrue(os.path.isfile(infile1))
        self.assertTrue(os.path.isfile(infile2))
        seq_db = BioSeqs.from_seqfile(infile1, 'fasta')
        seq_db.include(infile2, 'phylip')
        # Check the sequence data
        indict1 = SeqIO.to_dict(SeqIO.parse(infile1, 'fasta'))
        indict2 = SeqIO.to_dict(SeqIO.parse(infile2, 'phylip'))
        self.assertEqual(len(indict1) + len(indict2), len(seq_db))
        for key, value in viewitems(indict1) :
            self.assertEqual(str(value.seq), str(seq_db.data[key].seq))
        for key, value in viewitems(indict2) :
            self.assertEqual(str(value.seq), str(seq_db.data[key].seq))
        # Check the report information
        self.assertIn('local', seq_db._report[0][1])
        self.assertIn('Tests/Fasta/f001.fasta', seq_db._report[0][2])
        self.assertIn('fasta', seq_db._report[0][3])
        self.assertIn('local', seq_db._report[1][1])
        self.assertIn('Tests/Phylip/f003.phylip', seq_db._report[1][2])
        self.assertIn('phylip', seq_db._report[1][3])


    def test_join ( self ) :
        """
        Test BioSeqs.join() method.
        """
        infile1 = 'Fasta/f001.fasta'
        infile2 = 'Phylip/f003.phylip'
        self.assertTrue(os.path.isfile(infile1))
        self.assertTrue(os.path.isfile(infile2))
        seq_db = BioSeqs.from_seqfile(infile1, 'fasta')
        extra_db = BioSeqs.from_seqfile(infile2, 'phylip')
        seq_db.join(extra_db)
        # Check the sequence data
        indict1 = SeqIO.to_dict(SeqIO.parse(infile1, 'fasta'))
        indict2 = SeqIO.to_dict(SeqIO.parse(infile2, 'phylip'))
        self.assertEqual(len(indict1) + len(indict2), len(seq_db))
        for key, value in viewitems(indict1) :
            self.assertEqual(str(value.seq), str(seq_db.data[key].seq))
        for key, value in viewitems(indict2) :
            self.assertEqual(str(value.seq), str(seq_db.data[key].seq))
        # Check the report information
        self.assertIn('local', seq_db._report[0][1])
        self.assertIn('Tests/Fasta/f001.fasta', seq_db._report[0][2])
        self.assertIn('fasta', seq_db._report[0][3])
        self.assertIn('local', seq_db._report[1][1])
        self.assertIn('Tests/Phylip/f003.phylip', seq_db._report[1][2])
        self.assertIn('phylip', seq_db._report[1][3])


    @unittest.skipIf(OFFLINE_MODE, 'Offline mode activated')
    def test_update ( self ) :
        """
        Test BioSeqs.update() method.
        """
        query = '"homo sapiens"[porgn] AND mitochondrion[Filter] AND ' \
                'mRNA[Filter]'
        seq_db = BioSeqs.from_entrez(entrez_db='nuccore', email='eg@test.com',
                                     query=query, max_fetch=10)
        # Check the number of sequences fetched
        self.assertEqual(len(seq_db), 10)
        # Update the database fetching all the available sequences
        seq_db.update('mevolib.test@zaramit.org')
        # Check the number of sequences fetched
        handle = Entrez.esearch(db='nuccore', term=query, rettype='count')
        num_seqs = int(Entrez.read(handle)['Count'])
        self.assertEqual(len(seq_db), num_seqs)


    def test_statistics ( self ) :
        """
        Test BioSeqs.statistics() method.
        """
        infile = 'Fasta/f001.fasta'
        self.assertTrue(os.path.isfile(infile))
        seq_db = BioSeqs.from_seqfile(infile, 'fasta')
        nseqs, mean, stdev, minimum, maximum = seq_db.statistics()
        # Check the resultant values
        self.assertEqual(nseqs, len(seq_db))
        self.assertEqual(mean, 579.02)
        self.assertEqual(stdev, 0.14)
        self.assertEqual(minimum, 579)
        self.assertEqual(maximum, 580)


#-------------------------------------------------------------------------------

if ( __name__ == '__main__' ) :
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)


#-------------------------------------------------------------------------------
