#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  test_Cluster_NaiveRows.py
# Last version :  v1.00 ( 19/Jul/2016 )
# Description :  Test suite for MEvoLib.Cluster module: Genes method.
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  19/Jul/2016
#   VERSION :  v1.00
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

import os
import sys
import unittest

from Bio import SeqIO

from MEvoLib import MissingExtDependencyError
from MEvoLib import Cluster
from MEvoLib._py3k import getoutput, viewvalues, viewitems


#-------------------------------------------------------------------------------

# Try to avoid problems when the OS is in another language
os.environ['LANG'] = 'C'

mafft_exe = None
if ( sys.platform != 'win32' ) :
    output = getoutput('mafft --help')
    if ( ('not found' not in output) and ('MAFFT' in output) ) :
        mafft_exe = 'mafft'


#-------------------------------------------------------------------------------

class ClusterTestCase ( unittest.TestCase ) :

    def setUp ( self ) :
        self.files_to_clean = set()


    def tearDown ( self ) :
        for filename in self.files_to_clean :
            if ( os.path.isfile(filename) ) :
                os.remove(filename)


    def add_file_to_clean ( self, filename ) :
        """
        Adds a file for deferred removal by the tearDown() routine.

        Arguments :
            filename  ( string )
                File name to remove by the tearDown() routine.
        """
        self.files_to_clean.add(filename)



class GenesTestCase ( ClusterTestCase ) :

    def test_default ( self ) :
        """
        Testing default procedure for the genes method.
        """
        infile = 'Genbank/f006.genbank'
        informat = 'genbank'
        # Check the input
        self.assertTrue(os.path.isfile(infile))
        self.assertEqual(len(list(SeqIO.parse(infile, informat))), 5)
        # Generate the gene clustering
        subset_dict = Cluster.get_subsets('genes', infile, informat)
        # Check the output
        self.assertEqual(len(subset_dict), 63)
        self.assertIn('unprocessable', subset_dict)
        self.assertEqual(len(subset_dict['unprocessable']), 0)


    def test_feature_filter ( self ) :
        """
        Testing the genes method with a feature filter.
        """
        infile = 'Genbank/f006.genbank'
        informat = 'genbank'
        # Check the input
        self.assertTrue(os.path.isfile(infile))
        self.assertEqual(len(list(SeqIO.parse(infile, informat))), 5)
        # Generate the gene clustering
        subset_dict = Cluster.get_subsets('genes', infile, informat,
                                          feature_filter=['CDS'])
        # Check the output
        self.assertEqual(len(subset_dict), 14)
        for key, subset in viewitems(subset_dict) :
            if ( key == 'unprocessable' ) :
                self.assertEqual(len(subset), 0)
            else :
                self.assertEqual(len(subset), 5)


    def test_log_file ( self ) :
        """
        Testing the genes method with the generation of a log file in a given
        path.
        """
        infile = 'Genbank/f006.genbank'
        informat = 'genbank'
        logfile = 'tmp_test.log'
        # Check the input
        self.assertTrue(os.path.isfile(infile))
        self.assertEqual(len(list(SeqIO.parse(infile, informat))), 5)
        self.assertFalse(os.path.isfile(logfile))
        self.add_file_to_clean('tmp_test.log')
        # Generate the gene clustering
        subset_dict = Cluster.get_subsets('genes', infile, informat,
                                          log_file=logfile)
        # Check the clustering output
        self.assertEqual(len(subset_dict), 63)
        self.assertIn('unprocessable', subset_dict)
        self.assertEqual(len(subset_dict['unprocessable']), 0)
        # Check the content of the log file
        self.assertTrue(os.path.isfile(logfile))
        with open(logfile, 'r') as flog :
            content = flog.readlines()
            for feature in ['> misc_feature\n', '> D-loop\n', '> rRNA\n',
                            '> tRNA\n', '> CDS\n', '> gene\n'] :
                self.assertIn(feature, content)


    @unittest.skipIf(not mafft_exe,
                     'MAFFT software tool required for this test')
    def test_alignment ( self ) :
        """
        Testing default procedure for the genes method with alignment assitance
        for sequences without biological information (FASTA input instead of
        GENBANK).
        """
        infile = 'Fasta/f006.fasta'
        informat = 'fasta'
        # Check the input
        self.assertTrue(os.path.isfile(infile))
        self.assertEqual(len(list(SeqIO.parse(infile, informat))), 5)
        # Generate the gene clustering without metadata
        subset_dict = Cluster.get_subsets('genes', infile, informat)
        # Check the output
        self.assertEqual(len(subset_dict), 1)
        self.assertEqual(len(subset_dict['unprocessable']), 5)
        # Generate the gene clustering with external metadata (from a reference
        # sequence)
        subset_dict = Cluster.get_subsets('genes', infile, informat,
            ref_seq='rCRS', alignment_bin=mafft_exe)
        # Check the output
        self.assertEqual(len(subset_dict), 98)
        self.assertNotIn('unprocessable', subset_dict)
        for key, value in viewitems(subset_dict) :
            self.assertNotEqual(len(value), 0)
            self.assertTrue(len(value) % 5 == 0)
        


#-------------------------------------------------------------------------------

if ( __name__ == '__main__' ) :
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)


#-------------------------------------------------------------------------------
