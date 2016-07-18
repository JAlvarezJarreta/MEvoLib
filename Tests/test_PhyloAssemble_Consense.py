#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  test_PhyloAssemble_Consense.py
# Last version :  v1.00 ( 16/Jul/2016 )
# Description :  Test suite of MEvoLib.PhyloAssemble module: Consense software
#       tool.
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  16/Jul/2016
#   VERSION :  v1.00
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

import os
import unittest

from Bio import Phylo

from MEvoLib import MissingExtDependencyError
from MEvoLib import PhyloAssemble
from MEvoLib._py3k import getoutput, viewkeys


#-------------------------------------------------------------------------------

# Try to avoid problems when the OS is in another language
os.environ['LANG'] = 'C'

consense_exe = None
if ( sys.platform == 'win32' ) :
    raise MissingExtDependencyError('Testing with Consense not implemented on' \
                                    ' Windows yet')
else :
    output = getoutput('type consense')
    if ( 'not found' not in output ) :
        consense_exe = 'consense'
if ( not consense_exe ) :
    raise MissingExtDependencyError('Install Consense if you want to use it ' \
        'from MEvoLib.PhyloAssemble.get_consensus_tree().')


#-------------------------------------------------------------------------------

class PhyloAssembleTestCase ( unittest.TestCase ) :

    def setUp ( self ) :
        self.files_to_clean = set()


    def tearDown ( self ) :
        for filename in self.files_to_clean :
            if ( os.path.isfile(filename) ) :
                os.remove(filename)


    def standard_test ( self, informat, outformat, params ) :
        """
        Standard testing procedure used by all tests.

        Arguments :
            informat  ( string )
                Input file format.
            outformat  ( string )
                Output file format.
            params  ( string )
                Arguments passed to the consensus tree tool.
        """
        infile = '{}/f002.trees.{}'.format(informat.capitalize(), informat)
        outfile = 'tmp_test.tree'
        self.add_file_to_clean(outfile)
        # Check the input
        self.assertTrue(os.path.isfile(infile))
        self.assertEqual(len(list(Phylo.parse(infile, informat))), 9)
        # Generate the consensus tree
        PhyloAssemble.get_consensus_tree(consense_exe, infile, informat,
            args=params, outfile=outfile, outfile_format=outformat)
        # Check the output
        self.assertTrue(os.path.isfile(outfile))


    def add_file_to_clean ( self, filename ) :
        """
        Adds a file for deferred removal by the tearDown() routine.

        Arguments :
            filename  ( string )
                File name to remove by the tearDown() routine.
        """
        self.files_to_clean.add(filename)



class ConsenseTestCase ( PhyloAssembleTestCase ) :

    def test_simple_phylo_assembly ( self ) :
        """
        Test of the consensus tree method for all the available configurations
        with supported input and output formats.
        """
        for keyword in viewkeys(PhyloAssemble.get_keywords(consense_exe)) :
            self.standard_test('newick', 'newick', keyword)


    def test_conversion_phylo_assembly ( self ) :
        """
        Test of the consensus tree method with unsupported input and output
        formats (internal conversion required).
        """
        self.standard_test('nexus', 'nexus', 'default')


#-------------------------------------------------------------------------------

if ( __name__ == '__main__' ) :
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)


#-------------------------------------------------------------------------------
