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
# Description :  Test suite for MEvoLib.Cluster module: padded-Recursive-DCM3
#       (PRD) decomposition software tool.
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
from MEvoLib._py3k import getoutput, viewvalues


#-------------------------------------------------------------------------------

# Try to avoid problems when the OS is in another language
os.environ['LANG'] = 'C'

dcm3_exe = None
if ( sys.platform == 'win32' ) :
    raise MissingExtDependencyError('Testing with PRD not implemented on ' \
                                    'Windows yet')
else :
    output = getoutput('dcm')
    if ( ('not found' not in output) and ('dcm' in output) ) :
        dcm3_exe = 'dcm'
if ( not dcm3_exe ) :
    raise MissingExtDependencyError('Install DCM if you want to use it from' \
                                    ' MEvoLib.Cluster.get_subsets().')


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



class PRDTestCase ( ClusterTestCase ) :

    def test_clustering ( self ) :
        """
        Testing procedure for the PRD method.
        """
        infile = 'Fasta/f007.fasta'
        informat = 'fasta'
        treefile = 'Newick/f007.newick'
        treeformat = 'newick'
        # Check the input
        self.assertTrue(os.path.isfile(infile))
        self.assertEqual(len(list(SeqIO.parse(infile, informat))), 100)
        self.assertTrue(os.path.isfile(treefile))
        # Generate the subset division
        subset_dict = Cluster.get_subsets('prd', infile, informat,
            tree_file=treefile, file_format=treeformat, subset_size=25,
            overlapping=4, binary=dcm3_exe)
        # Check the output
        self.assertEqual(len(subset_dict), 17)
        result = [len(value)  for value in viewvalues(subset_dict)]
        result.sort()
        self.assertEqual(result, [16, 17, 17, 18, 18, 19, 19, 19, 20, 20, 20,
                                  20, 22, 22, 23, 24, 25])


#-------------------------------------------------------------------------------

if ( __name__ == '__main__' ) :
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)


#-------------------------------------------------------------------------------
