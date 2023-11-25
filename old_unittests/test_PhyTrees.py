#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  test_PhyTrees.py
# Last version :  v1.00 ( 17/Jul/2016 )
# Description :  Test suite for MEvoLib.Fetch.PhyTrees module.
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

from Bio import Phylo

from mevolib.fetch.PhyTrees import PhyTrees
from mevolib._py3k import viewitems


#-------------------------------------------------------------------------------

class PhyTreesCreation ( unittest.TestCase ) :
    """
    Test PhyTrees reading and writing.
    """

    def setUp ( self ) :
        self.files_to_clean = set()


    def tearDown ( self ) :
        for filename in self.files_to_clean :
            if ( os.path.isfile(filename) ) :
                os.remove(filename)


    def test_treefile_source ( self ) :
        """
        Test PhyTrees.from_treefile() and PhyTrees.write() methods.
        """
        infile = 'Newick/f002.trees.newick'
        self.assertTrue(os.path.isfile(infile))
        tree_db = PhyTrees.from_treefile(infile, 'newick')
        outfile = 'tmp_test.newick'
        outrepfile = 'tmp_test.rep'
        self.files_to_clean.add(outfile)
        self.files_to_clean.add(outrepfile)
        tree_db.write(outfile)
        self.assertTrue(os.path.isfile(outfile))
        # Check the content of both sequence files
        self.assertEqual(len(list(Phylo.parse(infile, 'newick'))),
                         len(list(Phylo.parse(outfile, 'newick'))))
        # Check the content of the report file
        with open(outrepfile, 'r') as repfile :
            for line in repfile.readlines() :
                self.assertTrue(('Num. trees: 9' in line) or
                    ('History:' in line) or
                    (bool(re.match(r"""\d\d\d\d/\d\d/\d\d\ \d\d:\d\d:\d\d[ ]+
                                       [ ]+.*Tests/Newick/f002\.trees\.newick
                                       \ +newick""", line, re.VERBOSE))))


    def test_bioseqs_source ( self ) :
        """
        Test PhyTrees.from_bioseqs() method and len() property.
        """
        infile = 'PhyTrees/f002.trees.newick'
        inrepfile = 'PhyTrees/f002.trees.rep'
        self.assertTrue(os.path.isfile(infile))
        self.assertTrue(os.path.isfile(inrepfile))
        tree_db = PhyTrees.from_phytrees(infile)
        # Check the content of the PhyTrees' object
        self.assertEqual(len(tree_db), len(list(Phylo.parse(infile, 'newick'))))
        # Check the content of the PhyTrees' report
        with open(inrepfile, 'r') as repfile :
            line = repfile.readline().strip() # Num. trees: 9
            self.assertEqual(len(tree_db), int(line[-2:]))
            line = repfile.readline() # History:
            line = repfile.readline().strip() # [First source information]
            source_info = line.split('    ')
            self.assertEqual(tree_db._report[0], tuple(source_info))



class PhyTreesMethods ( unittest.TestCase ) :
    """
    Test PhyTrees methods.
    """

    def test_str ( self ) :
        """
        Test str() property.
        """
        infile = 'Newick/f002.trees.newick'
        self.assertTrue(os.path.isfile(infile))
        tree_db = PhyTrees.from_treefile(infile, 'newick')
        outstr = str(tree_db)
        # Check the output string
        lines = outstr.split('\n')
        self.assertEqual(lines[0], 'Num. trees: 9')
        self.assertEqual(lines[1], 'History:')
        self.assertTrue(bool(re.match(r"""\d\d\d\d/\d\d/\d\d\ \d\d:\d\d:\d\d[ ]+
            [ ]+.*Tests/Newick/f002\.trees\.newick[ ]+newick""", lines[2],
            re.VERBOSE)))


    def test_include ( self ) :
        """
        Test PhyTrees.include() method.
        """
        infile1 = 'Newick/f002.trees.newick'
        infile2 = 'Nexus/f005.trees.nexus'
        self.assertTrue(os.path.isfile(infile1))
        self.assertTrue(os.path.isfile(infile2))
        tree_db = PhyTrees.from_treefile(infile1, 'newick')
        tree_db.include(infile2, 'nexus')
        # Check the sequence data
        inlist1 = [tree  for tree in Phylo.parse(infile1, 'newick')]
        inlist2 = [tree  for tree in Phylo.parse(infile2, 'nexus')]
        self.assertEqual(len(inlist1) + len(inlist2), len(tree_db))
        # Check the report information
        self.assertIn('Tests/Newick/f002.trees.newick', tree_db._report[0][1])
        self.assertIn('newick', tree_db._report[0][2])
        self.assertIn('Tests/Nexus/f005.trees.nexus', tree_db._report[1][1])
        self.assertIn('nexus', tree_db._report[1][2])


    def test_join ( self ) :
        """
        Test PhyTrees.join() method.
        """
        infile1 = 'Newick/f002.trees.newick'
        infile2 = 'Nexus/f005.trees.nexus'
        self.assertTrue(os.path.isfile(infile1))
        self.assertTrue(os.path.isfile(infile2))
        tree_db = PhyTrees.from_treefile(infile1, 'newick')
        extra_db = PhyTrees.from_treefile(infile2, 'nexus')
        tree_db.join(extra_db)
        # Check the sequence data
        inlist1 = [tree  for tree in Phylo.parse(infile1, 'newick')]
        inlist2 = [tree  for tree in Phylo.parse(infile2, 'nexus')]
        self.assertEqual(len(inlist1) + len(inlist2), len(tree_db))
        # Check the report information
        self.assertIn('Tests/Newick/f002.trees.newick', tree_db._report[0][1])
        self.assertIn('newick', tree_db._report[0][2])
        self.assertIn('Tests/Nexus/f005.trees.nexus', tree_db._report[1][1])
        self.assertIn('nexus', tree_db._report[1][2])


    def test_statistics ( self ) :
        """
        Test PhyTrees.statistics() method.
        """
        infile = 'Newick/f002.trees.newick'
        self.assertTrue(os.path.isfile(infile))
        tree_db = PhyTrees.from_treefile(infile, 'newick')
        nseqs, mean, minimum, maximum = tree_db.statistics()
        # Check the resultant values
        self.assertEqual(nseqs, len(tree_db))
        self.assertEqual(mean, 10.0)
        self.assertEqual(minimum, 10)
        self.assertEqual(maximum, 10)


#-------------------------------------------------------------------------------

if ( __name__ == '__main__' ) :
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)


#-------------------------------------------------------------------------------
