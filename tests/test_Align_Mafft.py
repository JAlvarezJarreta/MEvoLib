#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  test_Align_Mafft.py
# Last version :  v1.00 ( 16/Jul/2016 )
# Description :  Test suite for MEvoLib.Align module: Mafft software tool.
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  16/Jul/2016
#   VERSION :  v1.00
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

import os
import sys
import unittest

from Bio import SeqIO

from mevolib import MissingExtDependencyError
from mevolib import align
from mevolib._py3k import getoutput, viewkeys, viewitems


#-------------------------------------------------------------------------------

# Try to avoid problems when the OS is in another language
os.environ['LANG'] = 'C'

mafft_exe = None
if ( sys.platform == 'win32' ) :
    raise MissingExtDependencyError('Testing with MAFFT not implemented on ' \
                                    'Windows yet')
else :
    output = getoutput('mafft --help')
    if ( ('not found' not in output) and ('MAFFT' in output) ) :
        mafft_exe = 'mafft'
if ( not mafft_exe ) :
    raise MissingExtDependencyError('Install MAFFT if you want to use it from' \
                                    ' MEvoLib.Align.get_alignment().')


#-------------------------------------------------------------------------------

class AlignTestCase ( unittest.TestCase ) :

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
                Arguments passed to the alignment tool.
        """
        infile = '{}/f001.{}'.format(informat.capitalize(), informat)
        outfile = 'tmp_test.aln'
        self.add_file_to_clean(outfile)
        # Check the input
        self.assertTrue(os.path.isfile(infile))
        self.assertEqual(len(list(SeqIO.parse(infile, informat))), 50)
        # Generate the alignment
        align.get_alignment(mafft_exe, infile, informat, args=params,
                            outfile=outfile, outfile_format=outformat)
        # Check the output
        self.assertTrue(os.path.isfile(outfile))
        out_align = SeqIO.to_dict(SeqIO.parse(outfile, outformat))
        prevfile = '{}/f001.mafft_{}.aln'.format(outformat.capitalize(), params)
        self.assertTrue(os.path.isfile(prevfile))
        prev_align = SeqIO.to_dict(SeqIO.parse(prevfile, outformat))
        self.assertEqual(len(viewkeys(out_align)), len(viewkeys(prev_align)))
        for key, value in viewitems(out_align) :
            self.assertEqual(str(value.seq), str(prev_align[key].seq))


    def add_file_to_clean ( self, filename ) :
        """
        Adds a file for deferred removal by the tearDown() routine.

        Arguments :
            filename  ( string )
                File name to remove by the tearDown() routine.
        """
        self.files_to_clean.add(filename)



class MafftTestCase ( AlignTestCase ) :

    def test_simple_alignment ( self ) :
        """
        Test of the alignment method for all the available configurations with
        supported input and output formats.
        """
        for keyword in viewkeys(align.get_keywords(mafft_exe)) :
            self.standard_test('fasta', 'fasta', keyword)


    def test_conversion_alignment ( self ) :
        """
        Test of the alignment method with unsupported input and output formats
        (internal conversion required).
        """
        self.standard_test('genbank', 'phylip', 'default')


#-------------------------------------------------------------------------------

if ( __name__ == '__main__' ) :
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)


#-------------------------------------------------------------------------------
