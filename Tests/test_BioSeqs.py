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
# Last version :  v1.0 ( 14/Jul/2016 )
# Description :  
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  14/Jul/2016
#   VERSION :  v1.0
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------
# Note :  The content of this file has been created using "run_tests.py" from
#   Biopython (http://www.biopython.org) as template.
#-------------------------------------------------------------------------------

import unittest

from MEvoLib.Fetch.BioSeqs import BioSeqs


#-------------------------------------------------------------------------------

class BioSeqsCreation ( unittest.TestCase ) :
    """

    """

    def test_valid_data ( self ) :
        pass


    def test_valid_report ( self ) :
        pass



class BioSeqsMethods ( unittest.TestCase ) :
    """
    Test BioSeqs methods.
    """

    def test_len ( self ) :
        pass


    def test_str ( self ) :
        pass


    def test_include ( self ) :
        pass


    def test_update ( self ) :
        pass


    def test_write ( self ) :
        pass


    def test_statistics ( self ) :
        pass


#-------------------------------------------------------------------------------

if ( __name__ == '__main__' ) :
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)


#-------------------------------------------------------------------------------
