#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  __init__.py
# Last version :  v1.01 ( 16/Jul/2016 )
# Description :  Molecular Evolution Library for Python.
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  16/Jul/2016
#   VERSION :  v1.01
#   AUTHOR(s) :  J. Alvarez-Jarreta
#   CHANGES :  * Several bugs fixed.
#              * Minor functionalities added: Align.__init__.py,
#                  Cluster.__init__.py, Data.__init__.py, Inference.__init__.py,
#                  PhyloAssembly.__init__.py
#              * Documentation errors corrected.
#
#   DATE :  13/Jan/2016
#   VERSION :  v1.00
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

__version__ = '1.01'


#-------------------------------------------------------------------------------

class MissingExtDependencyError ( Exception ) :
    """
    Missing an external dependency. Used for our unit tests to allow skipping
    tests with missing external dependencies, e.g. missing command line tools.
    """
    pass


#-------------------------------------------------------------------------------
