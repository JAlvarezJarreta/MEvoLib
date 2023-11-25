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
# Last version :  v1.10 ( 16/Jul/2016 )
# Description :  MEvoLib's Data library.
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  16/Jul/2016
#   VERSION :  v1.10
#   AUTHOR(s) :  J. Alvarez-Jarreta
#   CHANGES :  * Added get_refseqs() method to have an easy access to all the
#                  available reference sequences.
#
#   DATE :  05/Feb/2016
#   VERSION :  v1.00
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

from __future__ import absolute_import

import os


#-------------------------------------------------------------------------------

def get_refseqs ( ) :
    """
    Returns :
        list
            List containing the name of all the available reference sequences.
    """
    return ( [fname[:-3]  for fname in os.listdir(os.path.dirname(__file__))
                              if ( fname.endswith('.py') and
                                   (fname != '__init__.py') )] )


#-------------------------------------------------------------------------------
