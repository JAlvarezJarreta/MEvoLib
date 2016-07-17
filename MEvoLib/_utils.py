#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  _utils.py
# Last version :  v1.01 ( 20/Apr/2016 )
# Description :  Set of functions aimed to provide support for different needs
#       on the configuration and execution processes of MEvoLib.
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  20/Apr/2016
#   VERSION :  v1.01
#   AUTHOR(s) :  J. Alvarez-Jarreta
#   CHANGES :  * Added get_tempfile_path() method to generate a complete path
#                  for a temporary file without actually creating it.
#
#   DATE :  02/Feb/2016
#   VERSION :  v1.00
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

from __future__ import absolute_import

import os
import multiprocessing
import tempfile


#-------------------------------------------------------------------------------

NUMCORES = multiprocessing.cpu_count()


#-------------------------------------------------------------------------------

def get_abspath ( filename ) :
    """
    Return the absolute path of 'filename'. If 'filename' is already an absolute
    path, it is returned without changes. Otherwise, the current working
    directory is used as the file's absolute path.

    Arguments :
        filename  ( string )
            File name to get the absolute path from.

    Returns :
        string
            Absolute path of 'filename'.
    """
    if ( not os.path.isabs(filename) ) :
        return ( os.path.join(os.getcwd(), filename) )
    else :
        return ( filename )



def get_tempfile_path () :
    """
    Returns :
        string
            Path of a new temporary file name (without creating it).
    """
    return ( os.path.join(tempfile.gettempdir(), tempfile.gettempprefix() +\
                              next(tempfile._get_candidate_names())) )


#-------------------------------------------------------------------------------
