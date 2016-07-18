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
# Description :  Python 3 compatibility tools. The inclusion of this module
#       provides full support for this tools in Python 2.7 and Python 3.3, 3.4,
#       3.5. (Base idea taken from Bio._py3k module of biopython)
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  16/Jul/2016
#   VERSION :  v1.10
#   AUTHOR(s) :  J. Alvarez-Jarreta
#   CHANGES :  * Added getouput() method for Windows and Unix systems.
#
#   DATE :  13/Jan/2016
#   VERSION :  v1.00
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

import sys
import os


#-------------------------------------------------------------------------------

if ( sys.version_info[0] >= 3 ) :
    # Code for Python 3
    from builtins import open, zip, map, filter, range, input
    import codecs, io

    def _is_int_or_long ( value ) :
        """
        Arguments :
            value  ( numeric type )

        Returns :
            bool
                True if value is instance of interger, False otherwise.

        * Note there are no longs on Python 3.
        """
        return ( isinstance(value, int) )


    def viewkeys ( dictionary ) :
        """
        Arguments :
            dictionary  ( dict )

        Returns :
            dict view
                Keys of 'dictionary'.
        """
        return ( dictionary.keys() )


    def viewvalues ( dictionary ) :
        """
        Arguments :
            dictionary  ( dict )

        Returns :
            dict view
                Values of 'dictionary'.
        """
        return ( dictionary.values() )


    def viewitems ( dictionary ) :
        """
        Arguments :
            dictionary  ( dict )

        Returns :
            dict view
                (key, value) pairs of 'dictionary'.
        """
        return ( dictionary.items() )


    # On Python 3 urllib, urllib2, and urlparse were merged
    from urllib.request import urlopen, Request, urlretrieve, urlparse
    from urllib.parse import urlencode, quote
    from urllib.error import HTTPError
    # On Python 3 subprocess.DEVNULL exists
    from subprocess import DEVNULL
    #On Python 3, this will be a unicode StringIO
    from io import StringIO
    from tempfile import TemporaryDirectory
else: # sys.version_info[0] < 3
    # Code for Python 2
    from __builtin__ import open, basestring, unicode
    # Import Python 3 like iterator functions:
    from future_builtins import zip, map, filter
    from __builtin__ import xrange as range
    from __builtin__ import raw_input as input

    def _is_int_or_long ( value ) :
        """
        Arguments :
            value  ( numeric type )

        Returns :
            bool
                True if value is instance of interger or long, False otherwise.
        """
        return ( isinstance(value, (int, long)) )


    def viewkeys ( dictionary ) :
        """
        Arguments :
            dictionary  ( dict )

        Returns :
            dict view
                Keys of 'dictionary'.
        """
        return ( dictionary.viewkeys() )


    def viewvalues ( dictionary ) :
        """
        Arguments :
            dictionary  ( dict )

        Returns :
            dict view
                Values of 'dictionary'.
        """
        return ( dictionary.viewvalues() )


    def viewitems ( dictionary ) :
        """
        Arguments :
            dictionary  ( dict )

        Returns :
            dict view
                (key, value) pairs of 'dictionary'.
        """
        return ( dictionary.viewitems() )


    # Under urllib.request on Python 3:
    from urllib2 import urlopen, Request
    from urllib import urlretrieve
    from urlparse import urlparse
    # Under urllib.parse on Python 3:
    from urllib import urlencode, quote
    # Under urllib.error on Python 3:
    from urllib2 import HTTPError
    # On Python 2 subprocess.DEVNULL doesn't exist
    DEVNULL = open(os.path.devnull, 'w')
    # On Python 2 this will be a (bytes) string based handle.
    # Note this doesn't work as it is unicode based:
    # from io import StringIO
    try :
        from cStringIO import StringIO
    except ImportError :
        from StringIO import StringIO
    try :
        input = raw_input
    except NameError :
        pass
    # There is no TemporaryDirectory class in the temp library
    from .TemporaryDirectory import TemporaryDirectory


if ( sys.platform == "win32" ) :
    # Can't use commands.getoutput on Python 2, Unix only/broken:
    # http://bugs.python.org/issue15073
    # Can't use subprocess.getoutput on Python 3, Unix only/broken:
    # http://bugs.python.org/issue10197
    def getoutput ( cmd ) :
        import subprocess
        child = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 universal_newlines=True, shell=False)
        stdout, stderr = child.communicate()
        # Remove trailing "\n" to match the Unix function
        return ( stdout.rstrip("\n") )
elif ( sys.version_info[0] >= 3 ):
    # Use subprocess.getoutput on Python 3
    from subprocess import getoutput
else : # sys.version_info[0] <= 2
    # Use commands.getoutput on Python 2
    from commands import getoutput

#-------------------------------------------------------------------------------
