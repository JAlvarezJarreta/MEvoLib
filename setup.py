#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  setup.py
# Last version :  v1.10 ( 16/Jul/2016 )
# Description :  Distutils based setup script for MEvoLib.
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  16/Jul/2016
#   VERSION :  v1.10
#   AUTHOR(s) :  J. Alvarez-Jarreta
#   CHANGES :  * Test suite included.
#
#   DATE :  12/Feb/2016
#   VERSION :  v1.00
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------
# Note :  The content of this file has been created using "setup.py" from
#   Biopython (http://www.biopython.org) as template.
#-------------------------------------------------------------------------------

from __future__ import print_function

import sys
import os
from distutils.core import setup
from distutils.core import Command
from distutils.command.install import install
from distutils.command.build_py import build_py


#-------------------------------------------------------------------------------

_CHECKED = None


#-------------------------------------------------------------------------------

def can_import ( module_name ) :
    """
    Check whether the 'module_name' can be imported or not.
    
    Arguments :
        module_name  ( string )
            Name of the module to be imported.

    Returns :
        bool
            True if 'module_name' can be imported, False otherwise.
    """
    try :
        __import__(module_name)
    except ImportError:
        return ( False )
    else :
        return ( True )



def check_dependencies ( ) : 
    """
    Return whether the installation should continue.

    Returns :
        bool
            True if it can continue, False otherwise.
    """
    # Check if NumPy or Biopython are missing, as they are required for MEvoLib
    # to work properly
    if ( not can_import('numpy') ) :
        print('Numerical Python (NumPy) is not installed.\nThis package is ' \
              'required for Biopython and MEvoLib.\n\nYou can find NumPy at' \
              ' http://www.numpy.org')
        return ( False )
    if ( not can_import('Bio') ) :
        print('Biopython is not installed.\nThis package is required for' \
              'MEvoLib.\n\nYou can find Biopython at http://www.biopython.org')
        return ( False )
    if ( not can_import('dendropy') ) :
        print('Dendropy is not installed.\nThis package is required for' \
              'MEvoLib.\n\nYou can find Dendropy at http://www.dendropy.org/')
        return ( False )
    # Exit automatically if running as part of some script
    if ( not sys.stdout.isatty() ) :
        sys.exit(-1)
    return ( True )



def check_dependencies_once ( ) :
    """
    Call 'check_dependencies', caching the result to avoid subsequent calls.

    Returns :
        bool
            True if all the dependencies have been satisfied, False otherwise.
    """
    global _CHECKED
    if ( _CHECKED is None ) :
        _CHECKED = check_dependencies()
    return ( _CHECKED )



class install_mevolib ( install ) :
    """
    Override the standard install to check for dependencies.
    """

    def run(self):
        if ( check_dependencies_once() ) :
            # Run the normal install
            install.run(self)



class build_py_mevolib ( build_py ) :
    """
    Override the standard build to check for dependencies.
    """

    def run(self):
        if ( not check_dependencies_once() ) :
            return
        build_py.run(self)



class test_mevolib ( Command ):
    """
    Run all of the tests for MEvoLib. This is a automatic test run class to make
    distutils kind of act like perl. With this you can do:

        python setup.py build
        python setup.py install
        python setup.py test

    """
    description = "Automatically run the test suite for MEvoLib."
    user_options = []

    def initialize_options(self):
        pass


    def finalize_options(self):
        pass


    def run(self):
        this_dir = os.getcwd()
        # Change to the test dir and run the tests
        os.chdir('Tests')
        sys.path.insert(0, '')
        import run_tests
        run_tests.main([])
        # Change back to the current directory
        os.chdir(this_dir)


#-------------------------------------------------------------------------------

# Check that we have the right Python version
if ( sys.version_info[:2] < (2, 7) ) :
    print('MEvoLib requires Python 2.7 (or Python 3.3 or later). ' \
          'Python %d.%d detected'.format(sys.version_info[:2]))
    sys.exit(1)
elif ( (sys.version_info[0] == 3) and (sys.version_info[:2] < (3, 3)) ) :
    print('MEvoLib requires Python 3.3 or later (or Python 2.7). ' \
          'Python %d.%d detected' % sys.version_info[:2])
    sys.exit(1)

# We now define the MEvoLib version number in MEvoLib/__init__.py
__version__ = 'unknown'
for line in open('MEvoLib/__init__.py') :
    if ( line.startswith('__version__') ) :
        exec(line.strip())

# Simple trick to use the 2to3 converted source under Python 3, change the
# current directory before/after running setup. Note as a side effect there will
# be a build folder underneath the python3_source folder.
old_path = os.getcwd()
try:
    src_path = python3_source
except NameError :
    src_path = os.path.dirname(os.path.abspath(sys.argv[0]))
os.chdir(src_path)
sys.path.insert(0, src_path)

setup_args = { 'name': 'MEvoLib',
               'version': '1.0',
               'author': 'J. Alvarez-Jarreta',
               'author_email': 'jorgeal@unizar.es',
               'url': 'http://zaramit.org/mevolib/',
               'description': 'MEvoLib is molecuar evolution library of free ' \
                              'tools for Python 2.7 and Python 3.3 or newer',
               'download_url': 'http://zaramit.org/downloads/' \
                               'MEvoLib_v1.0.tar.gz',
               'cmdclass': { 'install': install_mevolib,
                             'build_py': build_py_mevolib,
                             'test': test_mevolib, },
               'packages': [ 'MEvoLib',
                             'MEvoLib.Align',
                             'MEvoLib.Cluster',
                             'MEvoLib.Data',
                             'MEvoLib.Fetch',
                             'MEvoLib.Inference',
                             'MEvoLib.PhyloAssemble',
                             'MEvoLib._py3k' ],
               'package_data': {'MEvoLib.Data': ['*.gb']}, }

try:
    setup(**setup_args)
finally:
    del sys.path[0]
    os.chdir(old_path)

