#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  run_tests.py
# Last version :  v1.00 ( 17/Jul/2016 )
# Description :  Run a set of PyUnit-based regression tests.
#       This will find all modules whose name is "test_*.py" in the test
#       directory, and run them. Various command line options provide additional
#       facilities. By default, all tests are run.
#
#       Command line options:
#       --help        -- show usage info
#       -v;--verbose  -- run tests with higher verbosity (does not affect our
#                        print-and-compare style unit tests).
#       <test_name>   -- supply the name of one (or more) tests to be run.
#                        The .py file extension is optional.
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  17/Jul/2016
#   VERSION :  v1.00
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------
# Note :  The content of this file has been created using "run_tests.py" from
#   Biopython (http://www.biopython.org) as template.
#-------------------------------------------------------------------------------
"""
Run a set of PyUnit-based regression tests.

Command line options:
--help        -- show usage information
-v;--verbose  -- run tests with higher verbosity (does not affect our
                 print-and-compare style unit tests).
<test_name>   -- supply the name of one (or more) tests to be run. The .py file
                 extension is optional.
"""
#-------------------------------------------------------------------------------

from __future__ import print_function

import os
import sys
import time
import getopt
import unittest
import distutils.util
import traceback

# If we want to be able to call "run_tests.py" BEFORE MEvoLib is installed, we
# can't use this:
#     from MEvoLib._py3k import StringIO
try:
    from StringIO import StringIO  # Python 2 (byte strings)
except ImportError:
    from io import StringIO  # Python 3 (unicode strings)


#-------------------------------------------------------------------------------

# The default verbosity (not verbose)
VERBOSITY = 0

# Cache the system language
SYSTEM_LANG = os.environ.get('LANG', 'C')


#-------------------------------------------------------------------------------

def main ( argv ) :
    """
    Run tests.

    Arguments :
        argv  ( list )
            Arguments used when calling to "run_test.py".

    Returns :
        int
            Number of failures.
    """
    # Insert MEvoLib's paths in sys.path: '../build/lib.*' and '..'
    test_path = sys.path[0] or '.'
    source_path = os.path.abspath('{}/..'.format(test_path))
    sys.path.insert(1, source_path)
    build_path = os.path.abspath('{}/../build/lib.{}-{}'.format(
        test_path, distutils.util.get_platform(), sys.version[:3]))
    if ( os.access(build_path, os.F_OK) ) :
        sys.path.insert(1, build_path)
    # Get the command line options
    try:
        opts, args = getopt.getopt(argv, 'v', ['verbose', 'help', 'offline'])
    except getopt.error as msg :
        print(msg)
        print(__doc__)
        return ( 2 )
    verbosity = VERBOSITY
    # Process the options
    for option, values in opts :
        if ( option == '--help' ) :
            print(__doc__)
            return ( 0 )
#        if ( option == '--offline' ) :
#            print("Skipping any test requiring internet access")
#            OFFLINE_MODE = True
        if ( option in ['-v', '--verbose'] ) :
            verbosity = 2
    # The 'args' variable should contain the names of the tests to run
    for index in range(len(args)):
        # Remove the ".py" extension if it was included
        if ( args[index].endswith('.py') ) :
            args[index] = args[index][:-3]
    print('Python version: {}'.format(sys.version))
    print('Operating system: {} {}'.format(os.name, sys.platform))
    # Run the tests
    runner = TestRunner(args, verbosity)
    return ( runner.run() )



class TestRunner ( unittest.TextTestRunner ) :

    if ( __name__ == '__main__' ) :
        file = sys.argv[0]
    else :
        file = __file__
    testdir = os.path.abspath(os.path.dirname(file) or os.curdir)


    def __init__ ( self, tests = (), verbosity = 0 ) :
        self.tests = tests
        # If no tests were specified, run them all
        if ( not self.tests ) :
            # Make a list of all applicable test modules
            names = os.listdir(TestRunner.testdir)
            for name in names :
                if ( name.startswith('test_') and name.endswith('.py') ) :
                    self.tests.append(name[:-3])
            self.tests.sort()
        stream = StringIO()
        unittest.TextTestRunner.__init__(self, stream, verbosity = verbosity)


    def runTest ( self, name ):
        from mevolib import MissingExtDependencyError
        result = self._makeResult()
        output = StringIO()
        # Restore the language and thus default encoding (in case a prior
        # test changed this, e.g. to help with detecting command line tools)
        global SYSTEM_LANG
        os.environ['LANG'] = SYSTEM_LANG
        # Always run tests from the "Tests" folder where "run_tests.py"
        # should be located (as we assume this with relative paths, etc.)
        os.chdir(self.testdir)
        try:
            stdout = sys.stdout
            sys.stdout = output
            if ( name.startswith('test_') ) :
                sys.stderr.write('{} ... '.format(name))
                # It's either a unittest or a print-and-compare test
                loader = unittest.TestLoader()
                suite = loader.loadTestsFromName(name)
                if ( hasattr(loader, 'errors') and loader.errors ) :
                    # New in Python 3.5: we don't always get an exception.
                    # Instead this is a list of error messages as strings
                    for msg in loader.errors :
                        if ( 'MEvoLib.MissingExtDependencyError: ' in msg ) :
                            # Remove the traceback, etc.
                            msg = msg[msg.find('MEvoLib.Missing'):]
                            msg = msg[msg.find('Error: '):]
                            sys.stderr.write('skipping. {}\n'.format(msg))
                            return ( True )
                    # Looks like a real failure
                    sys.stderr.write('loading tests failed:\n')
                    for msg in loader.errors :
                        sys.stderr.write('{}\n'.format(msg))
                    return ( False )
            else : # wrong test name (they should all start with "test_")
                sys.stderr.write('ERROR: Wrong test name "{}"\n'.format(name))
            # Run all tests cases
            suite.run(result)
            if ( self.testdir != os.path.abspath('.') ) :
                sys.stderr.write('FAIL\n')
                result.stream.write(result.separator1 + '\n')
                result.stream.write('ERROR: {}\n'.format(name))
                result.stream.write(result.separator2 + '\n')
                result.stream.write('Current directory changed\n')
                result.stream.write('Was: {}\n'.format(self.testdir))
                result.stream.write('Now: {}\n'.format(os.path.abspath('.')))
                os.chdir(self.testdir)
                if ( not result.wasSuccessful() ) :
                    result.printErrors()
                return ( False )
            elif ( result.wasSuccessful() ) :
                sys.stderr.write('ok\n')
                return ( True )
            else : # something went wrong
                sys.stderr.write('FAIL\n')
                result.printErrors()
            return ( False )
        except MissingExtDependencyError as msg :
            # Seems this isn't always triggered on Python 3.5, exception
            # messages can be in 'loader.errors' instead
            sys.stderr.write('skipping. {}\n'.format(msg))
            return ( True )
        except Exception as msg :
            # This happened during the import
            sys.stderr.write('ERROR\n')
            result.stream.write(result.separator1 + '\n')
            result.stream.write('ERROR: {}\n'.format(name))
            result.stream.write(result.separator2 + '\n')
            result.stream.write(traceback.format_exc())
            return ( False )
        except KeyboardInterrupt as err :
            # We want to allow this and abort the test
            raise err
        except :
            return ( False )
        finally :
            sys.stdout = stdout


    def run ( self ) :
        """
        Run tests.

        Returns :
            int
                Number of failures.
        """
        failures = 0
        startTime = time.time()
        for test in self.tests :
            ok = self.runTest(test)
            if ( not ok ) :
                failures += 1
        total = len(self.tests)
        stopTime = time.time()
        timeTaken = stopTime - startTime
        sys.stderr.write(self.stream.getvalue())
        sys.stderr.write('-' * 70 + '\n')
        sys.stderr.write('Ran {} test{} in {:.3f} seconds\n'.format(total,
            (total != 1 and 's' or ''), timeTaken))
        sys.stderr.write('\n')
        if ( failures ) :
            sys.stderr.write('FAILED ({} failure{})\n'.format(failures,
                (failures != 1 and 's' or '')))
        return ( failures )


#-------------------------------------------------------------------------------

if ( __name__ == '__main__' ) :
    errors = main(sys.argv[1:])
    if ( errors ) :
        # Doing a sys.exit(...) isn't nice if run from IDLE...
        sys.exit(1)


#-------------------------------------------------------------------------------
