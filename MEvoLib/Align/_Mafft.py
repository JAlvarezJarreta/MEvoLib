#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  _Mafft.py
# Last version :  v1.01 ( 01/Feb/2016 )
# Description :  MEvoLib's variables to ease the usage of Mafft.
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  01/Feb/2016
#   VERSION :  v1.01
#   AUTHOR(s) :  J. Alvarez-Jarreta
#   CHANGES :  * The "-1" option for "--thread" argument doesn't work as
#                  expected, so we now include explicitly the number of cores.
#
#   DATE :  13/Jan/2016
#   VERSION :  v1.00
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

from __future__ import absolute_import

from MEvoLib._utils import NUMCORES


#-------------------------------------------------------------------------------

SPRT_INFILE_FORMATS = ['fasta']

INFILE_CMD = ''

KEYWORDS = { 'default':  ['--auto', '--quiet', '--thread', str(NUMCORES)],
             'linsi':    ['--localpair', '--maxiterate', '1000', '--quiet',
                          '--thread', str(NUMCORES)],
             'parttree': ['--parttree', '--retree', '2', '--partsize', '1000',
                          '--quiet', '--thread', str(NUMCORES)] }


#-------------------------------------------------------------------------------
