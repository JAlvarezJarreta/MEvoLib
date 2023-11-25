#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  _Muscle.py
# Last version :  v1.00 ( 13/Jan/2016 )
# Description :  MEvoLib's variables to ease the usage of Muscle.
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  13/Jan/2016
#   VERSION :  v1.00
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

SPRT_INFILE_FORMATS = ['fasta']

INFILE_CMD = '-in'

KEYWORDS = { 'default':  ['-quiet'],
             'fastdna':  ['-maxiters', '1', '-diags', '-quiet'],
             'fastprot': ['-maxiters', '1', '-diags', '-sv', '-distance1',
                          'kbit20_3', '-quiet'],
             'largein':  ['-maxiters', '2', '-quiet'] }


#-------------------------------------------------------------------------------
