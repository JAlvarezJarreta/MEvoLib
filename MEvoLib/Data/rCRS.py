#-------------------------------------------------------------------------------
#
#   MEvoLib  Copyright (C) 2016  J. Alvarez-Jarreta
#
#   This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#   This is free software, and you are welcome to redistribute it under certain
#   conditions; type `show c' for details.
#
#-------------------------------------------------------------------------------
# File :  rCRS.py
# Last version :  v1.0 ( 05/Feb/2016 )
# Description :  Human mitochondrial DNA revised Cambridge Reference Sequence
#       (rCRS)'s file path and Bio.SeqRecord object.
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  05/Feb/2016
#   VERSION :  v1.0
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

from __future__ import absolute_import

import os

from Bio import SeqIO


#-------------------------------------------------------------------------------

FILEPATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'rCRS.gb')

RECORD = SeqIO.read(FILEPATH, 'gb')


#-------------------------------------------------------------------------------
