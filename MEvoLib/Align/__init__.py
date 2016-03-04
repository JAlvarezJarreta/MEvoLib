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
# Last version :  v1.0 ( 26/Jan/2016 )
# Description :  Functions aimed to provide an easy interface to handle
#       different alignment tools.
#-------------------------------------------------------------------------------
# Historical report :
#
#   DATE :  26/Jan/2016
#   VERSION :  v1.0
#   AUTHOR(s) :  J. Alvarez-Jarreta
#
#-------------------------------------------------------------------------------

from __future__ import absolute_import

import os
import tempfile
import subprocess

from Bio import SeqIO, AlignIO

from . import _ClustalOmega
from . import _Mafft
from . import _Muscle

from MEvoLib._utils import get_abspath
from MEvoLib._py3k import viewitems, DEVNULL, StringIO


#-------------------------------------------------------------------------------

_TOOL_TO_LIB = { 'clustalo': _ClustalOmega,
                 'mafft': _Mafft,
                 'muscle': _Muscle }


#-------------------------------------------------------------------------------

def get_keywords ( tool ) :
    """
    Arguments :
        tool  ( string )
            Name of the alignment tool.

    Returns :
        dict
            Dictionary containing the keywords and their corresponding
            arguments.

    Raises :
        ValueError
            If the tool introduced isn't included in MEvoLib.Align.
    """
    tool = tool.lower()
    if ( tool not in _TOOL_TO_LIB ) :
        message = 'The alignment tool "{}" isn\'t included in ' \
                  'MEvoLib.Align'.format(tool)
        raise ValueError(message)
    # else : # tool in _TOOL_TO_LIB
    keyword_dict = {}
    for key, value in iter(viewitems(_TOOL_TO_LIB[tool].KEYWORDS)) :
        keyword_dict[key] = ' '.join(value)
    return ( keyword_dict )
        



def get_alignment ( binary, infile, infile_format, args = 'default',
                    outfile = None, outfile_format = 'fasta', **kwargs ) :
    """
    Align the sequences of the input file using the alignment tool and arguments
    given. The resultant alignment is returned as a Bio.Align.MultipleSeqAlign
    object and saved in the output file (if provided). If 'infile' or 'outfile'
    contain a relative path, the current working directory will be used to get
    the absolute path. If the output file already exists, the old file will be
    overwritten without any warning.

    The alignment tool might not be included in MEvoLib, but it can still be
    used passing in '**kwargs' the keys "informats" and "incmd" with the list of
    of supported infile formats and the infile argument, respectively.

    Arguments :
        binary  ( string )
            Name or path of the alignment tool.
        infile  ( string )
            Unaligned input sequence file.
        infile_format  ( string )
            Input file format.
        args  ( Optional[string] )
            Keyword or arguments to use in the call of the alignment tool,
            excluding infile and outfile arguments.  By default, 'default'
            arguments are used.
        outfile  ( Optional[string] )
            Alignment output file.
        outfile_format  ( Optional[string] )
            Output file format. By default, FASTA format.
        **kwargs  ( Optional[dict] )
            Keyworded arguments required to execute alignment tools not included
            in the current version of MEvoLib. It is necessary to pass a list
            of supported infile formats under "informats" key, and the infile
            argument (e.g. "-in") with "incmd" key.

    Returns :
        Bio.Align.MultipleSeqAlignment
            Resultant alignment.

    Raises :
        IOError
            If the input path or the input file provided doesn't exist.
        RuntimeError
            If the call to the alignment tool command raises an exception.

    * The input file format must be supported by Bio.SeqIO.
    * The output file format must be supported by Bio.AlignIO.
    """
    # Get the variables associated with the given alignment tool, or get those
    # values from **kwargs
    bin_path, bin_name = os.path.split(binary)
    bin_name = bin_name.lower()
    if ( bin_name in _TOOL_TO_LIB ) :
        tool_lib = _TOOL_TO_LIB[bin_name]
        sprt_infile_formats = tool_lib.SPRT_INFILE_FORMATS
        infile_cmd = tool_lib.INFILE_CMD
        keywords = tool_lib.KEYWORDS
    else : # bin_name not in _TOOL_TO_LIB
        # Include the required variables through **kwargs dictionary
        sprt_infile_formats = kwargs['informats']
        infile_cmd = kwargs['incmd']
        keywords = dict()
    # Get the command line to run in order to get the resultant alignment
    infile_path = get_abspath(infile)
    # If the input file format is not supported by the alignment tool, convert
    # it to a temporary supported file
    if ( infile_format.lower() not in sprt_infile_formats ) :
        tmpfile = tempfile.NamedTemporaryFile()
        SeqIO.convert(infile_path, infile_format, tmpfile.name,
                      sprt_infile_formats[0])
        infile_path = tmpfile.name
    # Get argument list from keyword dictionary or 'args' string
    if ( args in keywords ) :
        arg_list = keywords[args]
    else : # args not in keywords
        # Remove possible empty strings in the given arguments
        arg_list = [arg  for arg in args.split(' ')]
    # Create full command line list (removing empty elements)
    command = [x  for x in [binary] + arg_list + [infile_cmd, infile_path] if x]
    # Run the alignment process handling any Runtime exception
    try :
        output = subprocess.check_output(command, stderr=DEVNULL,
                                         universal_newlines=True)
    except subprocess.CalledProcessError as e :
        message = 'Running "{}" raised an exception'.format(' '.join(e.cmd))
        raise RuntimeError(message)
    else :
        alignment = AlignIO.read(StringIO(output), 'fasta')
        if ( outfile ) :
            # Save the resultant alignment in the given outfile and format
            outfile_path = get_abspath(outfile)
            AlignIO.write(alignment, outfile_path, outfile_format)
        # Return the resultant alignment as a Bio.Align.MultipleSeqAligment
        # object
        return ( alignment )


#-------------------------------------------------------------------------------
