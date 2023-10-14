from contextlib import nullcontext as does_not_raise
from pathlib import Path
from typing import ContextManager, Optional
import tempfile
import subprocess
import pytest

from mevolib.inference import _FastTree as Fast

class TestInferenceFastTree:
    
    tmp_dir: Path
    manifest_dir: Path

    @pytest.mark.parametrize(
            "format_list",
            [
                (["fasta","phylip"])
            ]
    )
    def test_sprt_infile_formats(self, format_list:list):
        assert len(Fast.SPRT_INFILE_FORMATS)==len(format_list)
        assert Fast.SPRT_INFILE_FORMATS==format_list

    @pytest.mark.parametrize(
            "keywords",
            [
                ({
                "default": ["-gtr", "-nt", "-nopr", "-quiet"],
                "GTR+CAT": ["-gtr", "-nt", "-nopr", "-quiet"],
                "GTR+G": ["-gtr", "-nt", "-nocat", "-gamma", "-nopr", "-quiet"],
                "JC+CAT": ["-nt", "-nopr", "-quiet"],
                "JC+G": ["-nt", "-nocat", "-gamma", "-nopr", "-quiet"],
                "JTT+CAT": ["-nopr", "-quiet"],
                "WAG+CAT": ["-wag", "-nopr", "-quiet"],
                }),
            ],
    )
    def test_keywords(self, keywords:dict):
        assert len(Fast.KEYWORDS)==len(keywords)
        for aKey, bKey in zip(Fast.KEYWORDS.keys(),keywords.keys()):
            assert aKey==bKey
        for key in Fast.KEYWORDS:
            assert Fast.KEYWORDS[key]==keywords[key]

    @pytest.mark.parametrize(
            "args, infile_path, bootstraps, log_tmpfile, arg_list",
            [
                (
                    "GTR+CAT",
                    "tests/Fasta/f001.mafft_default.aln",
                    1,
                    ["-gtr", "-nt", "-nopr", "-quiet", "-log", "tests/fasttree.log", "-boot", 1, "tests/Fasta/f001.mafft_default.aln"] ,
                    "tests/fasttree.log",
                    
                ),
                 
            ],

    )
    def test_gen_args(self, args: str, infile_path: str, bootstraps: int,
                        arg_list: list, log_tmpfile: str):
            
        assert arg_list==Fast.gen_args(args, infile_path, bootstraps,
                                        log_tmpfile)
    

    def test_myFunc(self): #command and output needed from __init__.py -> get_phylogeny(...) method to use
                           # in get_results(...) and cleanup(...) methods of _FastTree.py class
            
            #TODO: remove Fast.gen_args(...) and add fasttree -gtr -nt -nopr tests/Fasta/f001.mafft_default.aln
            #command = ["fasttree"]+ Fast.gen_args("-gtr, -nt, -nopr, -quiet -log", "my_test/test.py",0) 
            #subprocess.run(..., check=True, stdout=PIPE).stdout
            
            command = ["fasttree"]+ Fast.gen_args("GTR+CAT -log temp_file.txt", "output.newick",0) 
            # TODO: dont use genargs()
            print(command)
            output = subprocess.run(command, stderr=subprocess.DEVNULL,
                                            universal_newlines=True, check=True,
                                            stdout=subprocess.PIPE).stdout
            print(f"The output of the process is: {output}")     
            assert 0