from io import StringIO
from pathlib import Path
import random
import subprocess

from Bio import Phylo
from Bio.Phylo.Consensus import _BitString
import pytest

from mevolib.inference import _FastTree as Fast


class TestInferenceFastTree:
    tmp_dir: Path = Path("tests/fasttree_tmp_dir/").absolute()
    manifest_dir: Path
    # Couple of functions used to compare Phylo Trees:

    # store and return all _BitStrings
    def _bitstrs(self, tree):
        bitstrs = set()
        term_names = [term.name for term in tree.get_terminals()]
        term_names.sort()
        for clade in tree.get_nonterminals():
            clade_term_names = [term.name for term in clade.get_terminals()]
            boolvals = [name in clade_term_names for name in term_names]
            bitstr = _BitString("".join(map(str, map(int, boolvals))))
            bitstrs.add(bitstr)
        return bitstrs

    def compare(self, tree1, tree2):
        term_names1 = [term.name for term in tree1.get_terminals()]
        term_names2 = [term.name for term in tree2.get_terminals()]
        # false if terminals are not the same
        if set(term_names1) != set(term_names2):
            return False
        # true if _BitStrings are the same
        if self._bitstrs(tree1) == self._bitstrs(tree2):
            return True
        else:
            return False

    @pytest.mark.parametrize("format_list", [(["fasta", "phylip"])])
    def test_sprt_infile_formats(self, format_list: list):
        """
        Test function to check that FastTree module supports the required infile formats.
        Arguments :
            format_list: List of accepted infile formats for using FastTree tool.
        """
        assert len(Fast.SPRT_INFILE_FORMATS) == len(format_list)
        assert Fast.SPRT_INFILE_FORMATS == format_list

    @pytest.mark.parametrize(
        "keywords",
        [
            (
                {
                    "default": ["-gtr", "-nt", "-nopr", "-quiet"],
                    "GTR+CAT": ["-gtr", "-nt", "-nopr", "-quiet"],
                    "GTR+G": ["-gtr", "-nt", "-nocat", "-gamma", "-nopr", "-quiet"],
                    "JC+CAT": ["-nt", "-nopr", "-quiet"],
                    "JC+G": ["-nt", "-nocat", "-gamma", "-nopr", "-quiet"],
                    "JTT+CAT": ["-nopr", "-quiet"],
                    "WAG+CAT": ["-wag", "-nopr", "-quiet"],
                }
            ),
        ],
    )
    def test_keywords(self, keywords: dict):
        """
        Test function to check that FastTree module has the same keys and values for the keywords dictionary,
        that will be used later to pass arguments to a command.
        Arguments :
            keywords: Dictionary where keys are shortcuts for arguments combinations and values are lists of
            arguments a command may receive.
        """
        assert len(Fast.KEYWORDS) == len(keywords)
        for aKey, bKey in zip(Fast.KEYWORDS.keys(), keywords.keys()):
            assert aKey == bKey
        for key in Fast.KEYWORDS:
            assert Fast.KEYWORDS[key] == keywords[key]

    @pytest.mark.parametrize(
        "args, infile_path, bootstraps, log_tmpfile, arg_list",
        [
            (
                "GTR+CAT",
                "tests/Fasta/f001.mafft_default.aln",
                1,
                tmp_dir,
                [
                    "-gtr",
                    "-nt",
                    "-nopr",
                    "-quiet",
                    "-log",
                    tmp_dir,
                    "-boot",
                    "1",
                    "tests/Fasta/f001.mafft_default.aln",
                ],
            ),
            (
                "-nt -nopr",
                "tests/Fasta/f001.mafft_default.aln",
                0,
                tmp_dir,
                ["-nt", "-nopr", "-log", tmp_dir, "tests/Fasta/f001.mafft_default.aln"],
            ),
        ],
    )
    def test_gen_args(self, args: str, infile_path: str, bootstraps: int, log_tmpfile: str, arg_list: list):
        """
        Test function to ensure the construction of an argument list taking into account the parameters a
        command might need, and how it needs them. Ultimately, it ensures the correct execution of FastTree
        gen_args's function.
        Arguments :
            args: Keyword or arguments to use in the call of FastTree, excluding
                infile and outfile arguments.
            infile_path: Input alignment file path.
            bootstraps: Number of bootstraps to generate.
            log_tmpfile: Path of the directory we want to save the FastTree output data into (Only required while testing for
                reproducibility purposes or just if the user wants a specific folder to allocate the result of the execution
                into. Otherwise, a random one will be provided).
            arg_list: Handmade constructed argument list as the expected result the real function should return.
        """
        assert arg_list == Fast.gen_args(args, infile_path, bootstraps, self.tmp_dir)

    @pytest.mark.parametrize(
        "command, score, infile_path",
        [
            (
                [
                    "fasttree",
                    "-gtr",
                    "-nt",
                    "-nopr",
                    "-quiet",
                ],
                -1911.868,
                Path("tests/Fasta/f001.mafft_default.aln"),
            ),
            (
                [
                    "fasttree",
                    "-gtr",
                    "-nt",
                    "-nocat",
                    "-gamma",
                    "-nopr",
                    "-quiet",
                    "-boot",
                    "1",
                ],
                -2008.6733,
                Path("tests/Fasta/f001.mafft_default.aln"),
            ),
        ],
    )
    def test_get_results(self, command: list, score: float, infile_path: Path):
        """
        Test function to ensure the correct extraction of the phylogeny and it's associated score from a command and the
        resultant output.  Ultimately, it ensures the correct execution of FastTree get_results's function.
        Arguments :
         command: FastTree's command line executed.
         score: The associated score an inferenced phylogeny tree has.
         infile_path:  Input alignment file path.
        """
        # random temporary file path generation to save the results of the execution into
        r = random.randint(1, 999999)
        treefile_path = Path.joinpath(self.tmp_dir, "FastTree_bestTree." + str(r))

        command += ["-log", treefile_path, infile_path.absolute()]
        output = subprocess.run(
            command, stderr=subprocess.DEVNULL, universal_newlines=True, check=True, stdout=subprocess.PIPE
        ).stdout
        phylogeny = Phylo.read(StringIO(output), "newick")
        res_tree, res_score = Fast.get_results(command, output)

        assert self.compare(phylogeny, res_tree)
        assert score == res_score
