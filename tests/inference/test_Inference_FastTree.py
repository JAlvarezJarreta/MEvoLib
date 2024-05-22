from io import StringIO
import os
from pathlib import Path
import random

from Bio import Phylo
from Bio.Phylo.Consensus import _BitString
from Bio.Phylo import BaseTree
import pytest

from mevolib.inference import _FastTree as Fast


class MockStdOut:
    """
    Class defined to avoid making subprocess.run(command) every time we want to test the test_get_results function.
    As FastTree gives all the output in a whole file, we will store the results in a single file too, which will contain
    the resultant tree (after a command has been executed) and its associated score.
    Finally, we write this information in the matching output path, as the Inference get_results function expects a file
    to get some data from during its execution; and they have to be allocated in the temporary directory of the test
    class (and not in the data one) because every execution generates a couple of these two files with a random number
    in its path, so that storing and managing all of them could be quite expensive in memory terms.

    Arguments:
        tree_infile_path: Path where the resultant tree of the command execution is stored into.
        tree_outfile_path: Path where the resultant tree of the command execution and its score should be copied into; as the
        Inference module will need it to extract info from it.
        score: Associated phylogenetic score to the mentioned tree.
    """

    def __init__(self, tree_infile_path: Path, tree_outfile_path: Path, score: float):
        with open(tree_infile_path, "r") as fin:
            self.stdout = fin.read()
        with open(tree_outfile_path, "w") as fin:
            fin.write(self.stdout)
            fin.write(f"TreeLogLk	Placeholder	{score}")

    def get_mocked_output(self):
        return self.stdout


class TestInferenceFastTree:
    """
    Class made to ensure the correct operating of MEvoLib Inference FastTree module' functions.
    """

    tmp_dir: Path = Path("tests/fasttree_tmp_dir/").absolute()
    if not tmp_dir.exists():
        os.mkdir(tmp_dir)
    # Couple of functions used to compare Phylo Trees:

    """ 
    Divides a phylogenetic tree in function of the clades terminals, providing a way to compare trees.
    
    Arguments:
        tree: Phylo BaseTree object that wants to be transformed into bitstrings to get a comparison way.
    """

    def _bitstrs(self, tree: BaseTree):
        bitstrs = set()
        term_names = [term.name for term in tree.get_terminals()]
        term_names.sort()
        for clade in tree.get_nonterminals():
            clade_term_names = [term.name for term in clade.get_terminals()]
            boolvals = [name in clade_term_names for name in term_names]
            bitstr = _BitString("".join(map(str, map(int, boolvals))))
            bitstrs.add(bitstr)
        return bitstrs

    """ 
    Compares two phylogenetic trees and check they are "equal" (it is not 100% effective because of it 
    not being a char-by-char comparison; but still quite effective for the purposes of this library).
    
    Arguments:
        tree1: First Phylo BaseTree object that wants to be compared.
        tree2: Second Phylo BaseTree object that wants to be compared.
    """

    # Compare
    def compare(self, tree1: BaseTree, tree2: BaseTree):
        term_names1 = [term.name for term in tree1.get_terminals()]
        term_names2 = [term.name for term in tree2.get_terminals()]
        # false if terminals or BitStrings are not the same
        return set(term_names1) == set(term_names2) and self._bitstrs(tree1) == self._bitstrs(tree2)

    @pytest.mark.parametrize("format_list", [(["fasta", "phylip"])])
    def test_sprt_infile_formats(self, format_list: list):
        """
        Test function to check that FastTree module supports the required infile formats.

        Arguments:
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
        "args, infile_path, bootstraps, arg_list",
        [
            (
                "GTR+CAT",
                "tests/flatfiles/f001.mafft_default.aln",
                1,
                [
                    "-gtr",
                    "-nt",
                    "-nopr",
                    "-quiet",
                    "-log",
                    tmp_dir,
                    "-boot",
                    "1",
                    "tests/flatfiles/f001.mafft_default.aln",
                ],
            ),
            (
                "-nt -nopr",
                "tests/flatfiles/f001.mafft_default.aln",
                0,
                ["-nt", "-nopr", "-log", tmp_dir, "tests/flatfiles/f001.mafft_default.aln"],
            ),
        ],
    )
    def test_gen_args(self, args: str, infile_path: str, bootstraps: int, arg_list: list):
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
        "command, score, infile_path, expected_output",
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
                Path("f001.mafft_default.FastTree_output_1"),
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
                Path("f001.mafft_default.FastTree_output_2"),
            ),
        ],
    )
    def test_get_results(self, command: list, score: float, infile_path: Path, expected_output: Path):
        """
        Test function to ensure the correct extraction of the phylogeny and it's associated score from a command and the
        resultant output.  Ultimately, it ensures the correct execution of FastTree get_results's function.

        Arguments :
         command: FastTree's command line executed.
         score: The associated score an inferenced phylogeny tree has.
         infile_path: Input alignment file path.
         expected_output: Path of the output that the execution of the command would get.
        """
        # random temporary file path generation to save the results of the execution into
        r = random.randint(1, 999999)
        treefile_path = Path.joinpath(self.tmp_dir, "FastTree_bestTree." + str(r))

        command += ["-log", treefile_path, infile_path.absolute()]

        expected_output = pytest.data_dir / expected_output

        run_mocker = MockStdOut(expected_output, treefile_path, score)
        mocked_subprocess_output = run_mocker.get_mocked_output()

        phylogeny = Phylo.read(StringIO(mocked_subprocess_output), "newick")
        res_tree, res_score = Fast.get_results(command, mocked_subprocess_output)

        assert self.compare(phylogeny, res_tree)
        assert score == res_score

    @pytest.mark.parametrize(
        "command, tmp_file",
        [
            (
                [
                    "-log",
                    tmp_dir,
                ],
                None,
            ),
            (
                [
                    "-log",
                    tmp_dir,
                ],
                tmp_dir,
            ),
        ],
    )
    def test_cleanup(self, command: list, tmp_file: str):
        """
        Test function to ensure the correct removal of all the files of a given directory.

        Arguments :
            command: FastTree's command line executed.
            tmp_file: Path of the folder we want to delete (Just for testing purposes,
            because when called from __init__.py, the cleanup input is a bit different).
        """
        Fast.cleanup(command, tmp_file)
        assert not self.tmp_dir.exists()
