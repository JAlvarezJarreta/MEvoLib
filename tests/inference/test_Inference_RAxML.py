import os
from pathlib import Path
import random

from Bio import Phylo
from Bio.Phylo.Consensus import _BitString
from Bio.Phylo import BaseTree
import pytest

from mevolib.inference import _RAxML as Rax
from mevolib._utils import NUMCORES


class MockStdOut:
    """
    Class defined to avoid making subprocess.run(command) every time we want to test the test_get_results function.
    As RAxML(HPC) gives the output splitted into some files, we will store the results in two main files: the resultant
    tree itself (whose value is in stdout_tree) and the info associated to it (whose value is in stdout_info).
    Finally, we write this information in the matching output path, as the Inference get_results function expects a file
    to get some data from during its execution; and they have to be allocated in the temporary directory of the test
    class (and not in the data one) because every execution generates a couple of these two files with a random number
    in its path, so that storing and managing all of them could be quite expensive in memory terms.

    Arguments:
        tree_infile_path: Path where the resultant tree of the command execution is stored into.
        tree_outfile_path: Path where the resultant tree of the command execution should be copied into; as the
        Inference module will need it to extract info from it.
        info_infile_path: Path where the information of the resultant tree of the command execution is stored into.
        info_outfile_path: Path where the information of the resultant tree of the command execution should be copied
        into; as the Inference module will need it to extract info from it.
    """

    def __init__(
        self, tree_infile_path: Path, tree_outfile_path: Path, info_infile_path: Path, info_outfile_path: Path
    ):
        with open(tree_infile_path, "r") as fin:
            self.stdout_tree = fin.read()
        with open(tree_outfile_path, "w") as fin:
            fin.write(self.stdout_tree)
        with open(info_infile_path, "r") as fin:
            self.stdout_info = fin.read()
        with open(info_outfile_path, "w") as fin:
            fin.write(self.stdout_info)

    def get_mocked_tree_output(self):
        return self.stdout_tree

    def get_mocked_info_output(self):
        return self.stdout_info


class TestInferenceRAxML:
    """
    Class made to ensure the correct operating of MEvoLib Inference RAxML module' functions.
    """

    tmp_dir: Path = Path("tests/raxml_tmp_dir/").absolute()
    if not tmp_dir.exists():
        os.mkdir(tmp_dir)
    # Couple of functions used to compare Phylo Trees:

    # store and return all _BitStrings
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

    def compare(self, tree1: BaseTree, tree2: BaseTree):
        term_names1 = [term.name for term in tree1.get_terminals()]
        term_names2 = [term.name for term in tree2.get_terminals()]
        # false if terminals or BitStrings are not the same
        return set(term_names1) == set(term_names2) and self._bitstrs(tree1) == self._bitstrs(tree2)

    @pytest.mark.parametrize("format_list", [(["fasta", "phylip"])])
    def test_sprt_infile_formats(self, format_list: list):
        """
        Test function to check that RAxML module supports the required infile formats.

        Arguments :
            format_list: List of accepted infile formats for using RAxML tool.
        """
        assert len(Rax.SPRT_INFILE_FORMATS) == len(format_list)
        assert Rax.SPRT_INFILE_FORMATS == format_list

    @pytest.mark.parametrize(
        "keywords",
        [
            (
                {
                    "default": ["-m", "GTRCAT", "--silent"],
                    "JC+CAT": ["-m", "GTRCAT", "--JC69", "--silent"],
                    "JC+G": ["-m", "GTRGAMMA", "--JC69", "--silent"],
                    "K80+CAT": ["-m", "GTRCAT", "--K80", "--silent"],
                    "K80+G": ["-m", "GTRGAMMA", "--K80", "--silent"],
                    "HKY+CAT": ["-m", "GTRCAT", "--HKY85", "--silent"],
                    "HKY+G": ["-m", "GTRGAMMA", "--HKY85", "--silent"],
                    "GTR+CAT": ["-m", "GTRCAT", "--silent"],
                    "GTR+G": ["-m", "GTRGAMMA", "--silent"],
                    "JTT+CAT": ["-m", "PROTCATJTT", "--silent"],
                    "JTT+G": ["-m", "PROTGAMMAJTT", "--silent"],
                    "WAG+CAT": ["-m", "PROTCATWAG", "--silent"],
                    "WAG+G": ["-m", "PROTGAMMAWAG", "--silent"],
                }
            ),
        ],
    )
    def test_keywords(self, keywords: dict):
        """
        Test function to check that RAxML module has the same keys and values for the keywords dictionary,
        that will be used later to pass arguments to a command.

        Arguments :
            keywords: Dictionary where keys are shortcuts for arguments combinations and values are lists of
            arguments a command may receive.
        """
        assert len(Rax.KEYWORDS) == len(keywords)
        for aKey, bKey in zip(Rax.KEYWORDS.keys(), keywords.keys()):
            assert aKey == bKey
        for key in Rax.KEYWORDS:
            assert Rax.KEYWORDS[key] == keywords[key]

    @pytest.mark.parametrize(
        "args, infile_path, bootstraps, seed, arg_list",
        [
            (
                "JC+CAT",
                "tests/flatfiles/f001.mafft_linsi.aln",
                1,
                404,
                [
                    "-p",
                    "404",
                    "-n",
                    str(os.getpid()),
                    "-T",
                    str(NUMCORES),
                    "-w",
                    tmp_dir,
                    "-m",
                    "GTRCAT",
                    "--JC69",
                    "--silent",
                    "-N",
                    "1",
                    "-s",
                    "tests/flatfiles/f001.mafft_linsi.aln",
                ],
            ),
            (
                "-p 648 -m GTRGAMMA --HKY85",
                "tests/flatfiles/f001.mafft_linsi.aln",
                0,
                648,
                [
                    "-p",
                    "648",
                    "-m",
                    "GTRGAMMA",
                    "--HKY85",
                    "-s",
                    "tests/flatfiles/f001.mafft_linsi.aln",
                ],
            ),
        ],
    )
    def test_gen_args(self, args: str, infile_path: str, bootstraps: int, seed: int, arg_list: list):
        """
        Test function to ensure the construction of an argument list taking into account the parameters a
        command might need, and how it needs them. Ultimately, it ensures the correct execution of RAxML
        gen_args's function.

        Arguments :
            args: Keyword or arguments to use in the call of RAxML, excluding
                infile and outfile arguments.
            infile_path: Input alignment file path.
            bootstraps: Number of bootstraps to generate.
            tmpdir_path: Path of the directory we want to save the RAxML output data into (Only required while testing for
                reproducibility purposes or just if the user wants a specific folder to allocate the result of the execution
                into. Otherwise, a random one will be provided).
            seed: Specify a random number seed for the parsimony inferences (Only required for while testing for
                reproducibility purposes. Otherwise, a random one will be provided).
            arg_list: Handmade constructed argument list as the expected result the real function should return.
        """
        assert arg_list == Rax.gen_args(args, infile_path, bootstraps, self.tmp_dir, seed)

    @pytest.mark.parametrize(
        "command, score, infile_path, expected_inference_tree, expected_inference_info",
        [
            (
                [
                    "raxmlHPC",
                    "-p",
                    "404",
                    "-T",
                    str(NUMCORES),
                    "-m",
                    "GTRCAT",
                    "-N",
                    "1",
                ],
                -1975.511424,
                Path("tests/flatfiles/f001.mafft_linsi.aln"),
                Path("f001.mafft_linsi.RAxML_best_tree_1"),
                Path("f001.mafft_linsi.RAxML_info_1"),
            ),
            (
                [
                    "raxmlHPC",
                    "-p",
                    "404",
                    "-T",
                    str(NUMCORES),
                    "-m",
                    "GTRGAMMA",
                    "--HKY85",
                    "--silent",
                ],
                -1988.653009,
                Path("tests/flatfiles/f001.mafft_linsi.aln"),
                Path("f001.mafft_linsi.RAxML_best_tree_2"),
                Path("f001.mafft_linsi.RAxML_info_2"),
            ),
        ],
    )
    def test_get_results(
        self,
        command: list,
        score: float,
        infile_path: Path,
        expected_inference_tree: Path,
        expected_inference_info: Path,
    ):
        """
        Test function to ensure the correct extraction of the phylogeny and it's associated score from a command and the
        resultant output.  Ultimately, it ensures the correct execution of FastTree's get_results function.

        Arguments :
         command: FastTree's command line executed.
         score: The associated score an inferenced phylogeny tree has.
         infile_path:  Input alignment file path.
         expected_inference_tree: Path where the Phylo.BaseTree output of subprocess.run(command) is stored,
         to avoid the unnecessary execution of such an expensive function.
        expected_inference_info: Path where the information output (time, alignment patterns, score...) of subprocess.run(command)
            is stored, to avoid the unnecessary execution of such an expensive function.
        """
        # random temporary file path generation to save the results of the execution into
        r = random.randint(1, 999999)

        treefile_path = Path.joinpath(self.tmp_dir, "RAxML_bestTree." + str(r))
        infofile_path = Path.joinpath(self.tmp_dir, "RAxML_info." + str(r))

        command += ["-n", str(r), "-w", self.tmp_dir, "-s", infile_path.absolute()]

        expected_inference_tree = pytest.data_dir / expected_inference_tree
        expected_inference_info = pytest.data_dir / expected_inference_info

        run_mocker = MockStdOut(
            expected_inference_tree, treefile_path, expected_inference_info, infofile_path
        )
        mocked_subprocess_info = run_mocker.get_mocked_info_output()

        phylogeny = Phylo.read(treefile_path, "newick")
        res_tree, res_score = Rax.get_results(command, mocked_subprocess_info)

        assert score == res_score
        assert self.compare(phylogeny, res_tree)

        """
        Another testing way by reading the tree and comparing both strings, to avoid running the 
        command. However, it does not work, as format ends up rounding decimals, making the subsequent
        comparison to fail.
        
        str_tree=""
        with open(treefile_path, "r") as treefile:
            content = treefile.readlines()
            for line in content:
                str_tree += line
        phylogeny = Phylo.read(treefile_path, "newick")
        res_tree, res_score = Rax.get_results(command, output)
        assert score == res_score
        assert str_tree == res_tree.format("newick")
       """

   @pytest.mark.parametrize(
        "command, tmp_file",
        [
            (["-w", tmp_dir], None),
            (["-w", tmp_dir], tmp_dir),
        ],
    )
    def test_cleanup(self, command: list, tmp_file: str):
        """
        Test function to ensure the correct removal of all the files of a given directory.
        Arguments :
            command: RAxML's command line executed.
            tmp_file: Path of the folder we want to delete (Just for testing purposes,
            because when called from __init__.py, the cleanup input is a bit different).

        """
        Rax.cleanup(command, tmp_file)
        assert not self.tmp_dir.exists()
