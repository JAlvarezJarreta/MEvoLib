from pathlib import Path
import subprocess
import pytest
import random
import os

from Bio import Phylo
from Bio.Phylo.Consensus import _BitString

from mevolib.inference import _RAxML as Rax
from mevolib._utils import NUMCORES


class MockStdOut:
    def __init__(self, tree_infile_path, tree_outfile_path, info_infile_path, info_outfile_path):
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
    tmp_dir: Path = Path("tests/raxml_tmp_dir/").absolute()
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
        "args, infile_path, bootstraps, tmpdir_path, seed, arg_list",
        [
            (
                "JC+CAT",
                "tests/Fasta/f001.mafft_linsi.aln",
                1,
                tmp_dir,
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
                    "tests/Fasta/f001.mafft_linsi.aln",
                ],
            ),
            (
                "-p 648 -m GTRGAMMA --HKY85",
                "tests/Fasta/f001.mafft_linsi.aln",
                0,
                None,
                648,
                [
                    "-p",
                    "648",
                    "-m",
                    "GTRGAMMA",
                    "--HKY85",
                    "-s",
                    "tests/Fasta/f001.mafft_linsi.aln",
                ],
            ),
        ],
    )
    def test_gen_args(
        self, args: str, infile_path: str, bootstraps: int, tmpdir_path: str, seed: int, arg_list: list
    ):
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
        assert arg_list == Rax.gen_args(args, infile_path, bootstraps, tmpdir_path, seed)

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
                -1974.894224,
                Path("tests/Fasta/f001.mafft_linsi.aln"),
                Path("f001.mafft_linsi.RAxML_best_tree_1.aln"),
                Path("f001.mafft_linsi.RAxML_info_1.aln"),
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
                Path("tests/Fasta/f001.mafft_linsi.aln"),
                Path("f001.mafft_linsi.RAxML_best_tree_2.aln"),
                Path("f001.mafft_linsi.RAxML_info_2.aln"),
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
        resultant output.  Ultimately, it ensures the correct execution of FastTree get_results's function.
        Arguments :
         command: FastTree's command line executed.
         score: The associated score an inferenced phylogeny tree has.
         infile_path:  Input alignment file path.
         expected_inference_tree: Path where the Phylo.BaseTree output of subprocess.run(command) is stored, to avoid the unnecesary
            execution of such an expensive function.
        expected_inference_info: Path where the information output (time, alignment patterns, score...) of subprocess.run(command)
            is stored, to avoid the unnecesary execution of such an expensive function.
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
        mocked_subprocess_tree = run_mocker.get_mocked_tree_output()
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
