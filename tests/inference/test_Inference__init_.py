from contextlib import nullcontext as does_not_raise
from io import StringIO
import os
from pathlib import Path
import random
from typing import ContextManager
import subprocess

from Bio import Phylo
from Bio.Phylo import BaseTree
from Bio.Phylo.Consensus import _BitString
import pytest
from pytest import raises

from mevolib import inference
from mevolib.inference import _FastTree as Fast, _RAxML as Rax
from mevolib._utils import NUMCORES


class MockStdOut:
    """
    Class defined to avoid making subprocess.run(command) every time we want to test the test_get_results function.
    The process depends on the tool we are using, as FastTree gives all the output in a whole file, whilst RAxML(HPC)
    gives it splitted into some files; so in the second case we have to manage more documents.

    Arguments:
        tree_infile_path: Path where the resultant tree of the command execution is stored into.
        tree_outfile_path: Path where the resultant tree of the command execution and its score should be copied into; as the
        Inference module will need it to extract info from it.
        score: Associated phylogenetic score to the mentioned tree.
        info_infile_path: Path where the RAxML information of the resultant tree of the command execution is stored into.
        info_outfile_path: Path where the RAxML information of the resultant tree of the command execution should be copied
        into; as the Inference module will need it to extract info from it.
    """

    def __init__(self):
        pass

    # FastTree mocker file management
    def fastTreeAction(self, tree_infile_path: Path, tree_outfile_path: Path, score: float):
        with open(tree_infile_path, "r") as fin:
            self.stdout = fin.read()
        with open(tree_outfile_path, "w") as fin:
            fin.write(self.stdout)
            fin.write(f"TreeLogLk	Placeholder	{score}")

    def get_mocked_output(self):
        return self.stdout

    # RAxMl mocker file management
    def raxmlAction(
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


class TestInferenceInit:
    """
    Class made to ensure the correct operating of MEvoLib Inference module' functions.
    """

    tmp_dir: Path = Path("tests/inference_init_tmp_dir/").absolute()
    if not tmp_dir.exists():
        os.mkdir(tmp_dir)

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

    @pytest.mark.parametrize(
        "phylo, boots",
        [
            ({"fasttree": inference._FastTree, "raxml": inference._RAxML}, {}),
        ],
    )
    def test_phylo_to_lib(self, phylo: dict, boots: dict):
        """
        Test function to check that the manin Inference module supports the required infile formats.

        Arguments :
            phylo: Dictionary that contains the name of a phylogenetic tool as key and it's matching
            MEvoLib module as the value.
            Bootstrapping: Dictionary that contains MEvoLib's bootstrapping tools.
        """
        dicti = dict([("inference", list(phylo.keys())), ("bootstrap", list(boots.keys()))])
        assert dicti == inference.get_tools()

    @pytest.mark.parametrize(
        "tool, expected",
        [
            (
                "fasttree",
                does_not_raise(),
            ),
            (
                "raxml",
                does_not_raise(),
            ),
            (
                "motorola",
                raises(ValueError),
            ),
        ],
    )
    def test_get_keywords(self, tool: str, expected: ContextManager):
        """
        Test function to check that the main Inference module has the same tools in the dictionary.

        Arguments :
            tool: String that represents the name of a Inference tool.
            expected: Context Manager that indicates whether or not an exception should be raised.
        Raises:
            ValueError: if the tool is not included neither in the phylogenetic or the bootstrapping Inference's
            tools.
        """
        tool = tool.lower()
        keyword_dict = {}
        if tool in inference._PHYLO_TOOL_TO_LIB:
            tool_lib_dict = inference._PHYLO_TOOL_TO_LIB
        else:  # tool in _BOOTS_TOOL_TO_LIB
            tool_lib_dict = inference._BOOTS_TOOL_TO_LIB
        for key, value in tool_lib_dict.items():
            keyword_dict[key] = value
        with expected:
            assert keyword_dict == inference.get_keywords(tool)

    @pytest.mark.parametrize(
        "binary, infile, expected_output_fast, expected_inference_tree_rax, "
        + "expected_inference_info_rax, infile_format, args, outfile, "
        + "outfile_format,bootstraps,  seed, command, "
        + "score, expected",
        [
            (
                "FastTree",
                "tests/flatfiles/f001.mafft_default.aln",
                Path("f001.mafft_default.FastTree_output_1"),
                None,
                None,
                "fasta",
                "GTR+CAT",
                None,
                "newick",
                0,
                None,
                [
                    "fasttree",
                    "-gtr",
                    "-nt",
                    "-nopr",
                    "-quiet",
                ],
                -1911.868,
                does_not_raise(),
            ),
            (
                "RAxML",
                "tests/flatfiles/f001.mafft_linsi.aln",
                None,
                Path("f001.mafft_linsi.RAxML_best_tree_1"),
                Path("f001.mafft_linsi.RAxML_info_1"),
                "fasta",
                "default",
                None,
                "newick",
                1,
                404,
                [
                    "raxmlHPC",
                    "-p",
                    "404",
                    "-T",
                    str(NUMCORES),
                    "-m",
                    "GTRCAT",
                    "--silent",
                    "-N",
                    "1",
                ],
                -1974.894224,
                does_not_raise(),
            ),
            (
                "RAxMLh",
                "tests/flatfiles/f001.mafft_linsi.aln",
                None,
                Path("f001.mafft_linsi.RAxML_best_tree_1"),
                Path("f001.mafft_linsi.RAxML_info_1"),
                "fasta",
                "default",
                None,
                "newick",
                1,
                404,
                [
                    "raxmlHPC",
                    "-p",
                    "404",
                    "-T",
                    str(NUMCORES),
                    "-m",
                    "GTRCAT",
                    "--silent",
                    "-N",
                    "1",
                ],
                -1974.894224,
                pytest.raises(ValueError),
            ),
            (
                "FastTree",
                "tests/flatfiles/output_aln.stockholm",
                Path("f001.mafft_default.FastTree_output_1"),
                None,
                None,
                "stockholm",
                "GTR+CAT",
                "output",
                "newick",
                0,
                None,
                [
                    "fasttree",
                    "-gtr",
                    "-nt",
                    "-nopr",
                    "-quiet",
                ],
                -1911.868,
                does_not_raise(),
            ),
        ],
    )
    def test_get_phylogeny(
        self,
        binary: str,
        infile: str,
        expected_output_fast: Path,
        expected_inference_tree_rax,
        expected_inference_info_rax,
        infile_format: str,
        args: str,
        outfile: str,
        outfile_format: str,
        bootstraps: int,
        seed: int,
        command: list,
        score: float,
        expected: ContextManager,
    ):
        """
        Test function to ensure the correct extraction of the phylogeny and it's associated score from a command and the
        resultant output.  Ultimately, it ensures the correct execution of Inference's get_phylogeny function.

        Arguments:
            binary: String that represents the Inference tool that wants to be used.
            infile: Path of the alignment input file that we want to do the analysis from.
            expected_output_fast: In case we are using FastTree, path of the output that the execution of the command
            would get; and None otherwise.
            expected_inference_tree_rax: In case we are using RAxML, path of the output tree that the execution of the
            command would get; and None otherwise.
            expected_inference_info_rax: In case we are using RAxML, path where the information output (time,
            alignment patterns, score...) that the execution of the command would get; and None otherwise.
            infile_format: Extension of the input file.
            args: Keyword or arguments to use in the call of the phylogenetic
            inference tool, excluding infile and outfile arguments. By default,
            'default' arguments are used.
            outfile: Phylogenetic tree output file.
            outfile_format: Output file format. By default, NEWICK format.
            bootstraps: Number of bootstraps to generate. By default, 0 (only use the
            input alignment).
            seed: Specify a random number seed for the parsimony inferences (Only required for while testing for
            reproducibility purposes. Otherwise, a random one will be provided).
            command: List of parameters that compounds the whole command to do a custom Inference.
            score: Phylogenetic score associated to the inference of the input_file using the selected tool.
            expected: Context Manager that indicates whether or not an exception should be raised.
        """
        if not self.tmp_dir.exists():
            os.mkdir(self.tmp_dir)
        with expected:
            run_mocker = MockStdOut()

            if not (outfile is None):
                outfile = Path.joinpath(self.tmp_dir, outfile)

            r = random.randint(1, 999999)
            infile_path = Path(infile).absolute()

            if command[0] == "fasttree":
                treefile_path = Path.joinpath(self.tmp_dir, "FastTree_bestTree." + str(r)).absolute()
                command += ["-log", treefile_path, infile_path]
                tmpfile = treefile_path

                # FastTree command mocking management:
                expected_output_fast = pytest.data_dir / expected_output_fast

                run_mocker.fastTreeAction(expected_output_fast, treefile_path, score)
                mocked_subprocess_output = run_mocker.get_mocked_output()

                phylogeny = Phylo.read(StringIO(mocked_subprocess_output), "newick")

            else:
                command += ["-n", str(r), "-w", self.tmp_dir, "-s", infile_path]
                tmpfile = self.tmp_dir

                # RAxML command mocking management:
                treefile_path = Path.joinpath(self.tmp_dir, "RAxML_bestTree." + str(r)).absolute()
                infofile_path = Path.joinpath(self.tmp_dir, "RAxML_info." + str(r)).absolute()

                expected_inference_tree_rax = pytest.data_dir / expected_inference_tree_rax
                expected_inference_info_rax = pytest.data_dir / expected_inference_info_rax

                run_mocker.raxmlAction(
                    expected_inference_tree_rax, treefile_path, expected_inference_info_rax, infofile_path
                )

                mocked_subprocess_tree = run_mocker.get_mocked_tree_output()

                phylogeny = Phylo.read(StringIO(mocked_subprocess_tree), "newick")

            res_tree, res_score = inference.get_phylogeny(
                binary,
                infile,
                infile_format,
                args,
                outfile,
                outfile_format,
                bootstraps,
                tmpfile,
                self.tmp_dir,
                seed,
            )

            assert self.compare(phylogeny, res_tree)
            assert score == res_score
