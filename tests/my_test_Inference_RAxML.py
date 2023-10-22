from pathlib import Path
import subprocess
import pytest
import random
import os


from Bio import Phylo
from Bio.Phylo.Consensus import _BitString

from mevolib.inference import _RAxML as Rax
from mevolib._utils import NUMCORES


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
        "command, score, infile_path",
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
        treefile_path = Path.joinpath(self.tmp_dir, "RAxML_bestTree." + str(r))

        command += ["-n", str(r), "-w", self.tmp_dir, "-s", infile_path.absolute()]
        output = subprocess.run(
            command, stderr=subprocess.DEVNULL, universal_newlines=True, check=True, stdout=subprocess.PIPE
        ).stdout
        phylogeny = Phylo.read(treefile_path, "newick")
        res_tree, res_score = Rax.get_results(command, output)

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
