from contextlib import nullcontext as does_not_raise
from io import StringIO
from pathlib import Path
import random
from typing import ContextManager
import subprocess
import pytest
from pytest import raises

from Bio import Phylo
from Bio.Phylo import BaseTree
from Bio.Phylo.Consensus import _BitString

from mevolib import inference
from mevolib.inference import _FastTree as Fast, _RAxML as Rax
from mevolib._utils import NUMCORES


class TestInferenceInit:
    tmp_dir: Path = Path("tests/inference_init_tmp_dir/").absolute()
    manifest_dir: Path

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

    @pytest.mark.parametrize(
        "phylo, boots",
        [
            ({"fasttree": inference._FastTree, "raxml": inference._RAxML}, {}),
        ],
    )
    def test_phylo_to_lib(self, phylo: dict, boots: dict):
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
        "binary, infile, infile_format,args, outfile, "
        + "outfile_format,bootstraps, tmpfile, seed, command, "
        + "score, expected",
        [
            (
                "FastTree",
                "tests/Fasta/f001.mafft_default.aln",
                "fasta",
                "GTR+CAT",
                None,
                None,
                0,
                tmp_dir,
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
                "tests/Fasta/f001.mafft_linsi.aln",
                "fasta",
                "default",
                None,
                None,
                1,
                tmp_dir,
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
                "tests/Fasta/f001.mafft_linsi.aln",
                "fasta",
                "default",
                None,
                None,
                1,
                tmp_dir,
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
        ],
    )
    def test_get_phylogeny(
        self,
        binary: str,
        infile: str,
        infile_format: str,
        args: str,
        outfile: str,
        outfile_format: str,
        bootstraps: int,
        tmpfile: str,
        seed: int,
        command: list,
        score: float,
        expected: ContextManager,
    ):
        with expected:
            r = random.randint(1, 999999)
            infile_path = Path(infile).absolute()
            if command[0] == "fasttree":
                treefile_path = Path.joinpath(self.tmp_dir, "FastTree_bestTree." + str(r)).absolute()
                command += ["-log", treefile_path, infile_path]
                tmpfile = treefile_path
            else:
                treefile_path = Path.joinpath(self.tmp_dir, "RAxML_bestTree." + str(r)).absolute()
                command += ["-n", str(r), "-w", self.tmp_dir, "-s", infile_path]
                tmpfile = self.tmp_dir
            output = subprocess.run(
                command,
                stderr=subprocess.DEVNULL,
                universal_newlines=True,
                check=True,
                stdout=subprocess.PIPE,
            ).stdout
            if command[0] == "fasttree":
                phylogeny = Phylo.read(StringIO(output), "newick")
            else:
                phylogeny = Phylo.read(treefile_path, "newick")

            res_tree, res_score = inference.get_phylogeny(
                binary, infile, infile_format, args, outfile, outfile_format, bootstraps, tmpfile, seed
            )

            assert self.compare(phylogeny, res_tree)
            assert score == res_score
