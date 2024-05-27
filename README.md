# **M**olecular **Evo**lution **Lib**rary for Python

[![CI](https://github.com/JAlvarezJarreta/MEvoLib/actions/workflows/python-ci.yml/badge.svg?branch=main)](https://github.com/JAlvarezJarreta/MEvoLib/actions/workflows/python-ci.yml)
![Coverage](https://gist.githubusercontent.com/JAlvarezJarreta/0c3fa1e9db6f5dcbe69b6b4a3f4a0501/raw/badge.svg)

The MEvoLib is a library of freely available Python tools for molecular evolution.

The documentation file with further details on how to use our library is included in the release in "docs/manual.pdf".

This MEvoLib package is open source software made available under the Apache 2.0 License terms. Please see the LICENSE file for further details.

If you use MEvoLib in work contributing to a scientific publication, you can
refer to the latest paper:

> Alvarez-Jarreta, J., Ruiz-Pesini, E. [_MEvoLib v1.0: the first molecular evolution library for Python._](https://doi.org/10.1186/s12859-016-1303-3) BMC Bioinformatics 17, 436 (2016).


## Installation

MEvoLib is currently supported and tested for Python 3.8+.

Installing MEvoLib is really easy and straightforward:
```bash
git clone --branch version/2.0 https://github.com/JAlvarezJarreta/MEvoLib
pip install MEvoLib/.
```

If you want to install in development mode, you will need to include the `dev` tag:
```bash
pip install -e MEvoLib/.[dev]
```

## Additional requirements

Many modules require third party tools you should install, such as Mafft, PhyML or SuperFine.
