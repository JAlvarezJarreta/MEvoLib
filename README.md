
# **M**olecular **Evo**lution **Lib**rary for Python

The MEvoLib is a library of freely available Python tools for molecular evolution.

The documentation file with further details on how to use our library is included in the release in "docs/manual.pdf".

This MEvoLib package is open source software made available under the Apache 2.0 License terms. Please see the LICENSE file for further details.

If you use MEvoLib in work contributing to a scientific publication, you can
refer to the latest paper:

> Alvarez-Jarreta, J., Ruiz-Pesini, E. [_MEvoLib v1.0: the first molecular evolution library for Python._](https://doi.org/10.1186/s12859-016-1303-3) BMC Bioinformatics 17, 436 (2016).


## Installation

MEvoLib is currently supported and tested for Python 3.8+.

To build, test and install MEvoLib, download and decompress the source code, go
to the source directory, and type:

    python setup.py build
    python setup.py test
    sudo python setup.py install

If you need to do an additional configuration, e.g. changing the base directory,
type:

    python setup.py


## Additional requirements

Many modules require third party tools you should install, such as Mafft, PhyML or SuperFine.
