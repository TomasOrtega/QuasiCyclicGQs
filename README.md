# Quasi Cyclic Representation of Classical Generalized Quadrangles

This repository contains the Quasi Cyclic (QC) representation of classical [Generalized Quadrangles](https://en.wikipedia.org/wiki/Generalized_quadrangle) (GQ). 
These have been obtained using the methods described in "Practical Implementation of Geometric Quasi-Cyclic LDPC Codes" by Simeon Michael Ball and Tomàs Ortega.
See also [PCT/EP2023/062797](https://patentscope.wipo.int/search/en/detail.jsf?docId=WO2023218050&_cid=P11-LP8MC2-94041-1).

Standard [.alist](https://www.inference.org.uk/mackay/codes/alist.html) files are provided for LDPC simulations.

## Overview

[Generalized Quadrangles](https://en.wikipedia.org/wiki/Generalized_quadrangle) are incidence structures whose main feature is the lack of any triangles (yet they contain many quadrangles). 
This repository provides a practical Quasi Cyclic representation of their point-line incidence matrix.
These matrices can be seen as parity check matrices of error correcting codes, which are particularly useful for LDPC (Low-Density Parity-Check) applications.

## Contents

All representations are in the **[representations](representations)** folder, with subdirectories

- **[elliptic_quadrangle](representations/elliptic_quadrangle)** - representations for $Q(5,q)$ and its dual.
- **[hermitian_quadrangle](representations/hermitian_quadrangle)** - representations for $W(3,q)$ and its dual.
- **[symplectic_quadrangle](representations/symplectic_quadrangle)** - representations for $H(4,q^2)$ and its dual.
- **[others](representations/others)** - other representations, such as for the projective space $PG(k-1,q)$.

All folders contain the raw output from our generation scripts in .txt files, as well as [.alist](https://www.inference.org.uk/mackay/codes/alist.html) files for LDPC simulations. Files that start with "G_" contain the generator matrix for the corresponding code described by the parity check matrix.

## Usage

The [.alist](https://www.inference.org.uk/mackay/codes/alist.html) files provided are ready to use with any standard LDPC BER/FER simulator, such as [aff3ct](https://aff3ct.github.io/index.html).

The [utilities](utilities) folder contains Python scripts that convert raw outputs to .alist files. The input and output filenames for the conversion have to be edited in the code.

## Citation

If you find this work helpful, please consider citing

    @inproceedings{ball2024QC,
      title={Practical Implementation of Geometric Quasi-Cyclic LDPC Codes},
      author={Ball, Simeon and Ortega, Tomàs},
      booktitle={Discrete Mathematics Days},
      year={2024}
    }

## Contact

For questions or inquiries about this repository, feel free to email [tomas.ortega@uci.edu](mailto:tomas.ortega@uci.edu).

