# Chem Func

[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/chemfunc)](https://badge.fury.io/py/chemfunc)
[![PyPI version](https://badge.fury.io/py/chemfunc.svg)](https://badge.fury.io/py/chemfunc)
[![Downloads](https://pepy.tech/badge/chemfunc)](https://pepy.tech/project/chemfunc)
[![license](https://img.shields.io/github/license/swansonk14/chemfunc.svg)](https://github.com/swansonk14/chemfunc/blob/main/LICENSE.txt)

Useful functions and scripts for working with small molecules.

## Installation

Optionally, create a conda environment.
```bash
conda create -y -n chemfunc python=3.10
conda activate chemfunc
```

Install the latest version of Chem Func using pip.
```
pip install chemfunc
```

Alternatively, clone the repository and install the local version of the package.
```
git clone https://github.com/swansonk14/chemfunc.git
cd chemfunc
pip install -e .
```

If there are version issues with the required packages, create a conda environment with specific working versions of the packages as follows.
```bash
pip install -r requirements.txt
pip install -e .
```

**Note:** If you get the issue `ImportError: libXrender.so.1: cannot open shared object file: No such file or directory`, run `conda install -c conda-forge xorg-libxrender`.


## Features

Chem Func contains a variety of useful functions and scripts for working with small molecules.

Functions can be imported from the `chemfunc` package. For example:
```python
from pathlib import Path
from chemfunc.sdf_to_smiles import sdf_to_smiles

sdf_to_smiles(
    data_path=Path('molecules.sdf'),
    save_path=Path('molecules.csv')
)
```

Most modules can also be run as scripts from the command line using the `chemfunc` command along with the appropriate function name. For example:
```bash
chemfunc sdf_to_smiles \
    --data_path molecules.sdf \
    --save_path molecules.csv
```

To see a list of available scripts, run `chemfunc -h`.

For each script, run `chemfunc <script_name> -h` to see a description of the arguments for that script.


## Contents

Below is a list of the contents of the package.

[`canonicalize_smiles.py`](https://github.com/swansonk14/chemfunc/blob/main/chemfunc/canonicalize_smiles.py) (function, script)

Canonicalizes SMILES using RDKit canonicalization and optionally strips salts.

[`chemical_diversity.py`](https://github.com/swansonk14/chemfunc/blob/main/chemfunc/chemical_diversity.py) (function, script)

Computes the chemical diversity of a set of molecules in terms of Tanimoto distances.

[`cluster_molecules.py`](https://github.com/swansonk14/chemfunc/blob/main/chemfunc/cluster_molecules.py) (function, script)

Performs k-means clustering to cluster molecules based on Morgan fingerprints.

[`compute_property_distribution.py`](https://github.com/swansonk14/chemfunc/blob/main/chemfunc/compute_property_distribution.py) (function, script)

Computes one or more molecular properties for a set of molecules.

[`deduplicate_smiles.py`](https://github.com/swansonk14/chemfunc/blob/main/chemfunc/deduplicate_smiles.py) (function, script)

Deduplicate a CSV files by SMILES.

[`filter_molecules.py`](https://github.com/swansonk14/chemfunc/blob/main/chemfunc/filter_molecules.py) (function, script)

Filters molecules to those with values in a certain range.

[`measure_experimental_reproducibility.py`](https://github.com/swansonk14/chemfunc/blob/main/chemfunc/measure_experimental_reproducibility.py) (function, script)

Measures the experimental reproducibility of two biological replicates by using one replicate to predict the other.

[`molecular_fingerprints.py`](https://github.com/swansonk14/chemfunc/blob/main/chemfunc/molecular_fingerprints.py) (functions, script)

Contains functions to compute fingerprints for molecules. Parallelized for speed. The function `save_fingerprints` can be used as a script to compute fingerprints from a CSV file and save them as an NPZ file.

[`molecular_properties.py`](https://github.com/swansonk14/chemfunc/blob/main/chemfunc/molecular_properties.py) (functions)

Contains functions to compute molecular properties. Parallelized for speed.

[`molecular_similarities.py`](https://github.com/swansonk14/chemfunc/blob/main/chemfunc/molecular_similarities.py) (functions)

Contains functions to compute similarities between molecules. Parallelized for speed.

[`nearest_neighbor.py`](https://github.com/swansonk14/chemfunc/blob/main/chemfunc/nearest_neighbor.py) (function, script)

Given a dataset of molecules, computes the nearest neighbor molecule in a second dataset using one of several similarity metrics.

[`plot_property_distribution.py`](https://github.com/swansonk14/chemfunc/blob/main/chemfunc/plot_property_distribution.py) (function, script)

Plots the distribution of molecular properties of a set of molecules.

[`plot_tsne.py`](https://github.com/swansonk14/chemfunc/blob/main/chemfunc/plot_tsne.py) (function, script)

Runs a t-SNE on molecular fingerprints from one or more chemical libraries.

[`regression_to_classification.py`](https://github.com/swansonk14/chemfunc/blob/main/chemfunc/regression_to_classification.py) (function, script)

Converts regression data to classification data using given thresholds.

[`sample_molecules.py`](https://github.com/swansonk14/chemfunc/blob/main/chemfunc/sample_molecules.py) (function, script)

Samples molecules from a CSV file, either uniformly at random across the entire dataset or uniformly at random from each cluster within the data.

[`sdf_to_smiles.py`](https://github.com/swansonk14/chemfunc/blob/main/chemfunc/sdf_to_smiles.py) (function, script)

Converts an SDF file to a CSV file with SMILES.

[`select_from_clusters.py`](https://github.com/swansonk14/chemfunc/blob/main/chemfunc/select_from_clusters.py) (function, script)

Selects the best molecule from each cluster.

[`smiles_to_svg.py`](https://github.com/swansonk14/chemfunc/blob/main/chemfunc/smiles_to_svg.py) (function, script)

Converts a SMILES string to an SVG image of the molecule.

[`visualize_molecules.py`](https://github.com/swansonk14/chemfunc/blob/main/chemfunc/visualize_molecules.py)(function, script)

Converts a file of SMILES to images of molecular structures.

[`visualize_reactions.py`](https://github.com/swansonk14/chemfunc/blob/main/chemfunc/visualize_reactions.py) (function, script)

Converts a file of reaction SMARTS to images of chemical reactions.
