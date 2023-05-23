# Chem Utils

Useful functions and scripts for working with small molecules.

## Installation

Install the latest version of Chem Utils using pip.
```
pip install chem_utils
```

Alternatively, clone the repository and install the local version of the package.
```
git clone https://github.com/swansonk14/chem_utils.git
cd chem_utils
pip install -e .
```


## Features

Chem Utils contains a variety of useful functions and scripts for working with small molecules.

Functions can be imported from the `chem_utils` package. For example:
```python
from pathlib import Path
from chem_utils.sdf_to_smiles import sdf_to_smiles

sdf_to_smiles(
    data_path=Path('molecules.sdf'),
    save_path=Path('molecules.csv')
)
```

Most modules can also be run as scripts from the command line using the `chem_utils` command along with the appropriate function name. For example:
```bash
chem_utils sdf_to_smiles \
    --data_path molecules.sdf \
    --save_path molecules.csv
```

To see a list of available scripts, run `chem_utils -h`.

For each script, run `chem_utils <script_name> -h` to see a description of the arguments for that script.


## Contents

Below is a list of the contents of the package.

[`canonicalize_smiles.py`](chem_utils/canonicalize_smiles.py) (function, script)

Canonicalizes SMILES using RDKit canonicalization and optionally strips salts.

[`chemical_diversity.py`](chem_utils/chemical_diversity.py) (function, script)

Computes the chemical diversity of a set of molecules in terms of Tanimoto distances.

[`cluster_molecules.py`](chem_utils/cluster_molecules.py) (function, script)

Performs k-means clustering to cluster molecules based on Morgan fingerprints.

[`compute_property_distribution.py`](chem_utils/compute_property_distribution.py) (function, script)

Computes one or more molecular properties for a set of molecules.

[`deduplicate_smiles.py`](chem_utils/deduplicate_smiles.py) (function, script)

Deduplicate a CSV files by SMILES.

[`filter_molecules.py`](chem_utils/filter_molecules.py) (function, script)

Filters molecules to those with values in a certain range.

[`measure_experimental_reproducibility.py`](chem_utils/measure_experimental_reproducibility.py) (function, script)

Measures the experimental reproducibility of two biological replicates by using one replicate to predict the other.

[`molecular_fingerprints.py`](chem_utils/molecular_fingerprints.py) (functions)

Contains functions to compute fingerprints for molecules. Parallelized for speed.

[`molecular_properties.py`](chem_utils/molecular_properties.py) (functions)

Contains functions to compute molecular properties. Parallelized for speed.

[`molecular_similarities.py`](chem_utils/molecular_similarities.py) (functions)

Contains functions to compute similarities between molecules. Parallelized for speed.

[`nearest_neighbor.py`](chem_utils/nearest_neighbor.py) (function, script)

Given a dataset of molecules, computes the nearest neighbor molecule in a second dataset using one of several similarity metrics.

[`plot_property_distribution.py`](chem_utils/plot_property_distribution.py) (function, script)

Plots the distribution of molecular properties of a set of molecules.

[`plot_tsne.py`](chem_utils/plot_tsne.py) (function, script)

Runs a t-SNE on molecular fingerprints from one or more chemical libraries.

[`regression_to_classification.py`](chem_utils/regression_to_classification.py) (function, script)

Converts regression data to classification data using given thresholds.

[`sample_molecules.py`](chem_utils/sample_molecules.py) (function, script)

Samples molecules from a CSV file, either uniformly at random across the entire dataset or uniformly at random from each cluster within the data.

[`sdf_to_smiles.py`](chem_utils/sdf_to_smiles.py) (function, script)

Converts an SDF file to a CSV file with SMILES.

[`select_from_clusters.py`](chem_utils/select_from_clusters.py) (function, script)

Selects the best molecule from each cluster.

[`smiles_to_svg.py`](chem_utils/smiles_to_svg.py) (function, script)

Converts a SMILES string to an SVG image of the molecule.

[`visualize_molecules.py`](chem_utils/visualize_molecules.py)(function, script)

Converts a file of SMILES to images of molecular structures.

[`visualize_reactions.py`](chem_utils/visualize_reactions.py) (function, script)

Converts a file of reaction SMARTS to images of chemical reactions.
