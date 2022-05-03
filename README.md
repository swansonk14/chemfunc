# Chem Utils

Useful scripts and functions for working with small molecules.

## Installation

Install conda environment.
```
conda env create -f environment.yml
```

Activate conda environment.
```
conda activate chem_utils
```

## Scripts

- [`canonicalize_smiles.py`](https://github.com/swansonk14/chem_utils/blob/main/canonicalize_smiles.py): Canonicalizes SMILES using RDKit canonicalization and optionally strips salts.

- [`chemical_diversity.py`](https://github.com/swansonk14/chem_utils/blob/main/chemical_diversity.py): Computes the chemical diversity of a set of molecules in terms of Tanimoto distances.

- [`cluster_molecules.py`](https://github.com/swansonk14/chem_utils/blob/main/cluster_molecules.py): Performs k-means clustering to cluster molecules based on Morgan fingerprints.

- [`dimensionality_reduction.py`](https://github.com/swansonk14/chem_utils/blob/main/dimensionality_reduction.py): Visualizes molecules in 2D by performing dimensionality reduction (t-SNE or UMAP) on Morgan fingerprints.

- [`measure_experimental_reproducibility.py`](https://github.com/swansonk14/chem_utils/blob/main/measure_experimental_reproducibility.py): Measures the experimental reproducibility of two biological replicates by using one replicate to predict the other.

- [`molecular_weight_distribution.py`](https://github.com/swansonk14/chem_utils/blob/main/molecular_weight_distribution.py): Plots the distribution of molecular weights of a set of molecules.

- [`morgan_fingerprint.py`](https://github.com/swansonk14/chem_utils/blob/main/morgan_fingerprint.py): Contains functions that compute Morgan fingerprints. Parallelized for speed.

- [`nearest_neighbor_tanimoto.py`](https://github.com/swansonk14/chem_utils/blob/main/nearest_neighbor_tanimoto.py): Given a dataset of molecules, computes the nearest neighbor molecule based on Tanimoto similarity (with Morgan fingerprints) in a second dataset.

- [`regression_to_classification.py`](https://github.com/swansonk14/chem_utils/blob/main/regression_to_classification.py): Converts regression data to classification data using given thresholds.

- [`sdf_to_smiles.py`](https://github.com/swansonk14/chem_utils/blob/main/sdf_to_smiles.py): Converts an SDF file to a CSV file with SMILES.

- [`visualize_molecules.py`](https://github.com/swansonk14/chem_utils/blob/main/visualize_molecules.py): Converts a file of SMILES to images of molecular structures.

- [`visualize_reactions.py`](https://github.com/swansonk14/chem_utils/blob/main/visualize_reactions.py): Converts a file of reaction SMARTS to images of chemical reactions.
