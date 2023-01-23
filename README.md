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

Pip install package.
```
pip install -e .
```

## Scripts

- [`canonicalize_smiles.py`](https://github.com/swansonk14/chem_utils/blob/main/chem_utils/canonicalize_smiles.py): Canonicalizes SMILES using RDKit canonicalization and optionally strips salts.

- [`chemical_diversity.py`](https://github.com/swansonk14/chem_utils/blob/main/chem_utils/chemical_diversity.py): Computes the chemical diversity of a set of molecules in terms of Tanimoto distances.

- [`cluster_molecules.py`](https://github.com/swansonk14/chem_utils/blob/main/chem_utils/cluster_molecules.py): Performs k-means clustering to cluster molecules based on Morgan fingerprints.

- [`compute_property_distribution.py`](https://github.com/swansonk14/chem_utils/blob/main/chem_utils/plot_property_distribution.py): Computes one or more molecular properties for a set of molecules.

- [`dimensionality_reduction.py`](https://github.com/swansonk14/chem_utils/blob/main/chem_utils/dimensionality_reduction.py): Visualizes molecules in 2D by performing dimensionality reduction (t-SNE or UMAP) on Morgan fingerprints.

- [`filter_molecules.py`](https://github.com/swansonk14/chem_utils/blob/main/chem_utils/filter_molecules.py): Filters molecules to those with values in a certain range.

- [`measure_experimental_reproducibility.py`](https://github.com/swansonk14/chem_utils/blob/main/chem_utils/measure_experimental_reproducibility.py): Measures the experimental reproducibility of two biological replicates by using one replicate to predict the other.

- [`molecular_fingerprints.py`](https://github.com/swansonk14/chem_utils/blob/main/chem_utils/molecular_fingerprints.py): Contains functions to compute fingerprints for molecules. Parallelized for speed.

- [`molecular_properties.py`](https://github.com/swansonk14/chem_utils/blob/main/chem_utils/molecular_properties.py): Contains functions to compute molecular properties. Parallelized for speed.

- [`molecular_similarities.py`](https://github.com/swansonk14/chem_utils/blob/main/chem_utils/molecular_similarities.py): Contains functions to compute similarities between molecules. Parallelized for speed.

- [`nearest_neighbor.py`](https://github.com/swansonk14/chem_utils/blob/main/chem_utils/nearest_neighbor.py): Given a dataset of molecules, computes the nearest neighbor molecule in a second dataset using one of several similarity metrics.

- [`plot_property_distribution.py`](https://github.com/swansonk14/chem_utils/blob/main/chem_utils/plot_property_distribution.py): Plots the distribution of molecular properties of a set of molecules.

- [`regression_to_classification.py`](https://github.com/swansonk14/chem_utils/blob/main/chem_utils/regression_to_classification.py): Converts regression data to classification data using given thresholds.

- [`sample_molecules.py`](https://github.com/swansonk14/chem_utils/blob/main/chem_utils/sample_molecules.py): Samples molecules from a CSV file, either uniformly at random across the entire dataset or uniformly at random from each cluster within the data.

- [`sdf_to_smiles.py`](https://github.com/swansonk14/chem_utils/blob/main/chem_utils/sdf_to_smiles.py): Converts an SDF file to a CSV file with SMILES.

- [`smiles_to_svg.py`](https://github.com/swansonk14/chem_utils/blob/main/chem_utils/smiles_to-svg.py): Converts a SMILES string to an SVG image of the molecule.

- [`select_from_clusters.py`](https://github.com/swansonk14/chem_utils/blob/main/chem_utils/select_from_clusters.py): Selects the best molecule from each cluster.

- [`visualize_molecules.py`](https://github.com/swansonk14/chem_utils/blob/main/chem_utils/visualize_molecules.py): Converts a file of SMILES to images of molecular structures.

- [`visualize_reactions.py`](https://github.com/swansonk14/chem_utils/blob/main/chem_utils/visualize_reactions.py): Converts a file of reaction SMARTS to images of chemical reactions.
