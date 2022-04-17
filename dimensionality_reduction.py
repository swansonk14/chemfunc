"""Runs dimensionality reduction (t-SNE or UMAP) on molecular fingerprints from one or more chemical libraries."""
import time
from pathlib import Path
from typing import Literal, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearnex import patch_sklearn
patch_sklearn()
from sklearn.manifold import TSNE
from tap import Tap
from tqdm import tqdm
from umap import UMAP

from morgan_fingerprint import compute_morgan_fingerprints


class Args(Tap):
    data_paths: list[Path]  # Path to CSV files containing SMILES.
    method: Literal['t-SNE', 'UMAP']  # Dimensionality reduction method.
    save_path: Path  # Path to a PDF file where the dimensionality reduction plot will be saved.
    max_molecules: Optional[list[int]] = None  # Maximum number of molecules sampled in each dataset.
    smiles_columns: Optional[list[str]] = None  # Name of the columns in the smiles_paths files containing SMILES.
    """If just one SMILES column is provided, it is applied to all files. Defaults to 'smiles'."""
    data_names: Optional[list[str]] = None  # Names of the data files for labeling the plot.
    highlight_data_names: Optional[set[str]] = None  # Names of the data files to highlight in the plot.


def dimesionality_reduction(data_paths: list[Path],
                            method: Literal['t-SNE', 'UMAP'],
                            save_path: Path,
                            max_molecules: Optional[list[int]] = None,
                            smiles_columns: Optional[list[Path]] = None,
                            data_names: Optional[list[str]] = None,
                            highlight_data_names: Optional[set[str]] = None):
    """Runs dimensionality reduction (t-SNE or UMAP) on molecular fingerprints from one or more chemical libraries.

    :param data_paths: Path to CSV files containing SMILES.
    :param method: Dimensionality reduction method.
    :param save_path: Path to a PDF file where the dimensionality reduction plot will be saved.
    :param max_molecules: Maximum number of molecules sampled in each dataset.
                          If just one is provided, it is applied to all files.
    :param smiles_columns: Name of the columns containing SMILES. If just one is provided, it is applied to all files.
                           Defaults to 'smiles'.
    :param data_names: Names of the data files for labeling the plot.
    :param highlight_data_names: Names of the data files to highlight in the plot.
    """
    # Validate max_molecules
    if max_molecules is None:
        max_molecules = [None] * len(data_paths)
    elif len(max_molecules) == 1:
        max_molecules = max_molecules * len(data_paths)
    elif len(max_molecules) != len(data_paths):
        raise ValueError('Number of max_molecules does not match number of data paths.')

    # Validate smiles_columns
    if smiles_columns is None:
        smiles_columns = ['smiles'] * len(data_paths)
    elif len(smiles_columns) == 1:
        smiles_columns = smiles_columns * len(data_paths)
    elif len(smiles_columns) != len(data_paths):
        raise ValueError('Number of SMILES columns does not match number of data paths.')

    # Validate data_names
    if data_names is None:
        data_names = [data_path.stem for data_path in data_paths]
    elif len(data_names) != len(data_paths):
        raise ValueError('Number of data names does not match number of data paths.')

    # Validate highlight_data_names
    if highlight_data_names is None:
        highlight_data_names = set()

    # Set colors
    if len(data_paths) <= 10:
        cmap = plt.get_cmap('tab10')
    elif len(data_paths) <= 20:
        cmap = plt.get_cmap('tab20')
    else:
        raise ValueError('Not enough colors for more than 20 data paths.')

    # Load data and subsample SMILES
    smiles, slices = [], []
    for data_path, smiles_column, data_name, max_mols in tqdm(zip(data_paths, smiles_columns, data_names, max_molecules),
                                                              total=len(data_paths), desc='Loading data'):
        # Load data
        new_smiles = list(pd.read_csv(data_path)[smiles_column])
        print(f'{data_name}: {len(new_smiles):,}')

        # Subsample if dataset is too large
        if max_mols is not None and len(new_smiles) > max_mols:
            print(f'Subsampling to {max_mols:,} molecules')
            np.random.seed(0)
            new_smiles = np.random.choice(new_smiles, size=max_mols, replace=False).tolist()

        slices.append(slice(len(smiles), len(smiles) + len(new_smiles)))
        smiles += new_smiles

    # Compute Morgan fingerprints
    morgans = compute_morgan_fingerprints(smiles)

    # Run dimensionality reduction
    if method == 't-SNE':
        reducer = TSNE(random_state=0, metric='jaccard', init='pca', n_jobs=-1, square_distances=True)
    elif method == 'UMAP':
        reducer = UMAP(random_state=0, metric='jaccard')
    else:
        raise ValueError(f'Dimensionality reduction method "{method}" is not supported.')

    print(f'Running {method}')
    start = time.time()
    X = reducer.fit_transform(morgans)
    print(f'time = {time.time() - start:.2f} seconds')

    print('Plotting')
    x_min, x_max = np.min(X, axis=0), np.max(X, axis=0)
    X = (X - x_min) / (x_max - x_min)

    plt.clf()
    plt.figure(figsize=(64, 48))
    plt.title(f'{method} using Morgan fingerprint with Jaccard similarity', fontsize=100)

    for index, (slc, data_name) in enumerate(zip(slices, data_names)):
        plt.scatter(
            X[slc, 0],
            X[slc, 1],
            s=600 if data_name in highlight_data_names else 200,
            color=cmap(index),
            label=data_name,
            marker='*' if data_name in highlight_data_names else '.'
        )

    plt.legend(fontsize=50)
    plt.xticks([]), plt.yticks([])

    print('Saving plot')
    save_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(save_path)


if __name__ == '__main__':
    dimesionality_reduction(**Args().parse_args().as_dict())