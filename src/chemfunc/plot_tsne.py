"""Runs a t-SNE on molecular fingerprints from one or more chemical libraries."""
import time
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
from tqdm import tqdm

from chemfunc.constants import SMILES_COLUMN
from chemfunc.molecular_fingerprints import compute_fingerprints


def plot_tsne(
        data_paths: list[Path],
        save_path: Path,
        metric: Literal['jaccard', 'euclidean'] = 'jaccard',
        embedder: Literal['morgan', 'file'] = 'morgan',
        max_molecules: list[int] | None = None,
        colors: list[str] | None = None,
        smiles_columns: list[str] | None = None,
        data_names: list[str] | None = None,
        highlight_data_names: set[str] | None = None,
        display_data_names: set[str] | None = None,
        point_size: int = 1000
) -> None:
    """Runs a t-SNE on molecular fingerprints from one or more chemical libraries.

    :param data_paths: Path to CSV files containing SMILES.
    :param save_path: Path to a PDF file where the t-SNE plot will be saved.
    :param metric: Metric to use to compared embeddings.
    :param embedder: Embedding to use for the molecules.
                     morgan: Computes Morgan fingerprint from the SMILES.
                     file: Uses all columns except the SMILES column from the data file as the embedding.
    :param max_molecules: Maximum number of molecules sampled in each dataset.
                          If just one is provided, it is applied to all files.
    :param colors: The colors to use for each dataset. If None, uses tab10 or tab20 depending on the number of datasets.
    :param smiles_columns: Name of the columns containing SMILES. If just one is provided, it is applied to all files.
                           Defaults to 'smiles'.
    :param data_names: Names of the data files for labeling the plot.
    :param highlight_data_names: Names of the data files to highlight in the plot.
    :param display_data_names: The names of the data files to display in the plot. If None, all are displayed.
    :param point_size: The size of the points in the plot.
    """
    # Validate max_molecules
    if max_molecules is None:
        max_molecules = [None] * len(data_paths)
    elif len(max_molecules) == 1:
        max_molecules = max_molecules * len(data_paths)
    elif len(max_molecules) != len(data_paths):
        raise ValueError('Number of max_molecules does not match number of data paths.')

    # Validate colors
    if colors is not None and len(colors) != len(data_paths):
        raise ValueError('Number of colors does not match number of data paths.')

    # Validate smiles_columns
    if smiles_columns is None:
        smiles_columns = [SMILES_COLUMN] * len(data_paths)
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
    if colors is None:
        if len(data_paths) <= 10:
            cmap = plt.get_cmap('tab10')
        elif len(data_paths) <= 20:
            cmap = plt.get_cmap('tab20')
        else:
            raise ValueError('Not enough colors for more than 20 data paths.')

        colors = [cmap(i) for i in range(len(data_paths))]

    # Load data and subsample SMILES
    smiles, slices, embeddings = [], [], []
    for data_path, smiles_column, data_name, max_mols in tqdm(
            zip(data_paths, smiles_columns, data_names, max_molecules),
            total=len(data_paths), desc='Loading data'):
        # Load data
        data = pd.read_csv(data_path)
        print(f'{data_name}: {len(data):,}')

        # Subsample if dataset is too large
        if max_mols is not None and len(data) > max_mols:
            print(f'Subsampling to {max_mols:,} molecules')
            data = data.sample(n=max_mols, replace=False, random_state=0)

        # Update slices
        slices.append(slice(len(smiles), len(smiles) + len(data)))

        # Update SMILES
        smiles += list(data[smiles_column])

        # Get molecule embeddings if provided
        if embedder == 'file':
            embedding_columns = list(data.columns)
            embedding_columns.remove(smiles_column)
            embeddings.append(data[embedding_columns].to_numpy())

    # Get/compute molecule embeddings
    if embedder == 'morgan':
        embeddings = compute_fingerprints(smiles, fingerprint_type='morgan')
    elif embedder == 'file':
        embeddings = np.concatenate(embeddings)
    else:
        raise ValueError(f'Embedder "{embedder}" is not supported.')

    # Run t-SNE
    tsne = TSNE(random_state=0, metric=metric, init='pca', n_jobs=-1)

    print(f'Running t-SNE')
    start = time.time()
    X = tsne.fit_transform(embeddings)
    print(f'time = {time.time() - start:.2f} seconds')

    print('Plotting')
    x_min, x_max = np.min(X, axis=0), np.max(X, axis=0)
    X = (X - x_min) / (x_max - x_min)

    plt.clf()
    plt.figure(figsize=(64, 48))
    plt.title(f't-SNE using Morgan fingerprint with {metric.title()} similarity', fontsize=100)

    tsne_data = {}
    for index, (slc, data_name) in enumerate(zip(slices, data_names)):
        if display_data_names is None or data_name in display_data_names:
            plt.scatter(
                X[slc, 0],
                X[slc, 1],
                s=1.5 * point_size if data_name in highlight_data_names else point_size,
                color=colors[index],
                label=data_name,
                marker='*' if data_name in highlight_data_names else '.'
            )

            tsne_data[f'{data_name}_x'] = X[slc, 0]
            tsne_data[f'{data_name}_y'] = X[slc, 1]

    plt.legend(fontsize=50)
    plt.xticks([]), plt.yticks([])
    plt.tight_layout()

    print('Saving plot')
    save_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(save_path, bbox_inches='tight')
