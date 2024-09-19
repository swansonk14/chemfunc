"""Clusters molecules by Morgan fingerprint."""
from pathlib import Path

import pandas as pd
from sklearn.cluster import KMeans, MiniBatchKMeans

from chemfunc.constants import CLUSTER_COLUMN, SMILES_COLUMN
from chemfunc.molecular_similarities import compute_fingerprints


def cluster_molecules(
        data_path: Path,
        save_path: Path | None = None,
        smiles_column: str = SMILES_COLUMN,
        num_clusters: int = 50,
        mini_batch: bool = False
) -> None:
    """Clusters molecules by Morgan fingerprint.

    :param data_path: Path to CSV file containing SMILES.
    :param save_path: Path to CSV file where the clustering will be saved (defaults to data_path).
    :param smiles_column: Name of the column containing SMILES.
    :param num_clusters: Number of clusters.
    :param mini_batch: Whether to use mini batch k-means instead of standard k-means.
    """
    # Set save path
    if save_path is None:
        save_path = data_path

    print('Loading data')
    data = pd.read_csv(data_path)
    print(f'Data size = {len(data):,}')

    print('Computing Morgan fingerprints')
    morgans = compute_fingerprints(data[smiles_column], fingerprint_type='morgan')

    print('Clustering')
    if mini_batch:
        kmeans = MiniBatchKMeans(n_clusters=num_clusters, random_state=0)
    else:
        kmeans = KMeans(n_clusters=num_clusters, random_state=0)

    kmeans.fit(morgans)

    data[CLUSTER_COLUMN] = kmeans.labels_ + 1

    print('Saving data')
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)
