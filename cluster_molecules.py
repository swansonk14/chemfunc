"""Clusters molecules by Morgan fingerprint."""
from pathlib import Path

import pandas as pd
from sklearnex import patch_sklearn
patch_sklearn()
from sklearn.cluster import KMeans
from tap import Tap

from constants import SMILES_COLUMN
from morgan_fingerprint import compute_morgan_fingerprints


class Args(Tap):
    data_path: Path  # Path to CSV file containing SMILES.
    save_path: Path  # Path to CSV file where the clustering will be saved.
    smiles_column: str = SMILES_COLUMN  # Name of the column containing SMILES.
    num_clusters: int = 50  # Number of clusters.


def cluster_molecules(data_path: Path,
                      save_path: Path,
                      smiles_column: str = SMILES_COLUMN,
                      num_clusters: int = 50) -> None:
    """Clusters molecules by Morgan fingerprint.

    :param data_path: Path to CSV file containing SMILES.
    :param save_path: Path to CSV file where the clustering will be saved.
    :param smiles_column: Name of the column containing SMILES.
    :param num_clusters: Number of clusters.
    """
    print('Loading data')
    data = pd.read_csv(data_path)
    print(f'Data size = {len(data):,}')

    print('Computing Morgan fingerprints')
    morgans = compute_morgan_fingerprints(data[smiles_column])

    print('Clustering')
    kmeans = KMeans(n_clusters=num_clusters, random_state=0).fit(morgans)
    data['cluster_label'] = kmeans.labels_ + 1

    print('Saving data')
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)


if __name__ == '__main__':
    cluster_molecules(**Args().parse_args().as_dict())
