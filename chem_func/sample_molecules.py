"""Samples molecules, either uniformly at random across the entire dataset or uniformly at random from each cluster."""
from pathlib import Path

import pandas as pd
from tqdm import tqdm


def sample_molecules(
        data_path: Path,
        save_path: Path,
        num_molecules: int,
        cluster_column: str | None = None
) -> None:
    """Samples molecules, either uniformly at random across the entire dataset or uniformly at random from each cluster.

    :param data_path: Path to CSV file containing SMILES.
    :param save_path: Path to CSV file where the selected molecules will be saved.
    :param num_molecules: Number of molecules to select.
    :param cluster_column: Name of the column containing cluster labels.
                           If None, molecules are selected uniformly at random.
                           If provided, molecules are selected uniformly at random in
                           each cluster with num_molecules per cluster.
    """
    print('Loading data')
    data = pd.read_csv(data_path)
    print(f'Data size = {len(data):,}')

    if cluster_column is not None:
        print(f'Selecting {num_molecules:,} molecules per cluster')
        sampled = pd.concat([
            data[data[cluster_column] == cluster_label].sample(n=num_molecules, random_state=0)
            for cluster_label in tqdm(sorted(data[cluster_column].unique()))
        ])
    else:
        print(f'Selection {num_molecules:,} molecules')
        sampled = data.sample(n=num_molecules, random_state=0)

    print('Saving data')
    save_path.parent.mkdir(parents=True, exist_ok=True)
    sampled.to_csv(save_path, index=False)
