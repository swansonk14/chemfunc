"""Computes the chemical diversity of a list of molecules in terms of Tanimoto distances."""
from pathlib import Path

import numpy as np
import pandas as pd
from tap import Tap

from constants import SMILES_COLUMN
from nearest_neighbor_tanimoto import compute_pairwise_tanimoto_distances


class Args(Tap):
    data_path: Path  # Path to CSV file containing molecule SMILES.
    smiles_column: str = SMILES_COLUMN  # Name of the column containing SMILES.


def chemical_diversity(data_path: Path,
                       smiles_column: str = SMILES_COLUMN) -> None:
    """Computes the chemical diversity of a list of molecules in terms of Tanimoto distances.

    Note: Does NOT remove duplicate SMILES before computing pairwise distances.

    :param data_path: Path to CSV file containing molecule SMILES.
    :param smiles_column: Name of the column containing SMILES.
    """
    # Load data
    data = pd.read_csv(data_path)
    print(f'Number of molecules = {len(data):,}')

    # Compute pairwise Tanimoto distances
    pairwise_tanimoto_distances = compute_pairwise_tanimoto_distances(mols_1=data[smiles_column])

    # Compute average Tanimoto distance between pairs
    unique_distances_mask = np.triu(np.ones(pairwise_tanimoto_distances.shape), k=1).astype(bool)
    unique_pairwise_tanimoto_distances = pairwise_tanimoto_distances[unique_distances_mask]
    average_tanimoto_distance = np.mean(unique_pairwise_tanimoto_distances)
    std_tanimoto_distance = np.std(unique_pairwise_tanimoto_distances)

    print(f'Average Tanimoto distance between pairs = '
          f'{average_tanimoto_distance:.3f} +/- {std_tanimoto_distance:.3f}')

    # Compute average minimum Tanimoto distance between each molecule and the other molecules
    np.fill_diagonal(pairwise_tanimoto_distances, np.inf)
    min_tanimoto_distances = np.min(pairwise_tanimoto_distances, axis=1)
    average_min_tanimoto_distance = np.mean(min_tanimoto_distances)
    std_min_tanimoto_distance = np.std(min_tanimoto_distances)

    print(f'Average minimum Tanimoto distance from each molecule to the rest = '
          f'{average_min_tanimoto_distance:.3f} +/- {std_min_tanimoto_distance:.3f}')


if __name__ == '__main__':
    chemical_diversity(**Args().parse_args().as_dict())
