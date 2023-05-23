"""Computes the chemical diversity of a list of molecules in terms of chemical distances."""
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd

from chemfunc.constants import SMILES_COLUMN
from chemfunc.molecular_similarities import get_similarity_function


def chemical_diversity(
        data_path: Path,
        smiles_column: str = SMILES_COLUMN,
        similarity_type: Literal['tanimoto', 'tversky'] = 'tanimoto'
) -> None:
    """Computes the chemical diversity of a list of molecules in terms of chemical distances.

    Note: Does NOT remove duplicate SMILES before computing pairwise distances.

    :param data_path: Path to CSV file containing molecule SMILES.
    :param smiles_column: Name of the column containing SMILES.
    :param similarity_type: Type of similarity to use when computing distances.
    """
    # Load data
    data = pd.read_csv(data_path)
    print(f'Number of molecules = {len(data):,}')

    # Compute pairwise distances
    pairwise_distances = 1 - get_similarity_function(similarity_type)(mols_1=data[smiles_column])

    # Compute average distance between pairs
    unique_distances_mask = np.triu(np.ones(pairwise_distances.shape), k=1).astype(bool)
    unique_pairwise_distances = pairwise_distances[unique_distances_mask]
    average_distance = np.mean(unique_pairwise_distances)
    std_distance = np.std(unique_pairwise_distances)

    print(f'Average {similarity_type.title()} distance between pairs = '
          f'{average_distance:.3f} +/- {std_distance:.3f}')

    # Compute average minimum distance between each molecule and the other molecules
    np.fill_diagonal(pairwise_distances, np.inf)
    min_distances = np.min(pairwise_distances, axis=1)
    average_min_distance = np.mean(min_distances)
    std_min_distance = np.std(min_distances)

    print(f'Average minimum {similarity_type.title()} distance from each molecule to the rest = '
          f'{average_min_distance:.3f} +/- {std_min_distance:.3f}')
