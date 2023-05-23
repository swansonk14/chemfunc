"""Selects the best molecule from each cluster."""
from pathlib import Path

import pandas as pd

from chemfunc.constants import CLUSTER_COLUMN


def select_from_clusters(
        data_path: Path,
        save_path: Path,
        value_column: str,
        cluster_column: str = CLUSTER_COLUMN,
        ascending: bool = False
) -> None:
    """Selects the best molecule from each cluster.

    :param data_path: Path to CSV file containing molecules.
    :param save_path: Path to CSV file where selected molecules will be saved.
    :param value_column: Name of the column containing values to sort molecules by.
    :param cluster_column: Name of the column containing cluster labels.
    :param ascending: Sorts molecules in each cluster from lowest to highest rather than highest to lowest.
    """
    # Load data
    data = pd.read_csv(data_path)
    print(f'Original data size = {len(data):,}')

    # Select molecules
    data.sort_values(by=value_column, ascending=ascending, inplace=True)
    data.drop_duplicates(subset=[cluster_column], inplace=True)
    print(f'Final data size = {len(data):,}')

    # Save data
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)
