"""Selects the best molecule from each cluster."""
from pathlib import Path

import pandas as pd
from tap import Tap

from constants import CLUSTER_COLUMN


class Args(Tap):
    data_path: Path  # Path to CSV file containing molecules.
    save_path: Path  # Path to CSV file where selected molecules will be saved.
    value_column: str  # Name of the column containing values to sort molecules by.
    cluster_column: str = CLUSTER_COLUMN  # Name of the column containing cluster labels.
    descending: bool = False  # Sorts molecules in each cluster from highest to lowest rather than lowest to highest.


def filter_molecules(data_path: Path,
                     save_path: Path,
                     value_column: str,
                     cluster_column: str = CLUSTER_COLUMN,
                     descending: bool = False) -> None:
    """Selects the best molecule from each cluster.

    :param data_path: Path to CSV file containing molecules.
    :param save_path: Path to CSV file where selected molecules will be saved.
    :param value_column: Name of the column containing values to sort molecules by.
    :param cluster_column: Name of the column containing cluster labels.
    :param descending: Sorts molecules in each cluster from highest to lowest rather than lowest to highest.
    """
    # Load data
    data = pd.read_csv(data_path)
    print(f'Original data size = {len(data):,}')

    # Select molecules
    data.sort_values(by=value_column, ascending=not descending, inplace=True)
    data.drop_duplicates(by=cluster_column, inplace=True)
    print(f'Final data size = {len(data):,}')

    # Save data
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)


if __name__ == '__main__':
    filter_molecules(**Args().parse_args().as_dict())
