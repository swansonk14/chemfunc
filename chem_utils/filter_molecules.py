"""Filters molecules to those with values in a certain range."""
from pathlib import Path

import pandas as pd


def filter_molecules(data_path: Path,
                     save_path: Path,
                     filter_column: str,
                     min_value: float = -float('inf'),
                     max_value: float = float('inf')) -> None:
    """Filters molecules to those with values in a certain range.

    :param data_path: Path to CSV file containing molecules.
    :param save_path: Path to CSV file where filtered molecules will be saved.
    :param filter_column: Name of the column to use to filter molecules.
    :param min_value: Minimum value (inclusive) for the filter_column.
    :param max_value: Maximum value (inclusive) for the filter_column.
    """
    # Load data
    data = pd.read_csv(data_path)
    print(f'Original data size = {len(data):,}')

    # Filter molecules
    data = data[(data[filter_column] >= min_value) & (data[filter_column] <= max_value)]
    print(f'Filtered data size = {len(data):,}')

    # Save data
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)


if __name__ == '__main__':
    from tap import Tap

    class Args(Tap):
        data_path: Path  # Path to CSV file containing molecules.
        save_path: Path  # Path to CSV file where filtered molecules will be saved.
        filter_column: str  # Name of the column to use to filter molecules.
        min_value: float = -float('inf')  # Minimum value (inclusive) for the filter_column.
        max_value: float = float('inf')  # Maximum value (inclusive) for the filter_column.

    filter_molecules(**Args().parse_args().as_dict())
