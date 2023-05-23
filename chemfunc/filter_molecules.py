"""Filters molecules to those with values in a certain range."""
from pathlib import Path

import numpy as np
import pandas as pd


def filter_molecules(
        data_path: Path,
        save_path: Path,
        filter_column: str,
        min_value: float | None = None,
        max_value: float | None = None,
        bottom_proportion: float | None = None,
        top_proportion: float | None = None
) -> None:
    """Filters molecules to those with values in a certain range.

    :param data_path: Path to CSV file containing molecules.
    :param save_path: Path to CSV file where filtered molecules will be saved.
    :param filter_column: Name of the column to use to filter molecules.
    :param min_value: Minimum value (inclusive) for the filter_column.
    :param max_value: Maximum value (inclusive) for the filter_column.
    :param bottom_proportion: Keeps this proportion of the molecules with the lowest filter_column values.
                              E.g., bottom_proportion 0.2 keeps the molecules with the lowest 20% of values.
                              Note: Applied *after* filtering by min_value and max_value
                                    and at the same time as top_proportion.
    :param top_proportion: Keeps this proportion of the molecules with the highest filter_column values.
                           E.g., top_proportion 0.2 keeps the molecules with the highest 20% of values.
                           Note: Applied *after* filtering by min_value and max_value
                                 and at the same time as bottom_proportion.
    """
    # Load data
    data = pd.read_csv(data_path)
    print(f'Original data size = {len(data):,}')

    # Filter molecules by min_value
    if min_value is not None:
        data = data[data[filter_column] >= min_value]
        print(f'Data size after filtering with min_value {min_value} = {len(data):,}')

    # Filter molecules by max_value
    if max_value is not None:
        data = data[data[filter_column] <= max_value]
        print(f'Data size after filtering with max_value {max_value} = {len(data):,}')

    # Filter molecules by bottom_proportion or top_proportion
    if bottom_proportion is not None or top_proportion is not None:
        # Get argsort of data by filter_column values
        argsort = np.argsort(data[filter_column])

        if bottom_proportion is not None:
            if not (0.0 <= bottom_proportion <= 1.0):
                raise ValueError(f'Bottom proportion must be between 0.0 and 1.0 but got {bottom_proportion}.')

            num_bottom = round(len(data) * bottom_proportion)
            bottom_indices = set(argsort[:num_bottom])
            print(f'Keeping bottom {len(bottom_indices):,} molecules.')
        else:
            bottom_indices = set()

        if top_proportion is not None:
            if not (0.0 <= top_proportion <= 1.0):
                raise ValueError(f'Top proportion must be between 0.0 and 1.0 but got {top_proportion}.')

            num_top = round(len(data) * top_proportion)
            top_indices = set(argsort[-num_top:])
            print(f'Keeping top {len(top_indices):,} molecules.')
        else:
            top_indices = set()

        # Determine overlap if any
        if bottom_proportion is not None and top_proportion is not None:
            intersection_indices = bottom_indices & top_indices
            num_in_bottom_and_top = len(intersection_indices)

            if num_in_bottom_and_top > 0:
                print(f'Number of molecules in both bottom and top sets = {num_in_bottom_and_top:,}')

        # Combine bottom and top molecules and select molecules
        keep_indices = sorted(bottom_indices | top_indices)
        data = data.iloc[keep_indices]
        print(f'Data size after selecting bottom and/or top molecules = {len(data):,}')

    # Save data
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)
    print(f'Final data size = {len(data):,}')
