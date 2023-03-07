"""Plots the distribution of property values over one or more sets of molecules."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm


def plot_property_distribution(data_paths: list[Path],
                               property_column: str,
                               save_dir: Path,
                               min_value: float = -float('inf'),
                               max_value: float = float('inf')) -> None:
    """Plots the distribution of property values over one or more sets of molecules

    :param data_paths: Path to CSV files containing SMILES.
    :param property_column: Name of the column containing the property values.
    :param save_dir: Path to a directory where the plot and data will be saved.
    :param min_value: Minimum property value to plot (removes outliers).
    :param max_value: Maximum property value to plot (removes outliers).
    """
    # Iterate over data paths
    prop_data = {}
    for data_path in tqdm(data_paths, desc='Data paths'):
        # Load data
        data = pd.read_csv(data_path)

        # Plot property values
        plt.hist(
            data[property_column][(data[property_column] >= min_value) & (data[property_column] <= max_value)],
            bins=100, density=True, label=data_path.stem, alpha=1 / len(data_paths)
        )

        # Save data with property values
        prop_data[data_path.stem] = data[property_column]

    # Label plot
    plt.xlabel(property_column)
    plt.ylabel('Density')
    plt.title(f'{property_column} Distribution')
    plt.legend()

    # Save plot
    save_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(save_dir / f'{property_column}.pdf', bbox_inches='tight')

    # Save data
    max_len = max(len(values) for values in prop_data.values())
    fig_data = pd.DataFrame({
        key: np.pad(values, (0, max_len - len(values)), constant_values=np.nan)
        for key, values in prop_data.items()
    })
    fig_data.to_csv(save_dir / f'{property_column}.csv', index=False)


if __name__ == '__main__':
    from tap import tapify

    tapify(plot_property_distribution)
