"""Converts regression data to classification data using provided thresholds."""
from pathlib import Path

import numpy as np
import pandas as pd


def regression_to_classification(
        data_path: Path,
        regression_column: str,
        classification_column: str,
        thresholds: list[float],
        save_path: Path | None = None,
        high_to_low: bool = False
) -> None:
    """Converts regression data to classification data using provided thresholds.

    :param data_path: Path to CSV file containing regression data.
    :param regression_column: Name of the column containing regression data.
    :param classification_column: Name of the column where the classification data will be stored.
    :param thresholds: Thresholds to use to convert the regression data to classification data.
    :param save_path: Path to CSV file where classification data will be saved. If None, uses data_path.
    :param high_to_low: Whether class indices should be assigned from highest regression value to lowest.
    """
    # Load regression data
    data = pd.read_csv(data_path)

    # Convert regression data to classification
    data[classification_column] = np.digitize(x=data[regression_column], bins=sorted(thresholds))

    # Modify order of class indices
    if high_to_low:
        data[classification_column] = max(data[classification_column]) - data[classification_column]

    # Print counts of each class
    print(data[classification_column].value_counts())

    # Save classification data
    if save_path is None:
        save_path = data_path

    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)
