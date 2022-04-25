"""Converts regression data to classification data using provided thresholds."""
from pathlib import Path

import numpy as np
import pandas as pd
from tap import Tap


class Args(Tap):
    data_path: Path  # Path to CSV file containing regression data.
    regression_column: str  # Name of the column containing regression data.
    classification_column: str  # Name of the column where the classification data will be stored.
    thresholds: list[float]  # Thresholds to use to convert the regression data to classification data.
    save_path: Path  # Path to CSV file where classification data will be saved.
    high_to_low: bool = False  # Whether class indices should be assigned from highest regression value to lowest.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def regression_to_classification(data_path: Path,
                                 regression_column: str,
                                 classification_column: str,
                                 thresholds: list[str],
                                 save_path: Path,
                                 high_to_low: bool = False) -> None:
    """Converts regression data to classification data using provided thresholds.

    :param data_path: Path to CSV file containing regression data.
    :param regression_column: Name of the column containing regression data.
    :param classification_column: Name of the column where the classification data will be stored.
    :param thresholds: Thresholds to use to convert the regression data to classification data.
    :param save_path: Path to CSV file where classification data will be saved.
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
    data.to_csv(save_path, index=False)


if __name__ == '__main__':
    regression_to_classification(**Args().parse_args().as_dict())
