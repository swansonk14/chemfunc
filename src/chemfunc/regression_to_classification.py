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
        high_to_low: bool = False,
        delete_class_indices: set[int] | None = None
) -> None:
    """Converts regression data to classification data using provided thresholds.

    :param data_path: Path to CSV file containing regression data.
    :param regression_column: Name of the column containing regression data.
    :param classification_column: Name of the column where the classification data will be stored.
    :param thresholds: Thresholds to use to convert the regression data to classification data.
    :param save_path: Path to CSV file where classification data will be saved. If None, uses data_path.
    :param high_to_low: Whether class indices should be assigned from highest regression value to lowest.
    :param delete_class_indices: Optional set of class indices to delete from the data.
        Class indices will be reassigned to fill in the gaps.
        Use case: Create a binary classification dataset with a gap between classes by specifying two thresholds
        and deleting the class in between.
        Ex: thresholds = [0.4, 0.6], delete_class_indices = {1} will label data < 0.4 as 0 and data >= 0.6 as 1
        (originally labeled 2) and will delete data in between 0.4 and 0.6 (originally labeled 1).
    """
    # Load regression data
    data = pd.read_csv(data_path)
    print(f"Data size = {len(data):,}")

    # Convert regression data to classification
    data[classification_column] = np.digitize(x=data[regression_column], bins=sorted(thresholds))

    # Optionally, modify order of class indices
    if high_to_low:
        data[classification_column] = max(data[classification_column]) - data[classification_column]

    # Optionally, delete class indices
    if delete_class_indices is not None:
        # Delete data with specified class indices
        print(f"Deleting indices: {delete_class_indices}")
        data = data[~data[classification_column].isin(delete_class_indices)]
        print(f"Data size after deleting class indices = {len(data):,}")

        # Reassign class indices
        old_to_new_class_indices = {old_index: new_index for new_index, old_index in
                                    enumerate(sorted(data[classification_column].unique()))}
        data[classification_column] = data[classification_column].replace(old_to_new_class_indices)

    # Print counts of each class
    print(data[classification_column].value_counts())

    # Save classification data
    if save_path is None:
        save_path = data_path

    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)
