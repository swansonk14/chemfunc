"""Measure experimental reproducibility by comparing biological replicates."""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import average_precision_score, mean_squared_error, r2_score, roc_auc_score


def measure_experimental_reproducibility(
        data_path: Path,
        rep_1_column: str,
        rep_2_column: str,
        threshold: float | None = None,
        plot: bool = False
) -> None:
    """Measure experimental reproducibility by comparing biological replicates.

    :param data_path: Path to CSV file containing two biological replicates.
    :param rep_1_column: Name of the column containing the first biological replicate.
    :param rep_2_column: Name of the column containing the second biological replicate.
    :param threshold: Threshold for binarizing the data.
    :param plot: Whether to plot the correlation between the two replicates.
    """
    # Load data
    data = pd.read_csv(data_path)
    print(f'Data size = {len(data):,}\n')

    # Get the two replicates
    rep_1 = data[rep_1_column]
    rep_2 = data[rep_2_column]

    # Plot the two replicates
    if plot:
        plt.scatter(rep_1, rep_2)
        plt.xlabel('Replicate 1')
        plt.ylabel('Replicate 2')
        plt.title('Replicate 1 vs 2')
        plt.show()

    # Measure experimental reproducibility for regression
    for score_name, score_function in [('R^2', r2_score), ('MSE', mean_squared_error)]:
        rep_1_to_rep_2 = score_function(rep_1, rep_2)
        rep_2_to_rep_1 = score_function(rep_2, rep_1)

        print(f'Replicate 1 to replicate 2 {score_name} = {rep_1_to_rep_2:.2f}')
        print(f'Replicate 2 to replicate 1 {score_name} = {rep_2_to_rep_1:.2f}\n')

    if threshold is not None:
        print(f'Using threshold = {threshold:.2f}\n')

        # Measure experimental reproducibility for classification
        rep_1_binary = rep_1 < threshold
        rep_2_binary = rep_2 < threshold

        print(f'Number of actives in replicate 1 = {sum(rep_1_binary):,}')
        print(f'Number of actives in replicate 2 = {sum(rep_2_binary):,}\n')

        for score_name, score_function in [('AUC', roc_auc_score), ('PRC-AUC', average_precision_score)]:
            # Note: negative needed since lower scores predict actives
            rep_1_predicts_rep_2 = score_function(rep_2_binary, -rep_1)
            rep_2_predicts_rep_1 = score_function(rep_1_binary, -rep_2)

            print(f'Replicate 1 predicts replicate 2 {score_name} = {rep_1_predicts_rep_2:.2f}')
            print(f'Replicate 2 predicts replicate 1 {score_name} = {rep_2_predicts_rep_1:.2f}\n')
