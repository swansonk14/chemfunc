"""Deduplicate a CSV file by SMILES."""
from pathlib import Path

import pandas as pd

from chemfunc.constants import SMILES_COLUMN


def deduplicate_smiles(
        data_path: Path,
        save_path: Path,
        smiles_column: str = SMILES_COLUMN
) -> None:
    """Deduplicate a CSV file by SMILES.

    :param data_path: Path to CSV file containing SMILES.
    :param save_path: Path to CSV file where deduplicated SMILES will be saved.
    :param smiles_column: Name of the column containing SMILES to deduplicate.
    """
    # Load data
    data = pd.read_csv(data_path)

    print(f'Data size = {len(data):,}')

    # Deduplicate by SMILES
    data.drop_duplicates(smiles_column, inplace=True)

    print(f'Data size after deduplication = {len(data):,}')

    # Save data
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)
