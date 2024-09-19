"""Computes one or more molecular properties a set of molecules."""
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from chemfunc.constants import SMILES_COLUMN
from chemfunc.molecular_properties import get_property_function


def compute_properties(data_path: Path,
                       properties: list[str],
                       save_path: Path | None = None,
                       smiles_column: str = SMILES_COLUMN) -> None:
    """Computes one or more molecular properties a set of molecules.

    :param data_path: Path to a CSV file containing SMILES.
    :param properties: The name of the properties to compute.
    :param save_path: Path to a CSV file where SMILES with properties is saved. If None, uses data_path.
    :param smiles_column: The name of the column in data_path containing SMILES.
    """
    # Load data
    data = pd.read_csv(data_path)

    # Iterate over properties
    for property_name in properties:
        # Get property function
        property_fn = get_property_function(property_name)

        # Compute property values
        with Pool() as pool:
            data[property_name] = list(
                tqdm(pool.imap(property_fn, data[smiles_column]),
                     total=len(data), desc=property_name)
            )

    # Set up save path
    if save_path is None:
        save_path = data_path

    save_path.parent.mkdir(parents=True, exist_ok=True)

    # Save data
    data.to_csv(save_path, index=False)
