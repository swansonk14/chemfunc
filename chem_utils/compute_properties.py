"""Computes one or more molecular properties a set of molecules."""
from multiprocessing import Pool
from pathlib import Path
from typing import Optional

import pandas as pd
from tqdm import tqdm

from chem_utils.constants import SMILES_COLUMN
from chem_utils.molecular_properties import get_available_property_functions, get_property_function


def compute_properties(data_path: Path,
                       properties: str,
                       save_path: Optional[Path] = None,
                       smiles_column: str = SMILES_COLUMN) -> None:
    """Computes one or more molecular properties a set of molecules.

    :param data_path: Path to a CSV file containing SMILES.
    :param properties: The name of the property to compute.
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


if __name__ == '__main__':
    from tap import Tap

    class Args(Tap):
        data_path: Path  # Path to a CSV file containing SMILES.
        properties: list[str]  # The names of the properties to compute.
        save_path: Optional[Path] = None  # Path to a CSV file where SMILES with properties is saved. If None, uses data_path.
        smiles_column: str = SMILES_COLUMN  # The name of the column in data_paths containing SMILES.

        def configure(self) -> None:
            self.add_argument('--properties', choices=get_available_property_functions())

    compute_properties(**Args().parse_args().as_dict())