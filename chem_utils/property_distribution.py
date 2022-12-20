"""Computes the distribution of molecular weights over a set of molecules."""
from multiprocessing import Pool
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.Descriptors import MolWt
from tqdm import tqdm

from chem_utils.constants import SMILES_COLUMN


def logP(smiles: str) -> float:
    """Computes the logP of a molecule.

    :param smiles: A SMILES representing a molecule.
    :return: The logP of the molecule.
    """
    mol = Chem.MolFromSmiles(smiles)
    logp = MolLogP(mol)

    return logp


def molecular_weight(smiles: str) -> float:
    """Computes the molecular weight of a molecule.

    :param smiles: A SMILES representing a molecule.
    :return: The molecular weight of the molecule.
    """
    mol = Chem.MolFromSmiles(smiles)
    mol_weight = MolWt(mol)

    return mol_weight


def property_distribution(data_paths: list[Path],
                          save_path: Path,
                          property: Literal['mol_weight', 'logp'],
                          smiles_column: str = SMILES_COLUMN,
                          min_value: float = -float('inf'),
                          max_value: float = float('inf')) -> None:
    """Computes the distribution of molecular weights over a set of molecules.

    :param data_paths: Path to CSV files containing SMILES.
    :param smiles_column: The name of the column in data_path containing SMILES.
    :param property: The name of the property to compute.
    :param save_path: Path to a PDF file where the plot of molecular weights will be saved.
    :param min_value: Minimum molecular weight to plot (removes outliers).
    :param max_value: Maximum molecular weight to plot (removes outliers).
    """
    # Select property function
    if property == 'mol_weight':
        property_fn = molecular_weight
    elif property == 'logp':
        property_fn = logP
    else:
        raise ValueError(f'Molecular property "{property}" is not supported.')

    # Iterate over data paths
    for data_path in tqdm(data_paths, desc='Data paths'):
        # Load data
        data = pd.read_csv(data_path)

        # Compute molecular weights
        with Pool() as pool:
            values = [
                value
                for value in tqdm(pool.imap_unordered(property_fn, data[smiles_column]),
                                       total=len(data), desc=property)
                if min_value <= value <= max_value
            ]

        # Plot molecular weights
        plt.hist(values, bins=100, density=True, label=data_path.stem, alpha=1 / len(data_paths))

    # Label plot
    plt.xlabel(property)
    plt.ylabel('Density')
    plt.title(f'{property} Distribution')
    plt.legend()

    # Save plot
    save_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(save_path, bbox_inches='tight')


if __name__ == '__main__':
    from tap import Tap

    class Args(Tap):
        data_paths: list[Path]  # Path to CSV files containing SMILES.
        save_path: Path  # Path to a PDF file where the plot of molecular weights will be saved.
        property: Literal['mol_weight', 'logp']  # The name of the property to compute.
        smiles_column: str = SMILES_COLUMN  # The name of the column in data_paths containing SMILES.
        min_value: float = -float('inf')  # Minimum property value to plot (removes outliers).
        max_value: float = float('inf')  # Maximum property value to plot (removes outliers).

    property_distribution(**Args().parse_args().as_dict())
