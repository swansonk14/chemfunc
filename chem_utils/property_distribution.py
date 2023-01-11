"""Computes the distribution of property values over a set of molecules."""
from multiprocessing import Pool
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
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
                          save_dir: Path,
                          prop: Literal['mol_weight', 'logp'],
                          smiles_column: str = SMILES_COLUMN,
                          min_value: float = -float('inf'),
                          max_value: float = float('inf'),
                          save_data: bool = False) -> None:
    """Computes the distribution of property values over a set of molecules.

    :param data_paths: Path to CSV files containing SMILES.
    :param save_dir: Path to a directory where the plot and data will be saved.
    :param prop: The name of the property to compute.
    :param smiles_column: The name of the column in data_path containing SMILES.
    :param min_value: Minimum property value to plot (removes outliers).
    :param max_value: Maximum property value to plot (removes outliers).
    :param save_data: Whether to save the property data as a CSV alongside the plot.
    """
    # Select property function
    if prop == 'mol_weight':
        prop_fn = molecular_weight
    elif prop == 'logp':
        prop_fn = logP
    else:
        raise ValueError(f'Molecular property "{prop}" is not supported.')

    # Iterate over data paths
    prop_data = {}
    for data_path in tqdm(data_paths, desc='Data paths'):
        # Load data
        data = pd.read_csv(data_path)

        # Compute property values
        with Pool() as pool:
            data[prop] = list(tqdm(pool.imap(prop_fn, data[smiles_column]), total=len(data), desc=prop))

        # Plot property values
        plt.hist(data[prop][(data[prop] >= min_value) & (data[prop] <= max_value)],
                 bins=100, density=True, label=data_path.stem, alpha=1 / len(data_paths))

        # Save data with property values
        prop_data[data_path.stem] = data[prop]

    # Label plot
    plt.xlabel(prop)
    plt.ylabel('Density')
    plt.title(f'{prop} Distribution')
    plt.legend()

    # Save plot
    save_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(save_dir / f'{prop}.pdf', bbox_inches='tight')

    # Save data
    if save_data:
        max_len = max(len(values) for values in prop_data.values())
        fig_data = pd.DataFrame({
            key: np.pad(values, (0, max_len - len(values)), constant_values=np.nan)
            for key, values in prop_data.items()
        })
        fig_data.to_csv(save_dir / f'{prop}.csv', index=False)


if __name__ == '__main__':
    from tap import Tap

    class Args(Tap):
        data_paths: list[Path]  # Path to CSV files containing SMILES.
        save_dir: Path  # Path to a directory where the plot and (optionally) data will be saved.
        prop: Literal['mol_weight', 'logp']  # The name of the property to compute.
        smiles_column: str = SMILES_COLUMN  # The name of the column in data_paths containing SMILES.
        min_value: float = -float('inf')  # Minimum property value to plot (removes outliers).
        max_value: float = float('inf')  # Maximum property value to plot (removes outliers).
        save_data: bool = False  # Whether to save the property data as a CSV alongside the plot.

    property_distribution(**Args().parse_args().as_dict())
