"""Computes the distribution of molecular weights over a set of molecules."""
from multiprocessing import Pool
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from tqdm import tqdm

from chem_utils.constants import SMILES_COLUMN


def molecular_weight(smiles: str) -> float:
    """Computes the molecular weight of a molecule.

    :param smiles: A SMILES representing a molecule.
    :return: The molecular weight of the molecule.
    """
    mol = Chem.MolFromSmiles(smiles)
    mol_weight = MolWt(mol)

    return mol_weight


def molecular_weight_distribution(data_path: Path,
                                  save_path: Path,
                                  smiles_column: str = SMILES_COLUMN) -> None:
    """Computes the distribution of molecular weights over a set of molecules.

    :param data_path: Path to a CSV file containing SMILES.
    :param smiles_column: The name of the column in data_path containing SMILES.
    :param save_path: Path to a PDF file where the plot of molecular weights will be saved.
    """
    # Load data
    data = pd.read_csv(data_path)

    # Compute molecular weights
    with Pool() as pool:
        mol_weights = list(tqdm(pool.imap_unordered(molecular_weight, data[smiles_column]),
                                total=len(data), desc='Molecular weight'))

    # Plot molecular weights
    plt.hist(mol_weights, bins=100, density=True)
    plt.xlabel('Molecular Weight')
    plt.ylabel('Density')
    plt.title('Molecular Weight Distribution')
    plt.tight_layout()

    # Save plot
    save_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(save_path)


if __name__ == '__main__':
    from tap import Tap

    class Args(Tap):
        data_path: Path  # Path to a CSV file containing SMILES.
        save_path: Path  # Path to a PDF file where the plot of molecular weights will be saved.
        smiles_column: str = SMILES_COLUMN  # The name of the column in data_path containing SMILES.

    molecular_weight_distribution(**Args().parse_args().as_dict())
