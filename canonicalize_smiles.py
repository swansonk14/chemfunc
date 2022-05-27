"""Canonicalizes SMILES using RDKit canonicalization and optionally strips salts."""
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from tap import Tap
from tqdm import tqdm

from constants import SMILES_COLUMN


class Args(Tap):
    data_path: Path  # Path to CSV file containing SMILES.
    save_path: Path  # Path where CSV file with canonicalized SMILES will be saved.
    smiles_column: str = SMILES_COLUMN  # Name of the column containing SMILES.
    remove_salts: bool = False  # Whether to remove salts from the SMILES.

    def process_args(self) -> None:
        self.save_path.parent.mkdir(parents=True, exist_ok=True)


def canonicalize_smiles(data_path: Path,
                        save_path: Path,
                        smiles_column: str = SMILES_COLUMN,
                        remove_salts: bool = False) -> None:
    """Canonicalizes SMILES using RDKit canonicalization and optionally strips salts.

    :param data_path: Path to CSV file containing SMILES.
    :param save_path: Path where CSV file with canonicalized SMILES will be saved.
    :param smiles_column: Name of the column containing SMILES.
    :param remove_salts: Whether to remove salts from the SMILES.
    """
    # Load data
    data = pd.read_csv(data_path)

    # Convert SMILES to mol
    mols = [Chem.MolFromSmiles(smiles) for smiles in tqdm(data[smiles_column], desc='SMILES to mol')]

    # Optionally remove salts
    if remove_salts:
        remover = SaltRemover()
        mols = [remover.StripMol(mol, dontRemoveEverything=True) for mol in tqdm(mols, desc='Stripping salts')]

    # Convert mol to SMILES (implicitly canonicalizes SMILES)
    data[smiles_column] = [Chem.MolToSmiles(mol) for mol in tqdm(mols, desc='Mol to SMILES')]

    # Save data
    data.to_csv(save_path, index=False)


if __name__ == '__main__':
    canonicalize_smiles(**Args().parse_args().as_dict())
