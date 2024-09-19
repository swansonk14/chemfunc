"""Canonicalizes SMILES using RDKit canonicalization and optionally strips salts."""
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from tqdm import tqdm

from chemfunc.constants import SMILES_COLUMN


def canonicalize_smiles(
        data_path: Path,
        save_path: Path | None = None,
        smiles_column: str = SMILES_COLUMN,
        canonical_smiles_column: str = SMILES_COLUMN,
        remove_salts: bool = False,
        delete_disconnected_mols: bool = False
) -> None:
    """Canonicalizes SMILES using RDKit canonicalization and optionally strips salts.

    :param data_path: Path to CSV file containing SMILES.
    :param save_path: Path where CSV file with canonicalized SMILES will be saved. If None, uses data_path.
    :param smiles_column: Name of the column containing SMILES.
    :param canonical_smiles_column: Name of the column where canonicalized SMILES will be saved.
    :param remove_salts: Whether to remove salts from the SMILES.
    :param delete_disconnected_mols: Whether to delete disconnected molecules, i.e., any molecule whose
                                     SMILES has a '.' in it. This is performed after (optionally) removing salts.
                                     This deletes the entire row from the Pandas DataFrame.
    """
    # Load data
    data = pd.read_csv(data_path)
    print(f'Data size = {len(data):,}')

    # Convert SMILES to mol
    mols = [Chem.MolFromSmiles(smiles) for smiles in tqdm(data[smiles_column], desc='SMILES to mol')]

    # Handle SMILES that RDKit cannot process
    valid_mols = [mol is not None for mol in mols]

    if not all(valid_mols):
        print(f'Found {len(valid_mols) - sum(valid_mols)} invalid molecules. Deleting.')
        data = data[valid_mols]
        mols = [mol for mol, valid_mol in zip(mols, valid_mols) if valid_mol]

    # Optionally remove salts
    if remove_salts:
        remover = SaltRemover()
        mols = [remover.StripMol(mol, dontRemoveEverything=True) for mol in tqdm(mols, desc='Stripping salts')]

    # Convert mol to SMILES (implicitly canonicalizes SMILES)
    data[canonical_smiles_column] = [Chem.MolToSmiles(mol) for mol in tqdm(mols, desc='Mol to SMILES')]

    # Optionally delete disconnected molecules
    if delete_disconnected_mols:
        connected_mols = ['.' not in smiles for smiles in data[canonical_smiles_column]]

        if not all(connected_mols):
            print(f'Found {len(connected_mols) - sum(connected_mols)} disconnected molecules. Deleting.')
            data = data[connected_mols]
            print(f'Remaining data size = {len(data):,}')

    # Save data
    if save_path is None:
        save_path = data_path

    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)
