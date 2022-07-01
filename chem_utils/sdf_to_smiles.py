"""Converts molecules in SDF format to a CSV with SMILES."""
from pathlib import Path
from typing import Optional

import pandas as pd
from rdkit import Chem
from tqdm import tqdm

from chem_utils.constants import SMILES_COLUMN


def sdf_to_smiles(data_path: Path,
                  save_path: Path,
                  properties: Optional[set[str]] = None) -> None:
    """Converts molecules in SDF format to a CSV with SMILES.

    :param data_path: Path to an SDF file.
    :param properties: Set of properties to extract from the SDF for each molecule.
    :param save_path: Path to a CSV file where SMILES strings will be saved.
    """
    # Default to empty set for properties
    if properties is None:
        properties = set()

    # Load SDF file
    num_skipped = 0
    mols = []

    with Chem.SDMolSupplier(str(data_path)) as suppl:
        for mol in tqdm(suppl, desc='Loading SDF'):
            if mol is None:
                num_skipped += 1
            else:
                mols.append(mol)

    print(f'Number of molecules = {len(mols):,}')
    print(f'Number skipped = {num_skipped:,}')

    # Put data in Pandas DataFrame
    data = pd.DataFrame(data=[
        {
            SMILES_COLUMN: Chem.MolToSmiles(mol),
            **{
                prop: mol.GetProp(prop)
                for prop in properties
            }
        } for mol in tqdm(mols, desc='Mol to SMILES')
    ])

    # Print stats
    print(f'Data size = {len(data):,}')
    print(f'Number of unique smiles = {data[SMILES_COLUMN].nunique():,}')

    for prop in properties:
        print(f'Number of unique {prop} = {data[prop].nunique():,}')

    # Save data as CSV
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)


if __name__ == '__main__':
    from tap import Tap

    class Args(Tap):
        data_path: Path  # Path to an SDF file.
        save_path: Path  # Path to a CSV file where SMILES strings will be saved.
        properties: Optional[set[str]] = None  # Set of properties to extract from the SDF for each molecule.

    sdf_to_smiles(**Args().parse_args().as_dict())
