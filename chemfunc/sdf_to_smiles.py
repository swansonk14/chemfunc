"""Converts molecules in SDF format to a CSV with SMILES."""
from pathlib import Path

import pandas as pd
from rdkit import Chem
from tqdm import tqdm

from chemfunc.constants import SMILES_COLUMN


def mol_to_smiles_with_properties(
        mol: Chem.Mol,
        smiles_column: str = SMILES_COLUMN,
        all_properties: bool = False,
        properties: list[str] | None = None,
) -> dict:
    """Converts an RDKit molecule to a dictionary containing SMILES and optionally properties.

    :param mol: RDKit molecule.
    :param smiles_column: The name of the column where SMILES will be saved.
    :param all_properties: Whether to extract all properties from the SDF.
    :param properties: List of properties to extract from the SDF for each molecule.
    :return: Dictionary containing SMILES and optionally properties.
    """
    # Get properties
    if all_properties:
        props = mol.GetPropsAsDict()
    elif properties is not None:
        props = {prop: mol.GetProp(prop) for prop in properties}
    else:
        props = {}

    return {
        smiles_column: Chem.MolToSmiles(mol),
        **props
    }


def sdf_to_smiles(
        data_path: Path,
        save_path: Path,
        smiles_column: str = SMILES_COLUMN,
        all_properties: bool = False,
        properties: list[str] | None = None,
        deduplicate: bool = False
) -> None:
    """Converts molecules in SDF format to a CSV with SMILES.

    :param data_path: Path to an SDF file.
    :param smiles_column: The name of the column where SMILES will be saved.
    :param all_properties: Whether to extract all properties from the SDF.
    :param properties: List of properties to extract from the SDF for each molecule.
    :param save_path: Path to a CSV file where SMILES strings will be saved.
    :param deduplicate: Whether to deduplicate SMILES.
    """
    # Validate arguments
    if all_properties and properties is not None:
        raise ValueError('Cannot specify both all_properties and properties')

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

    # Covert mols to SMILES and optionally extract properties
    data = pd.DataFrame([
        mol_to_smiles_with_properties(
            mol=mol,
            smiles_column=smiles_column,
            all_properties=all_properties,
            properties=properties
        ) for mol in tqdm(mols, desc='Mol to SMILES')
    ])

    # Optionally deduplicate
    if deduplicate:
        print('Deduplicating SMILES')
        data.drop_duplicates(subset=smiles_column, inplace=True)

    # Print stats
    print(f'Data size = {len(data):,}')
    print(f'Number of unique smiles = {data[smiles_column].nunique():,}')

    for prop in sorted(set(data.columns) - {smiles_column}):
        print(f'Number of {prop} (# unique) = {data[prop].notna().sum():,} ({data[prop].nunique():,})')

    # Save data as CSV
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)
