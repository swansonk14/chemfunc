"""Converts molecules in SDF format to a CSV with SMILES or SMARTS."""
from pathlib import Path
from typing import Literal

import pandas as pd
from rdkit import Chem
from tqdm import tqdm


def mol_to_properties(
        mol: Chem.Mol,
        molecule_column: str,
        molecule_type: Literal["smiles", "smarts"],
        all_properties: bool = False,
        properties: list[str] | None = None,
) -> dict:
    """Converts an RDKit molecule to a dictionary containing the molecule (SMILES or SMARTS) and optionally properties.

    :param mol: RDKit molecule.
    :param molecule_column: The name of the column where molecule (SMILES or SMARTS) will be saved.
    :param molecule_type: The type of molecule to save (either "smiles" or "smarts").
    :param all_properties: Whether to extract all properties from the SDF.
    :param properties: List of properties to extract from the SDF for each molecule.
    :return: Dictionary containing the molecule (SMILES or SMARTS) and optionally properties.
    """
    # Convert molecule to SMILES or SMARTS
    if molecule_type == "smiles":
        molecule = Chem.MolToSmiles(mol)
    elif molecule_type == "smarts":
        molecule = Chem.MolToSmarts(mol)
    else:
        raise ValueError(f"Invalid molecule_type: {molecule_type}")

    # Get properties
    if all_properties:
        props = mol.GetPropsAsDict()
    elif properties is not None:
        props = {prop: mol.GetProp(prop) for prop in properties}
    else:
        props = {}

    return {
        molecule_column: molecule,
        **props
    }


def convert_sdf(
        data_path: Path,
        save_path: Path,
        molecule_column: str,
        molecule_type: Literal["smiles", "smarts"],
        all_properties: bool = False,
        properties: list[str] | None = None,
        deduplicate: bool = False
) -> None:
    """Converts molecules in SDF format to a CSV with molecules (SMILES or SMARTS).

    :param data_path: Path to an SDF file.
    :param save_path: Path to a CSV file where the molecules (SMILES or SMARTS) will be saved.
    :param molecule_column: The name of the column where the molecules (SMILES or SMARTS) will be saved.
    :param molecule_type: The type of molecule to save (either "smiles" or "smarts").
    :param all_properties: Whether to extract all properties from the SDF.
    :param properties: List of properties to extract from the SDF for each molecule.
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
        mol_to_properties(
            mol=mol,
            molecule_column=molecule_column,
            molecule_type=molecule_type,
            all_properties=all_properties,
            properties=properties
        ) for mol in tqdm(mols, desc=f'Mol to {molecule_type}')
    ])

    # Optionally deduplicate
    if deduplicate:
        print(f'Deduplicating {molecule_type}')
        data.drop_duplicates(subset=molecule_column, inplace=True)

    # Print stats
    print(f'Data size = {len(data):,}')
    print(f'Number of unique {molecule_type} = {data[molecule_column].nunique():,}')

    for prop in sorted(set(data.columns) - {molecule_column}):
        print(f'Number of {prop} (# unique) = {data[prop].notna().sum():,} ({data[prop].nunique():,})')

    # Save data as CSV
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)
