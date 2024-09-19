"""Converts molecules in SDF format to a CSV with SMILES."""
from pathlib import Path

from chemfunc.constants import SMILES_COLUMN
from chemfunc.convert_sdf import convert_sdf


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
    convert_sdf(
        data_path=data_path,
        save_path=save_path,
        molecule_column=smiles_column,
        molecule_type='smiles',
        all_properties=all_properties,
        properties=properties,
        deduplicate=deduplicate,
    )
