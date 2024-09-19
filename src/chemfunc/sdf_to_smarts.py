"""Converts molecules in SDF format to a CSV with SMARTS."""
from pathlib import Path

from chemfunc.constants import SMARTS_COLUMN
from chemfunc.convert_sdf import convert_sdf


def sdf_to_smarts(
        data_path: Path,
        save_path: Path,
        smarts_column: str = SMARTS_COLUMN,
        all_properties: bool = False,
        properties: list[str] | None = None,
        deduplicate: bool = False
) -> None:
    """Converts molecules in SDF format to a CSV with SMARTS.

    :param data_path: Path to an SDF file.
    :param smarts_column: The name of the column where SMARTS will be saved.
    :param all_properties: Whether to extract all properties from the SDF.
    :param properties: List of properties to extract from the SDF for each molecule.
    :param save_path: Path to a CSV file where SMARTS strings will be saved.
    :param deduplicate: Whether to deduplicate SMARTS.
    """
    convert_sdf(
        data_path=data_path,
        save_path=save_path,
        molecule_column=smarts_column,
        molecule_type='smarts',
        all_properties=all_properties,
        properties=properties,
        deduplicate=deduplicate,
    )
