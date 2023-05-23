"""Converts molecule SMILES to grids of images."""
import math
from pathlib import Path
from typing import Literal

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from tqdm import trange, tqdm

from chemfunc.constants import SMILES_COLUMN


def visualize_molecules(
        data_path: Path,
        save_dir: Path,
        smiles_column: str = SMILES_COLUMN,
        num_rows: int = 5,
        mols_per_row: int = 10,
        num_molecules: int | None = None,
        image_format: Literal['PNG', 'SVG'] = 'PNG'
) -> None:
    """Converts molecule SMILES to grids of images.

    :param data_path: Path to CSV file containing SMILES.
    :param save_dir: Path to a directory where visualized molecules will be saved as images.
    :param smiles_column: Name of the column containing SMILES.
    :param num_rows: Number of rows of molecules/rationales per image.
    :param mols_per_row: Number of molecules/rationales per row.
    :param num_molecules: Number of molecules to visualize (if None, visualizes all molecules).
    :param image_format: Image format to use when saving images.
    """
    # Load data
    data = pd.read_csv(data_path)
    print(f'Data size = {len(data):,}')

    # Optionally subsample molecules
    if num_molecules is not None:
        data = data.iloc[:num_molecules]

    num_molecules = len(data)

    # Convert SMILES to Mols
    mols = [Chem.MolFromSmiles(s) for s in tqdm(data[smiles_column], desc='SMILES to Mol')]

    # Visualize molecules
    mols_per_image = num_rows * mols_per_row
    max_index = int(math.ceil(len(data) / mols_per_image))
    num_digits = len(str(max_index))
    print(f'Visualizing {num_molecules:,} molecules')

    save_dir.mkdir(parents=True, exist_ok=True)

    for i in trange(0, num_molecules, mols_per_image, desc='Saving images'):
        image_mols = mols[i:i + mols_per_image]

        img = Draw.MolsToGridImage(
            image_mols,
            molsPerRow=mols_per_row,
            subImgSize=(600, 600),
            legends=[f'Molecule {i + j + 1}' for j in range(len(image_mols))],
            useSVG=image_format == 'SVG'
        )

        save_path = save_dir / f'{str(i // mols_per_image + 1).zfill(num_digits)}.{image_format.lower()}'

        if image_format == 'SVG':
            with open(save_path, 'w') as f:
                f.write(img)
        elif image_format == 'PNG':
            img.save(save_path)
        else:
            raise ValueError(f'Image format "{image_format}" is not supported.')
