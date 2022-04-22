"""Converts molecule SMILES to grids of images."""
import math
from pathlib import Path
from typing import Optional

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from tap import Tap
from tqdm import trange, tqdm


class Args(Tap):
    data_path: Path  # Path to CSV file containing SMILES.
    save_dir: Path  # Path to a directory where visualized molecules will be saved as images.
    smiles_column: str = 'smiles'  # Name of the column containing SMILES.
    num_rows: int = 4  # Number of rows of molecules/rationales per image.
    mols_per_row: int = 8  # Number of molecules/rationales per row.
    num_molecules: Optional[int] = None  # Number of molecules to visualize (if None, visualizes all molecules).


def visualize_molecules(data_path: Path,
                        save_dir: Path,
                        smiles_column: str = 'smiles',
                        num_rows: int = 4,
                        mols_per_row: int = 8,
                        num_molecules: Optional[int] = None) -> None:
    """Converts molecule SMILES to grids of images.

    :param data_path: Path to CSV file containing SMILES.
    :param save_dir: Path to a directory where visualized molecules will be saved as images.
    :param smiles_column: Name of the column containing SMILES.
    :param num_rows: Number of rows of molecules/rationales per image.
    :param mols_per_row: Number of molecules/rationales per row.
    :param num_molecules: Number of molecules to visualize (if None, visualizes all molecules).
    """
    # Load data
    data = pd.read_csv(data_path)
    print(f'Data size = {len(data):,}')

    # Convert SMILES to Mols
    mols = [Chem.MolFromSmiles(s) for s in tqdm(data[smiles_column], desc='SMILES to Mol')]

    # Visualize molecules
    mols_per_image = num_rows * mols_per_row
    max_index = int(math.ceil(len(data) / mols_per_image))
    num_digits = len(str(max_index))
    num_molecules = num_molecules or len(data)
    print(f'Visualizing {num_molecules:,} molecules')

    save_dir.mkdir(parents=True, exist_ok=True)

    for i in trange(0, num_molecules, mols_per_image, desc='Saving images'):
        image_mols = mols[i:i + mols_per_image]

        img = Draw.MolsToGridImage(
            image_mols,
            molsPerRow=mols_per_row,
            subImgSize=(600, 600),
            legends=[f'Molecule {i + j + 1}' for j in range(len(image_mols))]
        )
        img.save(save_dir / f'{str(i // mols_per_image + 1).zfill(num_digits)}.png')


if __name__ == '__main__':
    visualize_molecules(**Args().parse_args().as_dict())
