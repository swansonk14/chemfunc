"""Convert molecule SMILES to grids of PNG images."""
import math
from pathlib import Path
from typing import Optional

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from tap import Tap
from tqdm import trange

from constants import (
    ACTIVITY_COLUMN,
    CANONICAL_SMILES_COLUMN,
    MEAN_COLUMN
)


class Args(Tap):
    data_path: Path  # Path to CSV file containing SMILES.
    smiles_column: str = 'smiles'  # Name of the column containing SMILES.
    value_column: Optional[str] = None  # Name of the column containing values to display.
    sort_by_value: bool = False  # Whether to sort the molecules by value from highest to lowest.
    activity: int = None  # If provided, only displays molecules with that activity.
    save_dir: Path  # Path to directory where where visualized rationales will be saved as PNGs.
    num_rows: int = 4  # Number of rows of molecules/rationales per image.
    mols_per_row: int = 8  # Number of molecules/rationales per row.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)

        if self.value_column is None:
            self.value_column = self.smiles_column


def visualize_molecules(args: Args) -> None:
    """Convert molecule SMILES to grids of PNG images."""
    # Load data
    data = pd.read_csv(args.data_path)
    print(f'Data size = {len(data):,}')

    if args.activity is not None:
        data = data[data[ACTIVITY_COLUMN] == args.activity]
        print(f'Data size after filtering by activity {args.activity}: {len(data):,}')

    if args.sort_by_value:
        data.sort_values(by=args.value_column, ascending=False, inplace=True)

    smiles = list(data[args.smiles_column])
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    values = list(data[args.value_column])

    # Visualize molecules
    mols_per_image = args.num_rows * args.mols_per_row
    max_index = int(math.ceil(len(data) / mols_per_image)) - 1
    num_digits = len(str(max_index))

    for i in trange(0, len(data), mols_per_image):
        img = Draw.MolsToGridImage(
            mols[i:i + mols_per_image],
            molsPerRow=args.mols_per_row,
            subImgSize=(600, 600),
            legends=[f'Molecule {i + j + 1}: {value:.4f}' for j, value in enumerate(values[i:i + mols_per_image])]
        )
        img.save(args.save_dir / f'{str(i // mols_per_image).zfill(num_digits)}.png')


if __name__ == '__main__':
    visualize_molecules(Args().parse_args())
