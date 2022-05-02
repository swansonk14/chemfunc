"""Converts reaction SMARTS to images."""
from pathlib import Path

import pandas as pd
from rdkit.Chem import AllChem, Draw
from tap import Tap
from tqdm import tqdm


class Args(Tap):
    data_path: Path  # Path to CSV file containing reaction SMARTS.
    save_dir: Path  # Path to a directory where visualized molecules will be saved as images.
    smarts_column: str = 'smarts'  # Name of the column containing reaction SMARTS.


def visualize_reactions(data_path: Path,
                        save_dir: Path,
                        smarts_column: str = 'smiles') -> None:
    """Converts reaction SMARTS to images

    :param data_path: Path to CSV file containing reaction SMARTS.
    :param save_dir: Path to a directory where visualized molecules will be saved as images.
    :param smarts_column: Name of the column containing reaction SMARTS.
    """
    # Load data
    data = pd.read_csv(data_path)
    print(f'Data size = {len(data):,}')

    # Convert SMARTS to reactions
    reactions = [AllChem.ReactionFromSmarts(s) for s in tqdm(data[smarts_column], desc='SMARTS to reaction')]

    # Visualize reactions
    num_digits = len(str(len(reactions)))
    save_dir.mkdir(parents=True, exist_ok=True)

    for i, reaction in enumerate(tqdm(reactions, desc='Saving images')):
        d2d = Draw.MolDraw2DCairo(800, 300)
        d2d.DrawReaction(reaction, highlightByReactant=True)
        png = d2d.GetDrawingText()

        with open(save_dir / f'{str(i + 1).zfill(num_digits)}.png', 'wb+') as f:
            f.write(png)


if __name__ == '__main__':
    visualize_reactions(**Args().parse_args().as_dict())
