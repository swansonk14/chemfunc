"""Converts reaction SMARTS to images."""
from pathlib import Path
from typing import Literal, Optional

import pandas as pd
from rdkit.Chem import AllChem, Draw
from tqdm import tqdm


def visualize_reactions(data_path: Path,
                        save_dir: Path,
                        save_format: Literal['png', 'svg'] = 'svg',
                        smarts_column: str = 'smarts',
                        name_column: Optional[str] = None) -> None:
    """Converts reaction SMARTS to images

    :param data_path: Path to CSV file containing reaction SMARTS.
    :param save_dir: Path to a directory where visualized molecules will be saved as images.
    :param save_format: Image format to use when saving images.
    :param smarts_column: Name of the column containing reaction SMARTS.
    :param name_column: Name of the column containing the reaction name to use when naming the image file.
    """
    # Load data
    data = pd.read_csv(data_path)
    print(f'Data size = {len(data):,}')

    # Convert SMARTS to reactions
    reactions = [AllChem.ReactionFromSmarts(s) for s in tqdm(data[smarts_column], desc='SMARTS to reaction')]

    # Get reaction names
    if name_column is not None:
        names = data[name_column]
    else:
        num_digits = len(str(len(reactions)))
        names = [str(i + 1).zfill(num_digits) for i in range(len(data))]

    # Visualize reactions
    save_dir.mkdir(parents=True, exist_ok=True)

    for i, (reaction, name) in enumerate(tqdm(zip(reactions, names), total=len(data), desc='Saving images')):
        if save_format == 'svg':
            d2d = Draw.MolDraw2DSVG(800, 300)
        else:
            d2d = Draw.MolDraw2DCairo(800, 300)

        d2d.DrawReaction(reaction, highlightByReactant=True)

        if save_format == 'svg':
            d2d.FinishDrawing()

        with open(save_dir / f'{name}.{save_format}', 'wb+' if save_format == 'png' else 'w') as f:
            f.write(d2d.GetDrawingText())


if __name__ == '__main__':
    from tap import Tap

    class Args(Tap):
        data_path: Path  # Path to CSV file containing reaction SMARTS.
        save_dir: Path  # Path to a directory where visualized molecules will be saved as images.
        save_format: Literal['png', 'svg'] = 'svg'  # Image format to use when saving images.
        smarts_column: str = 'smarts'  # Name of the column containing reaction SMARTS.
        name_column: Optional[str] = None  # Name of the column containing the reaction name to use when naming the image file.

    visualize_reactions(**Args().parse_args().as_dict())
