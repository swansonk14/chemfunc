"""Converts a SMILES string to an SVG image of the molecule."""
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2DSVG


def smiles_to_svg(smiles: str, save_path: Path) -> None:
    """Converts a SMILES string to an SVG image of the molecule.

    :param smiles: A SMILES string.
    :param save_path: Path to an SVG file where the molecule image will be saved.
    """
    # Convert SMILES to Mol
    mol = Chem.MolFromSmiles(smiles)

    # Convert Mol to SVG
    d = MolDraw2DSVG(600, 600)
    d.DrawMolecule(mol)
    d.FinishDrawing()
    svg = d.GetDrawingText()

    # Save SVG
    save_path.parent.mkdir(parents=True, exist_ok=True)
    save_path.write_text(svg)


if __name__ == '__main__':
    from tap import tapify

    tapify(smiles_to_svg)
