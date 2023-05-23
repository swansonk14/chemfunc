"""Converts a SMILES string to an SVG image of the molecule."""
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2DSVG


def smiles_to_svg(
        smiles: str,
        save_path: Path,
        width: int = 600,
        height: int = 600
) -> None:
    """Converts a SMILES string to an SVG image of the molecule.

    :param smiles: A SMILES string.
    :param save_path: Path to an SVG file where the molecule image will be saved.
    :param width: Width of the SVG image.
    :param height: Height of the SVG image.
    """
    # Convert SMILES to Mol
    mol = Chem.MolFromSmiles(smiles)

    # Convert Mol to SVG
    d = MolDraw2DSVG(width, height)
    d.DrawMolecule(mol)
    d.FinishDrawing()
    svg = d.GetDrawingText()

    # Save SVG
    save_path.parent.mkdir(parents=True, exist_ok=True)
    save_path.write_text(svg)
