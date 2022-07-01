"""Holds constant values."""
from typing import Union

from rdkit import Chem

Molecule = Union[str, Chem.Mol]

CLUSTER_COLUMN = 'cluster_label'
SMILES_COLUMN = 'smiles'
