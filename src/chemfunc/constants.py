"""Holds constant values."""
from rdkit.Chem import Mol

Molecule = str | Mol

CLUSTER_COLUMN = 'cluster_label'
SMARTS_COLUMN = 'smarts'
SMILES_COLUMN = 'smiles'
