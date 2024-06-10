"""Given a dataset of molecules, computes the nearest neighbor molecule in a second dataset using one of several similarity metrics."""
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd

from chemfunc.constants import SMILES_COLUMN
from chemfunc.molecular_similarities import get_similarity_function


def nearest_neighbor(
        data_path: Path,
        reference_data_path: Path,
        metric: Literal['tanimoto', 'mcs', 'tversky'] = 'tanimoto',
        save_path: Path | None = None,
        smiles_column: str = SMILES_COLUMN,
        reference_smiles_column: str | None = None,
        reference_name: str | None = None,
        match_valences: bool = False,
        ring_matches_ring_only: bool = False,
        complete_rings_only: bool = False,
        denominator: Literal['mol_1', 'mol_2', 'avg'] = 'mol_2',
) -> None:
    """Given a dataset, computes the nearest neighbor molecule by Tanimoto similarity in a second dataset.

    :param data_path: Path to CSV file containing data with SMILES whose neighbors are to be computed.
    :param reference_data_path: Path to CSV file containing reference SMILES which will be the neighbors of data_path.
    :param metric: Metric to use when computing similarity.
                   tanimoto: Tanimoto similarity using Morgan fingerprint.
                   mcs: Maximum common substructure as a proportion of the number of atoms in the reference molecule.
                   tversky: Tversky index with alpha = 0 and beta = 1, i.e., the proportion of reference substructures
                   in the molecule.
    :param save_path: Where the data with the neighbor info should be saved (defaults to data_path).
    :param smiles_column: Name of the column in data_path containing SMILES.
    :param reference_smiles_column: Name of the column in reference_data_path containing SMILES.
                                    If None, then smiles_column is used.
    :param reference_name: Name of the reference data when naming the new columns with neighbor info.
    :param match_valences: For MCS only, whether to match valences when computing the MCS.
    :param ring_matches_ring_only: For MCS only, whether to only match rings to rings when computing the MCS.
    :param complete_rings_only: For MCS only, whether to only match complete rings when computing the MCS.
    :param denominator: For MCS only, the denominator (molecule size) to use when computing the similarity.
                    mol_1: similarity = MCS_size / num_atoms_mol_1
                    mol_2: similarity = MCS_size / num_atoms_mol_2
                    avg: similarity = 0.5 * (MCS_size / num_atoms_mol_1 + MCS_size / num_atoms_mol_2)
    """
    # Set save path
    if save_path is None:
        save_path = data_path

    # Set reference smiles column
    if reference_smiles_column is None:
        reference_smiles_column = smiles_column

    print('Loading data')
    data = pd.read_csv(data_path)
    reference_data = pd.read_csv(reference_data_path)

    # Sort reference data and drop duplicates
    reference_data.drop_duplicates(subset=reference_smiles_column, inplace=True)
    reference_data.sort_values(by=reference_smiles_column, ignore_index=True, inplace=True)

    # Set up similarity function kwargs
    if metric == 'mcs':
        kwargs = {
            'match_valences': match_valences,
            'ring_matches_ring_only': ring_matches_ring_only,
            'complete_rings_only': complete_rings_only,
            'denominator': denominator,
        }
    else:
        kwargs = {}

    print(f'Computing similarities using {metric} metric')
    similarities = get_similarity_function(metric)(
        data[smiles_column],
        reference_data[reference_smiles_column],
        **kwargs
    )

    print('Finding minimum distance SMILES')
    prefix = f'{f"{reference_name}_" if reference_name is not None else ""}{metric}_'

    max_similarity_indices = np.argmax(similarities, axis=1)

    data[f'{prefix}nearest_neighbor'] = [
        reference_data[reference_smiles_column][max_similarity_index]
        for max_similarity_index in max_similarity_indices
    ]
    data[f'{prefix}nearest_neighbor_similarity'] = [
        similarities[i, max_similarity_index]
        for i, max_similarity_index in enumerate(max_similarity_indices)
    ]

    print('Saving')
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)
