"""Given a dataset of molecules, computes the nearest neighbor molecule in a second dataset using one of several similarity metrics."""
from itertools import product
from multiprocessing import Pool
from pathlib import Path
from typing import Literal, Iterable, Optional, Union

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.rdFMCS import FindMCS
from sklearnex import patch_sklearn
patch_sklearn()
from sklearn.metrics import pairwise_distances
from tqdm import tqdm

from chem_utils.constants import SMILES_COLUMN
from chem_utils.morgan_fingerprint import compute_fingerprints


def compute_pairwise_tanimoto_distances(mols_1: list[Union[str, Chem.Mol]],
                                        mols_2: Optional[list[Union[str, Chem.Mol]]] = None) -> np.ndarray:
    """
    Computes pairwise Tanimoto distances between the molecules in mols_1 and mols_2.

    :param mols_1: A list of molecules, either SMILES strings or RDKit molecules.
    :param mols_2: A list of molecules, either SMILES strings or RDKit molecules.
                   If None, copies mols_1 list.
    :return: A 2D numpy array of pairwise distances.
    """
    # Compute Morgan fingerprints
    fps_1 = np.array(compute_fingerprints(mols_1, fingerprint_type='morgan'), dtype=bool)
    fps_2 = np.array(compute_fingerprints(mols_2, fingerprint_type='morgan'), dtype=bool) if mols_2 is not None else fps_1

    # Compute pairwise Tanimoto distances
    tanimoto_distances = pairwise_distances(fps_1, fps_2, metric='jaccard', n_jobs=-1)

    return tanimoto_distances


def compute_mcs_size(mols: Iterable[Chem.Mol]) -> int:
    """
    Computes the size (number of atoms) of the maximum common substructure between molecules.

    :param mols: An iterable of molecules.
    :return: The size (number of atoms) of the maximum common substructure between molecules.
    """
    return FindMCS(mols).numAtoms


def compute_pairwise_mcs_similarities(mols_1: list[Union[str, Chem.Mol]],
                                      mols_2: Optional[list[Union[str, Chem.Mol]]] = None) -> np.ndarray:
    """
    Computes pairwise maximum common substructure (MCS) similarities between the molecules in mols_1 and mols_2.

    :param mols_1: A list of molecules, either SMILES strings or RDKit molecules.
    :param mols_2: A list of molecules, either SMILES strings or RDKit molecules.
                   If None, copies mols_1 list.
    :return: A 2D numpy array of pairwise similarities.
    """
    # Convert SMILES to RDKit molecules if needed
    mols_1 = [Chem.MolFromSmiles(mol) if isinstance(mol, str) else mol for mol in mols_1]

    if mols_2 is not None:
        mols_2 = [Chem.MolFromSmiles(mol) if isinstance(mol, str) else mol for mol in mols_2]
    else:
        mols_2 = mols_1

    # Compute pairwise MCS similarities
    with Pool() as pool:
        pairwise_mcs = np.array(list(tqdm(pool.imap(compute_mcs_size, product(mols_1, mols_2)),
                                          total=len(mols_1) * len(mols_2))))

    pairwise_mcs = pairwise_mcs.reshape(len(mols_1), len(mols_2))

    num_atoms_2 = np.array([mol.GetNumAtoms() for mol in mols_2])
    mcs_similarities = pairwise_mcs / num_atoms_2

    return mcs_similarities


def compute_pairwise_tversky_similarities(mols_1: list[Union[str, Chem.Mol]],
                                          mols_2: Optional[list[Union[str, Chem.Mol]]] = None) -> np.ndarray:
    """
    Computes pairwise Tversky similarities between the molecules in mols_1 and mols_2.

    Uses alpha = 0 and beta = 1 so that the similarity is the proportion of substructures in each reference
    molecule (from mols_2) that is present in the query molecule (from mols_1).

    :param mols_1: A list of molecules, either SMILES strings or RDKit molecules.
    :param mols_2: A list of molecules, either SMILES strings or RDKit molecules.
                   If None, copies mols_1 list.
    :return: A 2D numpy array of pairwise similarities.
    """
    # Compute Morgan fingerprints
    fps_1 = np.array(compute_fingerprints(mols_1, fingerprint_type='morgan'), dtype=int)
    fps_2 = np.array(compute_fingerprints(mols_2, fingerprint_type='morgan'), dtype=int) if mols_2 is not None else fps_1

    # Compute pairwise Tversky similarities
    intersection = fps_1 @ fps_2.transpose()
    size_2 = fps_2.sum(axis=1)
    tversky_similarities = intersection / size_2

    return tversky_similarities


def add_nearest_neighbors(data: pd.DataFrame,
                          similarities: np.ndarray,
                          reference_smiles: list[str],
                          prefix: str = '') -> None:
    """
    Adds nearest neighbors to a DataFrame.

    :param data: The Pandas DataFrame to which the nearest neighbors will be added.
    :param similarities: A NumPy matrix of similarities between the data SMILES (rows)
                         and the reference SMILES (columns).
    :param reference_smiles: The reference SMILES corresponding to the columns of similarities.
    :param prefix: The prefix to describe the nearest neighbors.
    """
    assert similarities.shape[1] == len(reference_smiles)

    max_similarity_indices = np.argmax(similarities, axis=1)

    data[f'{prefix}nearest_neighbor'] = [
        reference_smiles[max_similarity_index] for max_similarity_index in max_similarity_indices
    ]
    data[f'{prefix}nearest_neighbor_similarity'] = [
        similarities[i, max_similarity_index] for i, max_similarity_index in enumerate(max_similarity_indices)
    ]


def nearest_neighbor_tanimoto(data_path: Path,
                              reference_data_path: Path,
                              metrics: Literal['tanimoto', 'mcs', 'tversky'] = 'tanimoto',
                              save_path: Optional[Path] = None,
                              smiles_column: str = SMILES_COLUMN,
                              reference_smiles_column: Optional[str] = None,
                              reference_name: Optional[str] = None):
    """Given a dataset, computes the nearest neighbor molecule by Tanimoto similarity in a second dataset.

    :param data_path: Path to CSV file containing data with SMILES whose neighbors are to be computed.
    :param reference_data_path: Path to CSV file containing reference SMILES which will be the neighbors of data_path.
    :param metrics: Metrics to use when computing similarity.
                    tanimoto: Tanimoto similarity using Morgan fingerprint.
                    mcs: Maximum common substructure as a proportion of the number of atoms in the reference molecule.
                    tversky: Tversky index with alpha = 0 and beta = 1, i.e., the proportion of reference substructures
                    in the molecule.
    :param save_path: Where the data with the neighbor info should be saved (defaults to data_path).
    :param smiles_column: Name of the column in data_path containing SMILES.
    :param reference_smiles_column: Name of the column in reference_data_path containing SMILES.
                                    If None, then smiles_column is used.
    :param reference_name: Name of the reference data when naming the new columns with neighbor info.
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

    for metric in metrics:
        print(f'Computing similarities using {metric} metric')
        if metric == 'tanimoto':
            similarities = 1 - compute_pairwise_tanimoto_distances(
                mols_1=data[smiles_column],
                mols_2=reference_data[reference_smiles_column]
            )
        elif metric == 'mcs':
            similarities = compute_pairwise_mcs_similarities(
                mols_1=data[smiles_column],
                mols_2=reference_data[reference_smiles_column]
            )
        elif metric == 'tversky':
            similarities = compute_pairwise_tversky_similarities(
                mols_1=data[smiles_column],
                mols_2=reference_data[reference_smiles_column]
            )
        else:
            raise ValueError(f'Similarity metric "{metric}" is not supported.')

        print('Finding minimum distance SMILES')
        prefix = f'{f"{reference_name}_" if reference_name is not None else ""}{metric}_'
        add_nearest_neighbors(
            data=data,
            similarities=similarities,
            reference_smiles=reference_data[reference_smiles_column],
            prefix=prefix
        )

    print('Saving')
    save_path.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(save_path, index=False)


if __name__ == '__main__':
    from tap import Tap

    class Args(Tap):
        data_path: Path  # Path to CSV file containing data with SMILES whose neighbors are to be computed.
        reference_data_path: Path  # Path to CSV file containing reference SMILES which will be the neighbors of data_path.
        metrics: list[Literal['tanimoto', 'mcs', 'tversky']] = ['tanimoto']  # Metrics to use when computing similarity.
        """
        tanimoto: Tanimoto similarity using Morgan fingerprint.
        mcs: Maximum common substructure as a proportion of the number of atoms in the reference molecule.
        tversky: Tversky index with alpha = 0 and beta = 1, i.e., the proportion of reference substructures in the molecule.
        """
        save_path: Optional[Path] = None  # Where the data with the neighbor info should be saved (defaults to data_path).
        smiles_column: str = SMILES_COLUMN  # Name of the column in data_path containing SMILES.
        reference_smiles_column: Optional[str] = None  # Name of the column in reference_data_path containing SMILES.
        """If None, then smiles_column is used."""
        reference_name: Optional[str] = None  # Name of the reference data when naming the new columns with neighbor info.

    nearest_neighbor_tanimoto(**Args().parse_args().as_dict())
