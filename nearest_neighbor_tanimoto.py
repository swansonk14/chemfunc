""" Given a dataset, computes the nearest neighbor by Tanimoto similarity to molecules in a second dataset. """

from typing import List

import numpy as np
import pandas as pd
from tap import Tap



class Args(Tap):
    data_path: str  # Path to CSV file containing data with SMILES whose neighbors are to be computed.
    smiles_column: str = 'smiles'  # Name of the column in data_path containing SMILES.
    include_activity: bool = False  # Whether to add a column indicating the activity of the nearest neighbor.
    nearest_active: bool = False  # Whether to add columns with the SMILES and activity of the nearest active neighbor.
    reference_data_path: str  # Path to CSV file containing reference SMILES which will be the neighbors of data_path.
    reference_smiles_column: str = None  # Name of the column in reference_data_path containing SMILES.
    """If None, then smiles_column is used."""
    reference_target_columns: List[str] = None  # Names of the columns in reference_data_path containing targets.
    """If None, then all columns except reference_smiles_column will be used."""
    reference_name: str = None  # Name of the reference data for use in naming the new columns with neighbor info.
    save_path: str = None  # Where the data with the neighbor info should be saved (defaults to data_path).

    def process_args(self) -> None:
        if self.reference_smiles_column is None:
            self.reference_smiles_column = self.smiles_column

        if self.save_path is None:
            self.save_path = self.data_path


def compute_morgan_fingerprint(mol: Union[str, Chem.Mol],
                               radius: int = 2,
                               num_bits: int = 2048) -> np.ndarray:
    """
    Generates a binary Morgan fingerprint for a molecule.

    :param mol: A molecule (i.e., either a SMILES string or an RDKit molecule).
    :param radius: Morgan fingerprint radius.
    :param num_bits: Number of bits in Morgan fingerprint.
    :return: A 1D numpy array containing the binary Morgan fingerprint.
    """
    mol = Chem.MolFromSmiles(mol) if type(mol) == str else mol
    morgan_vec = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=num_bits)
    morgan_fp = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(morgan_vec, morgan_fp)
    morgan_fp = morgan_fp.astype(bool)

    return morgan_fp


def compute_morgan_fingerprints(mols: List[Union[str, Chem.Mol]],
                                radius: int = 2,
                                num_bits: int = 2048) -> np.ndarray:
    """
    Generates binary Morgan fingerprints for each molecule in a list of molecules (in parallel).

    :param mols: A list of molecules (i.e., either a SMILES string or an RDKit molecule).
    :param radius: Morgan fingerprint radius.
    :param num_bits: Number of bits in Morgan fingerprint.
    :return: A 2D numpy array containing the binary Morgan fingerprints.
    """
    morgan_fn = partial(compute_morgan_fingerprint, radius=radius, num_bits=num_bits)

    with Pool() as pool:
        morgan_fps = list(tqdm(pool.imap(morgan_fn, mols), total=len(mols)))

    return np.array(morgan_fps)


def compute_pairwise_tanimoto_distances(mols_1: List[Union[str, Chem.Mol]],
                                        mols_2: List[Union[str, Chem.Mol]]) -> np.ndarray:
    """
    Computes pairwise Tanimoto distances between the molecules in :attr:`mols_1` and :attr:`mols_1`.

    :param mols_1: A list of molecules, either SMILES strings or RDKit molecules.
    :param mols_2: A list of molecules, either SMILES strings or RDKit molecules.
    :param parallel: Whether to run the computation in parallel across multiple cores.
    :return: A 2D numpy array of pairwise distances.
    """
    # Compute Morgan fingerprints
    fps_1 = np.array(compute_morgan_fingerprints(mols_1), dtype=bool)
    fps_2 = np.array(compute_morgan_fingerprints(mols_2), dtype=bool)

    # Compute pairwise distances
    tanimoto_distances = pairwise_distances(fps_1, fps_2, metric='jaccard', n_jobs=-1)

    return tanimoto_distances


def add_nearest_neighbors(data: pd.DataFrame,
                          similarities: np.ndarray,
                          reference_smiles: List[str],
                          prefix: str = '') -> None:
    """
    Adds nearest neighbors to a DataFrame.

    :param data: The Pandas DataFrame to which the nearest neighbors will be added.
    :param similarities: A NumPy matrix of Tanimoto similarities between the data SMILES (rows)
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


def add_nearest_neighbor_activity(data: pd.DataFrame,
                                  reference_data: pd.DataFrame,
                                  reference_target_columns: List[str],
                                  prefix: str = '') -> None:
    """
    Adds activity information for the nearest neighbors.

    :param data: The Pandas DataFrame to which the nearest neighbor activity will be added.
    :param reference_data: The Pandas DataFrame containing the reference data, indexed by SMILES.
    :param reference_target_columns: The columns in reference_data containing targets.
    :param prefix: The prefix to describe the nearest neighbors.
    """
    data[f'{prefix}nearest_neighbor_activity'] = [
        '|'.join(
            target
            for target in reference_target_columns
            if reference_data.loc[nearest_neighbor][target] == 1
        )
        for nearest_neighbor in data[f'{prefix}nearest_neighbor']
    ]


def nearest_neighbor_tanimoto(args: Args):
    """Given a dataset, computes the nearest neighbor by Tanimoto similarity to molecules in a second dataset."""
    print('Loading data')
    data = pd.read_csv(args.data_path)
    reference_data = pd.read_csv(args.reference_data_path)

    # Sort reference data and drop duplicates
    reference_data.drop_duplicates(subset=args.reference_smiles_column, inplace=True)
    reference_data.sort_values(by=args.reference_smiles_column, ignore_index=True, inplace=True)

    print('Computing Morgan fingerprints')
    similarities = 1 - compute_pairwise_tanimoto_distances(
        mols_1=data[args.smiles_column],
        mols_2=reference_data[args.reference_smiles_column]
    )

    print('Finding minimum distance SMILES')
    prefix = f'{args.reference_name}_' if args.reference_name is not None else ''
    add_nearest_neighbors(
        data=data,
        similarities=similarities,
        reference_smiles=reference_data[args.reference_smiles_column],
        prefix=prefix
    )

    if args.include_activity or args.nearest_active:
        reference_data.set_index(args.reference_smiles_column, inplace=True)
        reference_target_columns = sorted(
            args.reference_target_columns if args.reference_target_columns is not None else reference_data.columns
        )

        if args.include_activity:
            add_nearest_neighbor_activity(
                data=data,
                reference_data=reference_data,
                reference_target_columns=reference_target_columns,
                prefix=prefix
            )

        if args.nearest_active:
            reference_activity = reference_data[reference_target_columns].any(axis=1)

            active_reference_smiles = reference_data.index[reference_activity]
            active_similarities = similarities[:, reference_activity]
            active_prefix = f'{prefix}active_'

            add_nearest_neighbors(
                data=data,
                similarities=active_similarities,
                reference_smiles=active_reference_smiles,
                prefix=active_prefix
            )

            if args.include_activity:
                add_nearest_neighbor_activity(
                    data=data,
                    reference_data=reference_data,
                    reference_target_columns=reference_target_columns,
                    prefix=active_prefix
                )

    print('Saving')
    makedirs(args.save_path, isfile=True)
    data.to_csv(args.save_path, index=False)


if __name__ == '__main__':
    nearest_neighbor_tanimoto(Args().parse_args())
