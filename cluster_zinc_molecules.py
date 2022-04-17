"""Clusters ZINC molecules."""

from pathlib import Path

from matplotlib import offsetbox
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from tap import Tap
from tqdm import tqdm

from constants import (
    ACTIVITY_COLUMN,
    ZINC_SMILES_COLUMN
)
from utils import compute_morgan_fingerprints


class Args(Tap):
    data_path: Path  # Path to CSV file containing SMILES strings and zinc indices.
    smiles_column: str = ZINC_SMILES_COLUMN  # SMILES column.
    num_clusters: int = 50  # Number of clusters.
    mols_per_row: int = 10  # Number of molecules in each row.
    num_rows: int = 5  # Maximum number of rows in images of each cluster.
    pred_column: str = ACTIVITY_COLUMN  # Column containing the predictions (will be sorted from highest to lowest).
    save_dir: Path  # Path to directory where the clustering and PNG images of molecules will be saved.

    def process_args(self) -> None:
        self.save_dir.mkdir(parents=True, exist_ok=True)


def cluster_zinc_molecules(args: Args) -> None:
    """Clusters ZINC molecules."""
    print('Loading data')
    data = pd.read_csv(args.data_path)
    print(f'Data size = {len(data):,}')
    data['mol'] = [Chem.MolFromSmiles(smiles) for smiles in data[args.smiles_column]]

    print('Sorting data')
    data.sort_values(by=args.pred_column, ascending=False, ignore_index=True, inplace=True)

    print('Computing Morgan fingerprints')
    morgans = compute_morgan_fingerprints(data['mol'])

    print('Clustering')
    kmeans = KMeans(n_clusters=args.num_clusters, random_state=0).fit(morgans)
    data['cluster_label'] = kmeans.labels_ + 1

    save_columns = list(data.columns)
    save_columns.remove('mol')

    print('Saving data')
    data.to_csv(args.save_dir / 'data.csv', columns=save_columns, index=False)

    print('Saving clusters')
    for cluster_label in tqdm(sorted(data['cluster_label'].unique())):
        # Save cluster data
        cluster_data = data[data['cluster_label'] == cluster_label]
        cluster_data.to_csv(args.save_dir / f'cluster_{cluster_label}.csv', columns=save_columns, index=False)

        # Only display top images
        cluster_data = cluster_data.iloc[:args.num_rows * args.mols_per_row]

        # Save cluster images
        img = Draw.MolsToGridImage(
            cluster_data['mol'],
            molsPerRow=args.mols_per_row,
            subImgSize=(600, 600),
            legends=[f'Score = {pred:.4f}' for pred in cluster_data[args.pred_column]]
        )
        img.save(args.save_dir / f'cluster_{cluster_label}.png')

    print('Saving top molecules')
    top_data = data.drop_duplicates(subset='cluster_label')  # mol with highest pred for each cluster since sorted
    top_data.to_csv(args.save_dir / 'top_molecules.csv', columns=save_columns, index=False)

    mols_per_row = int(np.ceil(np.sqrt(len(top_data))))
    img = Draw.MolsToGridImage(
        top_data['mol'],
        molsPerRow=mols_per_row,
        subImgSize=(600, 600),
        legends=[f'Cluster = {cluster}, Score = {pred:.4f}' for cluster, pred in zip(top_data['cluster_label'],
                                                                                     top_data[args.pred_column])]
    )
    img.save(args.save_dir / 'top_molecules.png', columns=save_columns, index=False)

    print('Running t-SNE')
    tsne = TSNE(n_components=2, init='pca', random_state=0, metric='jaccard')
    X = tsne.fit_transform(morgans)

    print('Plotting K-Means/t-SNE')
    x_min, x_max = np.min(X, axis=0), np.max(X, axis=0)
    X = (X - x_min) / (x_max - x_min)

    plt.figure(figsize=(6.4 * 40, 4.8 * 40))
    ax = plt.subplot(111)

    # Plot individual points
    for i, cluster_label in tqdm(enumerate(data['cluster_label'])):
        plt.text(X[i, 0], X[i, 1], str(cluster_label),
                 color=plt.cm.rainbow(cluster_label / args.num_clusters),
                 fontdict={'weight': 'bold', 'size': 80})

    # Show images of top molecules
    if hasattr(offsetbox, 'AnnotationBbox'):
        for i, smiles in tqdm(zip(top_data.index, top_data[args.smiles_column]), total=len(top_data)):
            img = Draw.MolsToGridImage([Chem.MolFromSmiles(smiles)], molsPerRow=1, subImgSize=(400, 400))
            imagebox = offsetbox.AnnotationBbox(offsetbox.OffsetImage(img, cmap=plt.cm.gray_r), X[i])
            ax.add_artist(imagebox)

    plt.xticks([]), plt.yticks([])
    plt.savefig(args.save_dir / 'kmeans_tsne.png')


if __name__ == '__main__':
    cluster_zinc_molecules(Args().parse_args())
