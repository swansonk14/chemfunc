"""Entry point to the chemfunc package."""
from tap import Tap

from chemfunc import *


FUNCTIONS = [
    canonicalize_smiles,
    chemical_diversity,
    cluster_molecules,
    compute_properties,
    deduplicate_smiles,
    filter_molecules,
    measure_experimental_reproducibility,
    nearest_neighbor,
    plot_property_distribution,
    plot_tsne,
    regression_to_classification,
    sample_molecules,
    save_fingerprints,
    sdf_to_smiles,
    select_from_clusters,
    smiles_to_svg,
    visualize_molecules,
    visualize_reactions
]
NAME_TO_FUNCTION = {
    function.__name__: function
    for function in FUNCTIONS
}


class ChemUtilsTap(Tap):
    chemfunc_command: str = ''

    def configure(self):
        self.add_subparsers(help='chemfunc functions', dest='chemfunc_command')

        for function in FUNCTIONS:
            tap_class = convert_to_tap(function)
            self.add_subparser(function.__name__, tap_class, help=tap_class.__doc__)


def main():
    """Entry point to the chemfunc package."""
    args = ChemUtilsTap().parse_args()

    function = NAME_TO_FUNCTION[args.chemfunc_command]

    args_dict = args.as_dict()
    del args_dict['chemfunc_command']

    function(**args_dict)
