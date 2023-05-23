"""Entry point to the chem_utils package."""
from tap import Tap

from chem_utils import *


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
    chem_utils_command: str = ''

    def configure(self):
        self.add_subparsers(help='chem_utils functions', dest='chem_utils_command')

        for function in FUNCTIONS:
            tap_class = convert_to_tap(function)
            self.add_subparser(function.__name__, tap_class, help=tap_class.__doc__)


def main():
    """Entry point to the chem_utils package."""
    args = ChemUtilsTap().parse_args()

    function = NAME_TO_FUNCTION[args.chem_utils_command]

    args_dict = args.as_dict()
    del args_dict['chem_utils_command']

    function(**args_dict)
