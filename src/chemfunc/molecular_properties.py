"""Contains functions to compute molecular properties."""
from functools import wraps
from typing import Callable

from rdkit import Chem
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.QED import qed
from rdkit.Chem.rdMolDescriptors import CalcTPSA

from chemfunc.constants import Molecule

PropertyFunction = Callable[[Molecule], float]
PROPERTY_FUNCTION_REGISTRY = {}


def register_property_function(property_type: str) -> Callable[[PropertyFunction], PropertyFunction]:
    """Creates a decorator which registers a property function in a global dictionary to enable access by name.

    :param property_type: The name to use to access the property function.
    :return: A decorator which will add a property function to the registry using the specified name.
    """

    def decorator(property_function: PropertyFunction) -> PropertyFunction:
        PROPERTY_FUNCTION_REGISTRY[property_type] = property_function
        return property_function

    return decorator


def get_property_function(property_type: str) -> PropertyFunction:
    """Gets a registered property function by name.

    :param property_type: The name of the property function.
    :return: The desired property function.
    """
    if property_type not in PROPERTY_FUNCTION_REGISTRY:
        raise ValueError(f'Property function "{property_type}" could not be found.')

    return PROPERTY_FUNCTION_REGISTRY[property_type]


def get_available_property_functions() -> list[str]:
    """Returns a list of names of available property functions."""
    return sorted(PROPERTY_FUNCTION_REGISTRY)


def smiles_to_mol_wrapper(property_function: PropertyFunction) -> PropertyFunction:
    """A decorator which converts a SMILES to an RDKit molecule before passing it to a property function."""

    @wraps(property_function)
    def wrapper(molecule: Molecule) -> float:
        if isinstance(molecule, str):
            molecule = Chem.MolFromSmiles(molecule)

        return property_function(molecule)

    return wrapper


@register_property_function('clogp')
@smiles_to_mol_wrapper
def compute_clogp(molecule: Molecule) -> float:
    """Computes the cLogP of a molecule.

    :param molecule: A molecule, either a SMILES string or an RDKit molecule.
    :return: The cLogP of the molecule.
    """
    return MolLogP(molecule)


@register_property_function('mol_weight')
@smiles_to_mol_wrapper
def compute_molecular_weight(molecule: Molecule) -> float:
    """Computes the molecular weight of a molecule.

    :param molecule: A molecule, either a SMILES string or an RDKit molecule.
    :return: The molecular weight of the molecule.
    """
    return MolWt(molecule)


@register_property_function('qed')
@smiles_to_mol_wrapper
def compute_qed(molecule: Molecule) -> float:
    """Computes the QED (quantitative estimate of drug-likeness) score of a molecule.

    :param molecule: A molecule, either a SMILES string or an RDKit molecule.
    :return: The QED score of the molecule.
    """
    return qed(molecule)


@register_property_function('pains_plus')
@smiles_to_mol_wrapper
def compute_pains_plus(molecule: Molecule) -> str:
    """Checks a molecule for PAINS (pan-assay interference compounds) and other unwanted substructures.

    :param molecule: A molecule, either a SMILES string or an RDKit molecule.
    :return: The PAINS and other unwanted substructure(s) found in the molecule.
    """
    # Set up all unwanted substructures filters
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.ALL)
    catalog = FilterCatalog(params)

    # Check for unwanted substructures
    matches = " | ".join(
        f"{entry.GetProp('FilterSet')}: {entry.GetDescription()}" for entry in catalog.GetMatches(molecule)
    )

    return matches


@register_property_function('tpsa')
@smiles_to_mol_wrapper
def compute_tpsa(molecule: Molecule) -> float:
    """Computes the topological polar surface area (TPSA) of a molecule.

    :param molecule: A molecule, either a SMILES string or an RDKit molecule.
    :return: The TPSA of the molecule.
    """
    return CalcTPSA(molecule)


try:
    import os
    import sys

    from rdkit.Chem import RDConfig

    sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))

    import sascorer


    @register_property_function('sa_score')
    @smiles_to_mol_wrapper
    def compute_sa_score(molecule: Molecule) -> float:
        """Computes the synthetic accessibility (SA) score of a molecule.

        :param molecule: A molecule, either a SMILES string or an RDKit molecule.
        :return: The SA score of the molecule.
        """
        return sascorer.calculateScore(molecule)
except ModuleNotFoundError:
    @register_property_function('sa_score')
    @smiles_to_mol_wrapper
    def compute_sa_score(molecule: Molecule) -> float:
        """Computes the synthetic accessibility (SA) score of a molecule.

        :param molecule: A molecule, either a SMILES string or an RDKit molecule.
        :return: The SA score of the molecule.
        """
        raise ModuleNotFoundError('The SA Score module is not installed properly with RDKit. '
                                  'Please install it to use this property.')
