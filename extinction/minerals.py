#
# minerals.py
# Created 2018.06.18 -- lia@astro.wisc.edu
#
# Some general atomic properties, useful for interpretting column densities
#-------------------------------------------------------------------------

from molecules import Molecule
from abundances import amu
import astropy.constants as c

class Mineral(Molecule):
    """
    Subclass of molecules.Molecule because it requires some of the same
    functionality.

    Inputs
    ------
        name : string (optional)
            A string describing the molecule. Defaults to None

        composition : dictionary
            A dictionary describing the number of atoms of each element type
            in the molecule

    Attributes
    ----------
    name : string (or None)

    composition : dictionary

    weight_amu : float
        The weight of the molecule in atomic mass units (amu)

    weight : astropy.units.Quantity
        The mass of the molecule

    elements : list
        Returns the elements that make up the keys for self.composition

    unit_mass : float
        Returns the mass of the mineral unit cell
    """
    def __init__(self, composition, name=None):
        Molecule.__init__(self, composition, name=name)

    @property
    def unit_mass(self):
        return sum([self.number(ele) * amu[ele] * c.u.to('g')
            for ele in self.elements])
