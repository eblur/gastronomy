#
# minerals.py
# Created 2018.06.18 -- lia@astro.wisc.edu
#
# Some general atomic properties, useful for interpretting column densities
#-------------------------------------------------------------------------

import astropy.units as u
import astropy.constants as c
from .abundances import Z, amu

class Molecule(object):
    def __init__(self, composition, name=None):
        """
        Inputs
        ------
        name : string (optional)
            A string describing the molecule. Defaults to None

        composition : dictionary
            A dictionary describing the number of atoms of each element type
            in the molecule

            For example, CO_2 would be {'C':1, 'O':2}

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

        Example
        -------

        >>> CO2 = Molecule({'C':1, 'O':2}, name='CO_2')
        >>> CO2.elements
        dict_keys(['C', 'O'])
        """
        self.name = name
        self.composition = composition

    @property
    def elements(self):
        return self.composition.keys()

    @property
    def weight_amu(self):
        result = 0.0
        for ele in self.composition.keys():
            result += self.composition[ele] * amu[ele]
        return result

    @property
    def weight(self):
        return self.weight_amu * c.u

    def number(self, ele):
        """
        Returns the number of atoms, in the molecule, for the
        specified elemental species

        Inputs
        ------
        ele : string
            Periodic table abbreviation for the element of interest

        Example
        -------

        >>> CO2 = Molecule({'C':1, 'O':2})
        >>> CO2.number('O')
        2
        """
        if ele not in self.composition.keys():
            return 0.0
        else:
            return self.composition[ele]
