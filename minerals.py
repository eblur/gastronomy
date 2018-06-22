#
# minerals.py
# Created 2018.06.18 -- lia@astro.wisc.edu
#
# Some general atomic properties, useful for interpretting column densities
#-------------------------------------------------------------------------

from molecules import Molecule
import astropy.constants as c

class Mineral(Molecule):
    """
    Subclass of molecules.Molecule because it requires some of the same
    functionality.
    """
    def __init__(self, composition, name=None):
        Molecule.__init__(self, composition, name=name)
