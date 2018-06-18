#
# minerals.py
# Created 2018.06.18 -- lia@astro.wisc.edu
#
# Some general atomic properties, useful for interpretting column densities
#-------------------------------------------------------------------------

import astropy.constants as c

# Weights in atomic mass units (amu's)
amu   = {'H':1.008,'He':4.0026,'C':12.011,'N':14.007,'O':15.999,'Ne':20.1797, \
         'Na':22.989,'Mg':24.305,'Al':26.981,'Si':28.085,'P':30.973,'S':32.06, \
         'Cl':35.45,'Ar':39.948,'Ca':40.078,'Ti':47.867,'Cr':51.9961,'Mn':54.938, \
         'Fe':55.845,'Co':58.933,'Ni':58.6934}

amu_g  = c.u.to('g').value # atomic mass unit (g)
mp     = c.m_p.to('g').value # proton mass (g)

class Mineral(object):
    def __init__(self, comp):
        self.composition = comp

    @property
    def weight_amu(self):
        result = 0.0
        for atom in self.composition.keys():
            result += self.composition[atom] * amu[atom]
        return result

    @property
    def weight_g(self):
        return self.weight_amu * amu_g
