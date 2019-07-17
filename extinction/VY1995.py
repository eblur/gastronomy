#
# VY1995.py
# Created 2018.06.18 -- lia@astro.wisc.edu
#
# Functions for loading and calling the Verner & Yaklovev 1995 photoionization
#   cross-sections (http://adsabs.harvard.edu/abs/1995A%26AS..109..125V)
#
#-----------------------------------------------------------------------------

import numpy as np

import astropy.units as u
from astropy.table import Table

VY_TABLE_FILE = 'VY1995_vizier.txt'

def __write_vy1995_table(filename=VY_TABLE_FILE):
    from astroquery.vizier import Vizier
    TABLE_NAME = 'J/A+AS/109/125/table1'
    Vizier.ROW_LIMIT = -1
    vy1995 = Vizier.get_catalogs(TABLE_NAME)[TABLE_NAME]
    vy1995.write(filename, format='ascii', overwrite=True)
    return

vy1995 = Table.read(VY_TABLE_FILE, format='ascii')

def get_element(z, ion=None, shell=None):
    """
    Returns
    -------
    Boolean array that filters the Verner & Yaklovev 1995 table of analytic formulae
    for atomic photoelectric absorption cross-sections.

    Inputs
    ------
    z : int
        Atomic number

    ion : int
        Ion state (0 = neutral, 1 = singly ionized, etc..)
        The default (None) will return a filter for any ion state of
        element with atomic number `z`.

    shell : int
        The Principle quantum number of the shell.
        The default (None) will return a filter for any shell with
        atomic number `z`.
    """
    result = (vy1995['Z'] == z)

    if ion is not None:
        assert (ion >= 0) and (ion < z)
        n_electrons = z - ion
        result = result & (vy1995['N'] == n_electrons)

    if shell is not None:
        result = result & (vy1995['n'] == shell)

    return result

def _compute_xsect_from_row(egrid, row):
    """
    Returns
    -------
    astropy.units.Quantity (array, units of 'Mbarn')
        The photoelectric absorption cross-section computed from the parameters
        given by a single row in the Verner & Yaklovev 1995 VizieR table

    Inputs
    ------
    egrid : astropy.units.Quantity (array)
        Energy values for which to compute the cross-sections

    row : astropy.table.Table row (or an iterable list of length 10)
        Parameters of the Verner & Yaklovev 1995 fit, in the following order:

        Z, N, n, l, E_th (eV), E_0 (eV), sigma_0 (Mbarn), y_a, P, y_w
    """
    Z, N, n, l, E_th, E_0, sigma_0, y_a, P, y_w = row
    #print(Z, N, n, l, E_th, E_0, sigma_0, y_a, P, y_w)

    # VY1995 formula
    y = egrid.to('eV') / (E_0 * u.eV)
    Q = 5.5 + l - 0.5 * P
    F = ((y-1.)**2 + y_w**2) * np.power(y, -Q) * np.power((1. + np.sqrt(y/y_a)), -P)
    formula = sigma_0 * u.Mbarn * F

    # everything below the threshold energy is zero
    efilt   = egrid.to('eV').value > E_th
    result  = np.zeros(len(egrid)) * u.Mbarn
    result[efilt] = formula[efilt]
    return result

def compute_xsect(egrid, z, ion=None, shell=None):
    """
    Returns
    -------
    astropy.units.Quantity (array, units of 'Mbarn')
        Returns the total photoelectric absorption cross-section computed for
        the user specified element, ion(s), and shell(s).

    Inputs
    ------
    egrid : astropy.units.Quantity (array)
        Energy values for which to compute the cross-sections

    z : int
        Atomic number

    ion : int
        Ion state (0 = neutral, 1 = singly ionized, etc..)
        The default (None) will return a filter for any ion state of
        element with atomic number `z`.

    shell : int
        The Principle quantum number of the shell.
        The default (None) will return a filter for any shell with
        atomic number `z`.
    """
    result = np.zeros(len(egrid)) * u.Mbarn
    filt   = get_element(z, ion=ion, shell=shell)
    rows   = vy1995[filt]
    for r in rows:
        result += _compute_xsect_from_row(egrid, r)
    return result
