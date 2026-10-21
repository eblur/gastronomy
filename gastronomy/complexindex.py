
import numpy as np
import astropy.units as u
import kkcalc as kk
from .VY1995 import mineral_abs_xsect
from .minerals import Mineral

# K_EDGES = {
#     'C': 0.284 * u.keV,
#     'N': 0.401 * u.keV,
#     'O': 0.532 * u.keV,
#     'Mg': 1.305 * u.keV,
#     'Si': 1.839 * u.keV,
#     'Fe': 7.112 * u.keV
#     }

# L_EDGES = {'Fe': 0.725 * u.keV}

# UPDATING K, L, and M edges from VY1995:
K_EDGES = {'Be': np.float64(0.119), 'B': np.float64(0.194), 'C': np.float64(0.291), 'N': np.float64(0.405), 'O': np.float64(0.538), 
           'F': np.float64(0.694), 'Ne': np.float64(0.87), 'Na': np.float64(1.08), 'Mg': np.float64(1.31), 'Al': np.float64(1.57), 
           'Si': np.float64(1.85), 'P': np.float64(2.15), 'S': np.float64(2.48), 'Cl': np.float64(2.83), 'Ar': np.float64(3.2), 
           'K': np.float64(3.61), 'Ca': np.float64(4.04), 'Sc': np.float64(4.49), 'Ti': np.float64(4.97), 'V': np.float64(5.48), 
           'Cr': np.float64(6.0), 'Mn': np.float64(6.55), 'Fe': np.float64(7.12), 'Co': np.float64(7.72), 'Ni': np.float64(8.35), 
           'Cu': np.float64(8.99), 'Zn': np.float64(9.67)}


L1_EDGES = {'Al': np.float64(0.126), 'Si': np.float64(0.156), 'P': np.float64(0.194), 'S': np.float64(0.235), 'Cl': np.float64(0.278), 
            'Ar': np.float64(0.326), 'K': np.float64(0.384), 'Ca': np.float64(0.442), 'Sc': np.float64(0.503), 'Ti': np.float64(0.569), 
            'V': np.float64(0.638), 'Cr': np.float64(0.703), 'Mn': np.float64(0.782), 'Fe': np.float64(0.857), 'Co': np.float64(0.94), 
            'Ni': np.float64(1.02), 'Cu': np.float64(1.11), 'Zn': np.float64(1.2)}

L23_EDGES = {'Si': np.float64(0.106), 'P': np.float64(0.14), 'S': np.float64(0.17), 'Cl': np.float64(0.209), 'Ar': np.float64(0.249), 
             'K': np.float64(0.301), 'Ca': np.float64(0.352), 'Sc': np.float64(0.405), 'Ti': np.float64(0.464), 'V': np.float64(0.527), 
             'Cr': np.float64(0.585), 'Mn': np.float64(0.655), 'Fe': np.float64(0.724), 'Co': np.float64(0.8), 'Ni': np.float64(0.876), 
             'Cu': np.float64(0.947), 'Zn': np.float64(1.04)}

M1_EDGES = {'Fe': np.float64(0.104), 'Co': np.float64(0.115), 'Ni': np.float64(0.125), 'Cu': np.float64(0.129), 'Zn': np.float64(0.145)}

__all__ = ['calculate_Xray_complex_index_of_refraction', 
           'construct_egrid', 
           'mineral_to_formula', 
           'formula_to_mineral'
           ]

def calculate_Xray_complex_index_of_refraction(mineral, emin=0.1 * u.keV, emax=10.0 * u.keV, loggrid=True, **kwargs):
    """
    Calculate the complex index of refraction for a given mineral over a specified energy range.

    Inputs
    ------
    mineral : Mineral
        A gastronomy mineral object
    emin : astropy Quantity (default=0.1 * u.keV)
        The minimum energy of the grid
    emax : astropy Quantity (default=10.0 * u.keV)
        The maximum energy of the grid
    loggrid : bool (default=True)
        Whether to use a logarithmic grid spacing (True) or linear grid spacing (False)
    **kwargs : additional keyword arguments
        Additional arguments for `construct_energy_grid`

    Returns
    -------
    full_egrid : astropy Quantity array, astropy Quantity array, astropy Quantity array
        The energy grid output from kkcalc (includes much broader range than input)
    re_part : numpy array
         The refractive index (real part of the complex index of refraction).
    im_part : numpy array
         The absorption index (imaginary part of the complex index of refraction).
    """
    egrid = construct_egrid(mineral.elements, emin, emax, loggrid=loggrid, **kwargs)
    abs_xs = mineral_abs_xsect(egrid, mineral)
    full_egrid, delta, beta = _kkcalc_optical_constants(mineral, egrid, abs_xs)
    re_part = 1 - delta
    im_part = beta
    return full_egrid, re_part, im_part

def construct_egrid(element_list, emin, emax, loggrid=True,
                     ne=1000, edge_step=1.0 * u.eV, edge_margin_eV=30):
    """
    Construct an energy grid for X-ray calculations, with finer sampling around absorption edges.

    Inputs
    ------
    element_list : list of strings
        A list of element symbols (e.g. ['Fe', 'O']) for which to include absorption edges in the grid.
    emin : astropy Quantity
        The minimum energy of the grid (e.g. 0.1 * u.keV)
    emax : astropy Quantity
        The maximum energy of the grid (e.g. 10.0 * u.keV)
    loggrid : bool (default=True)
        Whether to use a logarithmic grid spacing (True) or linear grid spacing (False)
    ne : int (default=1000)
        The number of points in the initial broad grid (before adding edge regions)
    edge_step : astropy Quantity (default=1.0 * u.eV)
        The step size for the fine grid around absorption edges
    edge_margin_eV : float (default=30)
        The margin in eV around each edge to define the fine grid region
    
    Returns
    -------
    astropy Quantity array
        An energy grid that includes finer sampling around the absorption edges of the specified elements.
    """
    # Construct the initial broad-band grid
    emin_eV = emin.to('eV').value
    emax_eV = emax.to('eV').value
    if loggrid:
        egrid_eV = np.logspace(np.log10(emin_eV), np.log10(emax_eV), ne)
    else:
        egrid_eV = np.linspace(emin_eV, emax_eV, ne)
    
    # Gather the edge energies for the elements in the mineral
    edge_energies = []
    for element in element_list:
        if element in K_EDGES:
            edge_energies.append(K_EDGES[element] * u.keV)
        if element in L1_EDGES:
            edge_energies.append(L1_EDGES[element] * u.keV)
        if element in L23_EDGES:
            edge_energies.append(L23_EDGES[element] * u.keV)
        if element in M1_EDGES:
            edge_energies.append(M1_EDGES[element] * u.keV)
    
    # insert a fine grid around each edge
    final_grid = np.copy(egrid_eV)
    for edge in edge_energies:
        if (emin_eV < edge.to('eV').value) and (edge.to('eV').value < emax_eV):
            edge_eV = edge.to('eV').value
            step_eV = edge_step.to('eV').value
            fine_egrid = np.arange(edge_eV - edge_margin_eV, 
                                   edge_eV + edge_margin_eV, 
                                   step_eV)
            final_grid = np.hstack([final_grid[final_grid < fine_egrid[0]],
                                    fine_egrid,
                                    final_grid[final_grid > fine_egrid[-1]]])
    # return with the correct units
    return final_grid * u.eV

def mineral_to_formula(mineral):
    """
    Convert a gastronomy mineral to a chemical formula string

    Input
    ------
    mineral : Mineral
        A gastronomy mineral object
    
    Returns
    -------
    str : A chemical formula string representing the mineral

    Examples
    --------
    >>> mineral_to_formula(Mineral({'Fe': 2, 'O': 3}))
        "Fe2O3"
    >>> mineral_to_formula(Mineral({'Mg': 1, 'Si': 1, 'O': 3}))
        "MgSiO3"
    """
    result = ''
    for el in mineral.elements:
        result += f"{el}{mineral.number(el)}"
    return result

def formula_to_mineral(formula):
    """
    Convert a chemical formula string to a gastronomy mineral object

    Input
    ------
    formula : str
        A chemical formula string, e.g. "Fe2O3"

    Returns
    -------
    Mineral : A gastronomy mineral object with the corresponding elements and their counts
    
    Examples
    --------
    >>> formula_to_mineral("Fe2O3")
        Mineral with elements: {'Fe': 2, 'O': 3}
    >>> formula_to_mineral("MgSiO3")
        Mineral with elements: {'Mg': 1, 'Si': 1, 'O': 3}    
    """
    # This is a very naive implementation and should be improved
    import re
    pattern = r'([A-Z][a-z]?)(\d*)'
    matches = re.findall(pattern, formula)
    elements = {}
    for el, num in matches:
        if num == '':
            num = 1
        else:
            num = int(num)
        elements[el] = num
    return Mineral(elements)

def _kkcalc_optical_constants(mineral, egrid, abs_xs):
    ## THIS IS THE SAME AS kkcalc.kk_calculate_real
    """
    Calculate optical constants with kkcalc for a given compound.

    Inputs
    ------
    mineral : gastronomy Mineral object
        The mineral for which to calculate optical constants.
    egrid : astropy Quantity array
        The energy grid over which the absorption was computed.
    abs_xs : astropy Quantity array
        The absorption cross section data for the compound.
    
    Returns
    -------
    full_energy : astropy Quantity array
        The full energy grid over which optical constants are computed.
    delta : numpy array
        The refractive index decrement.
    beta : numpy array
        The absorption index.

    The returned model is n = 1 - delta + i*beta, where n is the complex index of refraction.
    """
    # Information needed to interpolate to an extended energy range
    stoichiometry = kk.data.ParseChemicalFormula(mineral_to_formula(mineral))
    # A relativistic correction, special to kkcalc library
    relativistic_correction = kk.calc_relativistic_correction(stoichiometry)
    # I assume this gets the template energy range and spectrum for the compound
    full_energy, imaginary_spectrum = kk.data.calculate_asf(stoichiometry)

    # kkcalc iput data (our absorption cross-section)
    xs_egrid_eV = egrid.to('eV').value
    xs_abs_m2 = abs_xs.to('m^2').value

    # Convert the cross-section data to atomic scattering factor (ASF)
    asf_lab = kk.data.convert_data(np.vstack((xs_egrid_eV, xs_abs_m2)).T, FromType='NEXAFS', ToType='asf')
    full_energy, imaginary_spectrum = kk.data.merge_spectra(asf_lab, full_energy, imaginary_spectrum, 
        merge_points=None, add_background=False, fix_distortions=False)
    # print(full_energy.shape, imaginary_spectrum.shape)

    # Compute the real part
    real_spectrum = kk.KK_PP(full_energy, full_energy, imaginary_spectrum, relativistic_correction)
    # I'm not sure what this does -- I guess it adds a value to the end so that the arrays are the same size?
    imaginary_spectrum_values = kk.data.coeffs_to_ASF(full_energy, np.vstack((imaginary_spectrum, imaginary_spectrum[-1])))

    # Convert the ASF to refractive index components delta and beta
    r_0 = 2.81794029957951365441605230194258e-15 * u.m # Classical electron radius in meters (from kkcalc.data)
    f1 = real_spectrum
    f2 = imaginary_spectrum_values
    Nq = (mineral.density / mineral.weight).to('cm^-3')
    const = r_0 * (full_energy * u.eV).to('cm', equivalencies=u.spectral())**2 / (2. * np.pi)
    delta = (const * Nq * f1).to('')
    beta = (const * Nq * f2).to('')

    return full_energy * u.eV, delta, beta
