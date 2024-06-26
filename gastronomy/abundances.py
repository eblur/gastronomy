#OA
# abundances.py -- a library of abundance dictionaries and atomic properties
# Created 2018.06.18 -- lia@astro.wisc.edu
#
# Some typical abundance tables
#-------------------------------------------------------------------------

import numpy as np

# Atomic numbers (Z)
Z = {'H':1, 'He':2, 'Li':3, 'Be':4, 'B':5, 'C':6, 'N':7, 'O':8, 'F':9, 'Ne':10, 
     'Na':11, 'Mg':12, 'Al':13, 'Si':14, 'P':15, 'S':16, 'Cl':17, 'Ar':18, 'K':19, 'Ca':20, 
     'Sc':21, 'Ti':22, 'V':23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30}

# Weights in atomic mass units (amu's)
amu = {'H':1.008,'He':4.0026,'Li':6.941,'Be':9.012182,'B':10.811, 'C':12.011,'N':14.007, 
       'O':15.999,'F':18.998403,'Ne':20.1797,'Na':22.989,'Mg':24.305,'Al':26.981, 
       'Si':28.085,'P':30.973,'S':32.06, 'Cl':35.45,'Ar':39.948,'K':39.0983, 
       'Ca':40.078,'Sc':44.95591,'Ti':47.867,'V':50.9415,'Cr':51.9961,'Mn':54.938, 
       'Fe':55.845,'Co':58.933,'Ni':58.6934,'Cu':63.546,'Zn':65.38}

# ISM abundances from Wilms, Allen, McCray (2000)
wilms_ISM = {'H':12.0, 'He':10.99, 'C':8.38, 'N':7.88, 'O':8.69, 'Ne':7.94, 
    'Na':6.16, 'Mg':7.40, 'Al':6.33, 'Si':7.27, 'P':5.42, 'S':7.09, 'Cl':5.12, 
    'Ar':6.41, 'Ca':6.20, 'Ti':4.81, 'Cr':5.51, 'Mn':5.34, 'Fe':7.43, 
    'Co':4.92, 'Ni':6.05} # 12 + log A_z

wilms_ISM_gas_ratio = {'H':1.0, 'He':1.0, 'C':0.5, 'N':1.0, 'O':0.6, 'Ne':1.0,
    'Na':0.25, 'Mg':0.2, 'Al':0.02, 'Si':0.1, 'P':0.6, 'S':0.6, 'Cl':0.5, 
    'Ar':1.0, 'Ca':0.003, 'Ti':0.002, 'Cr':0.03, 'Mn':0.07, 'Fe':0.3, 
    'Co':0.05, 'Ni':0.04}

# Put all iron in dust
gas_ratio_alt = {'H':1.0, 'He':1.0, 'C':0.5, 'N':1.0, 'O':0.6, 'Ne':1.0, 
    'Na':0.25, 'Mg':0.2, 'Al':0.02, 'Si':0.1, 'P':0.6, 'S':0.6, 'Cl':0.5, 
    'Ar':1.0, 'Ca':0.003, 'Ti':0.002, 'Cr':0.03, 'Mn':0.07, 'Fe':0.0, 
    'Co':0.05, 'Ni':0.04}

## ISM abundances from Lodders, Palme, Gaile (2009), Table 4
## https://ui.adsabs.harvard.edu/abs/2009LanB...4B..712L/abstract

lodders = {'H':12.0, 'He':10.925, 'Li':1.10, 'Be':1.38,  # using solar photosphere values
           'B':2.70, 'Ne':8.05, 'Ar':6.50, # using Solar Photosphere for these values
           'C':7.41, 'N':6.28, 'O':8.42, 'F':4.44, 'Na':6.29, 
           'Mg':7.55, 'Al':6.45, 'Si':7.53, 'P':5.45, 'S':7.17, 
           'Cl':5.25, 'K':5.10, 'Ca':6.31, 'Sc':3.07, 'Ti':4.93, 
           'V':3.98, 'Cr':5.66, 'Mn':5.50, 'Fe':7.47, 'Co':4.89, 
           'Ni':4.27, 'Cu':4.27, 'Zn':4.65} # I skipped elements heavier than this


# ISM abundances from Hensley & Draine (2021) as outlined in Table 1 of Draine & Hensley (2022)
def ppm_to_A(ppm):
    # Convert (X/H) values in ppm to np.log10(X/H) + 12.0
    return np.log10(ppm * 1.e-6) + 12.0

# He/H value is approximate
HD21_ISM = {'H':12.0, 'He':11.0,
            'C':ppm_to_A(324.),
            'O':ppm_to_A(682.),
            'Mg':ppm_to_A(52.9),
            'Al':ppm_to_A(3.48),
            'Si':ppm_to_A(44.6),
            'S':ppm_to_A(17.2),
            'Ca':ppm_to_A(3.25),
            'Fe':ppm_to_A(43.7),
            'Ni':ppm_to_A(2.09)}

HD21_ISM_gas_ratio = {'H':1.0, 'He':1.0, 'Li':1.0, 'N':1.0, 'Ne':1.0, 'Ar':1.0, # inserted by me
                      'C':198./324.,
                      'O':434./682.,
                      'Mg':7.1/52.9,
                      'Al':0.07/3.48,
                      'Si':6.6/44.6,
                      'S':9.6/17.2,
                      'Ca':0.07/3.25,
                      'Fe':0.9/43.7,
                      'Ni':0.04/2.09}

HD21_dust = {'C':ppm_to_A(126.),
             'O':ppm_to_A(248.),
             'Mg':ppm_to_A(46.),
             'Al':ppm_to_A(3.4),
             'Si':ppm_to_A(38.0),
             'S':ppm_to_A(7.6),
             'Ca':ppm_to_A(3.2),
             'Fe':ppm_to_A(43.),
             'Ni':ppm_to_A(2.1)}

def get_total_abund(elem, abund_table=wilms_ISM):
    """
    Returns
    -------
    float
        N_z (number of atoms of element with atomic number z)
        per H (hydrogen atom)

    Inputs
    ------
    elem : string
        Element name from the periodic table

    abund_table : dict
        A dictionary describing the abundance of Z per H (12 + log A_z).
        The default is to use the Wilms, Allens, & McCray (2000) abundances
        described in Table 2 of their paper.
    """
    assert type(elem) == str
    assert type(abund_table) == dict
    return np.power(10.0, abund_table[elem] - 12.0)

def get_dust_abund(elem, abund_table=wilms_ISM, gas_ratio=wilms_ISM_gas_ratio):
    """
    Returns
    -------
    float
        N_z (number of atoms of element with atomic number z)
        per H (hydrogen atom) expected to be locked in dust grains.

    Inputs
    ------
    elem : string
        Element name from the periodic table

    abund_table : dict
        A dictionary describing the abundance of Z per H (12 + log A_z).
        The default is to use the Wilms, Allens, & McCray (2000) abundances
        described in Table 2 of their paper.

    gas_ratio : dict
        A dictionary describing the fraction of element Z suspected to remain
        in the gas phase (where 1 - gas_ratio is the 'depletion', or fraction
        locked up in dust grains).
    """
    assert type(elem) == str
    assert type(abund_table) == dict
    assert type(gas_ratio) == dict
    return get_total_abund(elem, abund_table) * (1.0 - gas_ratio[elem])

def get_gas_abund(elem, abund_table=wilms_ISM, gas_ratio=wilms_ISM_gas_ratio):
    """
    Returns
    -------
    float
        N_z (number of atoms of element with atomic number z)
        per H (hydrogen atom) expected to be in the gas phase.

    Inputs
    ------
    elem : string
        Element name from the periodic table

    abund_table : dict
        A dictionary describing the abundance of Z per H (12 + log A_z).
        The default is to use the Wilms, Allens, & McCray (2000) abundances
        described in Table 2 of their paper.

    gas_ratio : dict
        A dictionary describing the fraction of element Z suspected to remain
        in the gas phase (where 1 - gas_ratio is the 'depletion', or fraction
        locked up in dust grains).
    """
    assert type(elem) == str
    assert type(abund_table) == dict
    assert type(gas_ratio) == dict
    return get_total_abund(elem, abund_table) * gas_ratio[elem]
