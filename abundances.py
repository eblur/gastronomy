#
# abundances.py
# Created 2018.06.18 -- lia@astro.wisc.edu
#
# Some typical abundance tables
#-------------------------------------------------------------------------

# ISM abundances from Wilms, Allen, McCray (2000)
wilms_ISM = {'H':12.0, 'He':10.99, 'C':8.38, 'N':7.88, 'O':8.69, 'Ne':7.94, \
    'Na':6.16, 'Mg':7.40, 'Al':6.33, 'Si':7.27, 'P':5.42, 'S':7.09, 'Cl':5.12, \
    'Ar':6.41, 'Ca':6.20, 'Ti':4.81, 'Cr':5.51, 'Mn':5.34, 'Fe':7.43, \
    'Co':4.92, 'Ni':6.05} # 12 + log A_z

wilms_ISM_gas_ratio = {'H':1.0, 'He':1.0, 'C':0.5, 'N':1.0, 'O':0.6, 'Ne':1.0,
    'Na':0.25, 'Mg':0.2, 'Al':0.02, 'Si':0.1, 'P':0.6, 'S':0.6, 'Cl':0.5, \
    'Ar':1.0, 'Ca':0.003, 'Ti':0.002, 'Cr':0.03, 'Mn':0.07, 'Fe':0.3, \
    'Co':0.05, 'Ni':0.04}

# Put all iron in dust
gas_ratio_alt = {'H':1.0, 'He':1.0, 'C':0.5, 'N':1.0, 'O':0.6, 'Ne':1.0, \
    'Na':0.25, 'Mg':0.2, 'Al':0.02, 'Si':0.1, 'P':0.6, 'S':0.6, 'Cl':0.5, \
    'Ar':1.0, 'Ca':0.003, 'Ti':0.002, 'Cr':0.03, 'Mn':0.07, 'Fe':0.0, \
    'Co':0.05, 'Ni':0.04}

def get_ISM_abund(elem, abund_table=wilms_ISM):
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
    return get_ISM_abund(elem, abund_table) * (1.0 - gas_ratio[elem])

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
    return get_ISM_abund(elem, abund_table) * gas_ratio[elem]
