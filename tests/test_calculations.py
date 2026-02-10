import pytest
import numpy as np
import astropy.units as u
from gastronomy import abundances
from gastronomy import VY1995

def test_H_abundance():
    # The H/H abundance should always be equal to 1
    assert abundances.get_total_abund('H') == 1
    assert abundances.get_gas_abund('H') == 1


# An arbitrary check on the ppm_to_A function
def test_ppm_to_A():
    test = abundances.ppm_to_A(324.)
    assert abundances.HD21_ISM['C'] == pytest.approx(test)

# An arbitrary check to show comparison of two arrays
def test_interpolate():
    # Compute a cross-section
    egrid = np.linspace(1.0, 2.0, 100) * u.keV
    Oxygen_xsect_0 = VY1995.compute_xsect(egrid, 8).value
    # Compute the same cross-section on a coarser grid
    egrid2 = np.linspace(1.0, 2.0, 20) * u.keV
    Oxygen_xsect2 = VY1995.compute_xsect(egrid2, 8).value
    # Compute the cross-section for coarser grid using interpolation instead
    Oxygen_xsect_test = np.interp(egrid2.value, egrid.value, Oxygen_xsect_0)
    # Check that they are about the same to within 0.1%
    assert Oxygen_xsect_test == pytest.approx(Oxygen_xsect2, rel=1.e-3)

