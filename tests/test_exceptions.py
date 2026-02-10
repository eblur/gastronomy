import pytest
import numpy as np
from gastronomy import VY1995

def test_VY1995_quantities():
    egrid = np.linspace(1.0, 2.0, 10) # not a Quantity
    with pytest.raises(TypeError, match='Expected an astropy Quantity'):
        VY1995.compute_xsect(egrid, 10)
    
