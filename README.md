# gastronomy
Convenience functions and scripts for calculating gas properties from a few simple abundance tables.

## Install instructions

Download the source code using git clone or download and unpack the tar-ball from Github.

Navigate to the top "gastronomy" folder. Then run,

```bash
pip install -e .
```

Exit out of the directory before you try running it.

## Optional dependencies

If you plan to use the `complexindex` module (complex refractive index calculations), then the package `kkcalc` is required. Install it together with `gastronomy` from the source tree:

```bash
pip install .[complexindex]
```

### Note for Mac users

Zsh (the default terminal environment) treats `.[complexindex]` as a glob pattern. Use quotes instead, for example:

```bash
pip install '.[complexindex]'
```

<!-- From PyPI (published package):

```bash
pip install gastronomy[complexindex]
``` -->

## Quick start

To calculate an absorption cross-section of a compound from Verner & Yakovlev (1995) tables:

```python
import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
import gastronomy

egrid = np.linspace(0.1, 3.0, 1000000) * u.keV

compound = gastronomy.minerals.Mineral({'C':1, 'O':2}) # not a mineral, but this will work for any chemical formula

abs_xsect = gastronomy.VY1995.mineral_abs_xsect(egrid, compound)  # this will have units of area

plt.plot(egrid, abs_xsect.to('Mbarn'))
plt.loglog()
plt.xlabel(egrid.unit)
plt.ylabel('Mbarn per unit compound')
plt.show()
```

To compute the absorption cross-section for one ion from the Verner & Yakovlev (1995) tables:
```python
import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
import gastronomy

egrid = np.linspace(0.1, 3.0, 1000000) * u.keV

Z_Fe = 26 # atomic number for Iron

abs_xsect_FeI = gastronomy.VY1995.compute_xsect(egrid, Z_Fe, ion=0)  # this will have units of area
abs_xsect_FeIV = gastronomy.VY1995.compute_xsect(egrid, Z_Fe, ion=3) # He-like iron

plt.plot(egrid, abs_xsect_FeI.to('Mbarn'), label='Neutral Fe')
plt.plot(egrid, abs_xsect_FeIV.to('Mbarn'), label='Fe IV')
plt.loglog()
plt.legend()
plt.xlabel(egrid.unit)
plt.ylabel('Mbarn per unit compound')
plt.show()
```
