"""
Define constants and units.
Reexport the ones already defined in CRPropa.
"""

from crpropa import *

EPl = 1.9561e9 # J
me2 = 6.72236e-27 # J^2
me2 = (mass_electron * c_squared) ** 2  # squared electron mass [J^2/c^4]
sigmaThomson = sigma_thomson  # Thomson cross section [m^2]
alpha = alpha_finestructure  # fine structure constant
