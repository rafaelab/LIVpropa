import os
import sys
import numpy as np

from scipy.integrate import cumulative_trapezoid, romb, quad
from scipy.integrate import romb, quad

# import from CRPropa3-data folder
sys.path.append('CRPropa3-data/')
import gitHelp as gh
from interactionRate import *

# import units and constants from CRPropa
from crpropa import eV, Mpc

# units and constants
from constantsUnits import *

# import LIVParameters object
from scenarioLIV import *




def calc_rate_s_liv(s_kin, xs, E, field, sThr, parametersLIV = LIVParameters(order = 0, eta = 0.), z = 0, cdf = False):
	"""
	Calculate the interaction rate for given tabulated cross sections against an isotropic photon background.
	The tabulated cross sections need to be of length n = 2^i + 1 and the tabulation points log-linearly spaced.

	# Input
	. s_kin : tabulated (s - m**2) for cross sections [J^2]
	. xs    : tabulated cross sections [m^2]
	. E     : (array of) cosmic ray energies [J]
	. field : photon background, see photonField.py
	. z     : redshift
	. cdf   : calculate cumulative differential rate
	. liv   : an object of type `LIVParameter`
	. sThr  : thresholds for inner/outer integrals (this order)

	# Output
	. interaction rate 1/lambda(gamma) [1/Mpc] or
	. cumulative differential rate d(1/lambda)/d(s_kin) [1/Mpc/J^2]
	"""
	xi = parametersLIV.getXi()
	nLIV =  parametersLIV.getNLIV()
	eta = parametersLIV.getEta()
	sign = np.sign(eta)

	if cdf:
		# precalculate the field integral if it not exists and load it afterwards
		calculateDensityIntegral(field)
		file = 'temp/fieldDensity/' + field.name + '.txt'
		densityIntegral = np.loadtxt(file)

		# interpolate
		I = np.zeros((len(E), len(s_kin)))
		for j in range(len(E)):
			smin = sThr[1][j]
			npw = np.where(s_kin > smin)
			I[j, npw] = np.interp((s_kin[npw] - xi * E[j] ** nLIV / EPl ** (nLIV - 2)) / 4. / E[j], densityIntegral[:, 0], densityIntegral[:, 1])

		# calculate cdf
		y = np.array([xs * (s_kin - xi * E[j] ** nLIV / EPl ** (nLIV - 2) ) for j in range(len(E))]) * I

		if len(s_kin) == 0:
			return np.zeros(len(E))
		else:
			return cumulative_trapezoid(y = y, x = s_kin, initial = 0) / 8. / np.expand_dims(E, -1) ** 2 * Mpc    
		
	
	else:
		F = cumulative_trapezoid(x = s_kin, y = s_kin * xs, initial = 0)

		# new agorithm
		F1 = cumulative_trapezoid(x = s_kin, y = xs, initial = 0)
		Fmatrix = np.zeros((len(E), len(s_kin)))
		F1matrix = np.zeros((len(E), len(s_kin)))

		for ind, s in enumerate(sThr[0]):
			sbool = s_kin > s
			if np.any(sbool):
				npw = np.where(sbool)
				Fmatrix[ind, npw] = F[npw] - F[np.amin(npw)]
				F1matrix[ind, npw] = F1[npw] - F1[np.amin(npw)]

		Fmatrix = Fmatrix - xi * np.reshape(E, (-1, 1)) ** nLIV / EPl ** (nLIV - 2) * F1matrix
		sE = np.subtract.outer(-xi * E ** nLIV / EPl ** (nLIV - 2), -s_kin)

		n = field.getDensity(sE / np.reshape(4. * E, (-1, 1)), z)
		idx = np.where(np.isnan(n))
		n[idx] = 0.

		y = n * Fmatrix / sE ** 2 * s_kin

		# sthr = np.maximum(4 * me2 + eta / 2 ** (nLIV - 2) * E ** nLIV / EPl ** (nLIV - 2), np.zeros(len(E)))
		for ind, s in enumerate(sThr[1]):
			sbool = s_kin < s
			if np.any(sbool):
				npw = np.where(sbool)
				y[ind, npw] = 0.

		ds = mean_log_spacing(s_kin)
		return romb(y, dx = ds) / 2. / E * Mpc

