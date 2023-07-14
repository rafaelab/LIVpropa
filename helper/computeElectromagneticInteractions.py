import os
import sys
import numpy as np
import warnings

warnings.filterwarnings('ignore')

# import from CRPropa3-data folder
sys.path.append('CRPropa3-data/')
import gitHelp as gh
import photonField

# import LIV parameters definitions
from scenarioLIV import LIVParameters

# import modified interaction rate
from interactionRateLIV import *

# units and constants
from constantsUnits import *


def sigmaPP(s):
	""" 
	Pair production cross section (Breit-Wheeler), see Lee 1996 
	"""
	smin = 4 * me2
	if (s < smin):
		return 0.

	b = np.sqrt(1. - smin / s)
	return sigmaThomson * 3 / 16 * (1 - b * b) * ((3 - b ** 4) * (np.log1p(b) - np.log1p(-b)) - 2. * b * (2 - b * b))

def sigmaDPP(s):
	""" 
	Double-pair production cross section, see R.W. Brown eq. (4.5) with k^2 = q^2 = 0 
	"""
	smin = 16. * me2
	if (s < smin):
		return 0.

	return 6.45e-34 * (1 - smin / s) ** 6

def sigmaICS(s):
	"""
	Inverse Compton scattering cross sections, see Lee 1996
	"""
	smin = me2
	if (s < smin):  # numerically unstable close to smin
		return 0

	# note: formula unstable for (s - smin) / smin < 1E-5
	b = (s - smin) / (s + smin)
	A = 2. / b / (1 + b) * (2 + 2 * b - b * b - 2 * b * b * b)
	B = (2 - 3 * b * b - b * b * b) / b ** 2 * (np.log1p(b) - np.log1p(-b))
	return sigmaThomson * 3. / 8. * smin / s / b * (A - B)

def sigmaTPP(s):
	""" 
	Triplet-pair production cross section, see Lee 1996 
	"""
	beta = 28. / 9. * np.log(s / me2) - 218. / 27.
	if beta < 0:
		return 0.
	
	return sigmaThomson * 3. / 8. / np.pi * alpha * beta

def getTabulatedXS(sigma, skin):
	""" 
	Get cross section for tabulated s_kin 
	"""
	if sigma in (sigmaPP, sigmaDPP):  # photon interactions
		return np.array([sigma(s) for s in skin])
	if sigma in (sigmaTPP, sigmaICS):  # electron interactions
		return np.array([sigma(s) for s in skin + me2])
	return False

def getSmin(sigma, parametersLIV = LIVParameters(order = 0, eta = 0.), Erange = (1e9 * eV, 1e23 * eV)):
	""" 
	Return minimum required s_kin = s - (mc^2)^2 for interaction.
	LIV parameters are defined as keyword arguments.
	Note that LIV is not implemented for DPP and TPP.

	The last parameter (`Erange`) is provided to prevent NaNs.
	"""
	xi = parametersLIV.getXi()
	nLIV =  parametersLIV.getNLIV()
	eta = parametersLIV.getEta()
	sign = np.sign(eta)

	smin = {}
	smin[sigmaDPP] = 16 * me2
	smin[sigmaTPP] = np.exp((218 / 27) / (28 / 9)) * me2 - me2
	smin[sigmaPP] = 4 * me2
	smin[sigmaICS] = 1e-40 * me2 # arbitrary
	if order > 0:
		# adjust conventions
		xiGamma = sign * EPl / energyQG
		xiElectron = sign * EPl / energyQG
		nLIV = order + 2
		eta = 0.

		EminLIV = Erange[0] ** nLIV / EPl ** (nLIV - 2)
		EmaxLIV = Erange[1] ** nLIV / EPl ** (nLIV - 2)

		smin[sigmaPP] = np.amax([np.amin([xi * EminLIV, xi * EmaxLIV]), smin[sigmaPP]])
		smin[sigmaICS] = np.amax([np.amin([xi * EminLIV, xi * EmaxLIV]), smin[sigmaICS]])
		
	return smin[sigma]

def getSthrOuter(sigma, E, parametersLIV):
	"""
	Calculate s_thr for the outer integral.
	This is defined only for PP and ICS.
	"""
	xi = parametersLIV.getXi()
	nLIV =  parametersLIV.getNLIV()
	smin0 = getSmin(sigma, parametersLIV = LIVParameters(order = 0, eta = 0))

	if sigma == sigmaPP:
		sthr = np.amax([xi * E ** nLIV / EPl ** (nLIV - 2), smin0])
		sthr = np.amax([smin0 + xi / 2 ** (nLIV - 2) * E ** nLIV / EPl ** (nLIV - 2), sthr])
		return sthr
	elif sigma == sigmaICS:
		sthr = np.amax([xi * E ** nLIV / EPl ** (nLIV - 2), smin0])
		sthr = np.amax([me2 + xi * E ** nLIV / EPl ** (nLIV - 2), sthr])
	else:
		return smin0

def getSthrInner(sigma, E, parametersLIV):
	"""
	Calculate s_thr for the inner integral.
	"""
	xi = parametersLIV.getXi()
	nLIV =  parametersLIV.getNLIV()
	smin0 = getSmin(sigma, parametersLIV = LIVParameters(order = 0, eta = 0))

	if sigma == sigmaPP or sigma == sigmaICS:
		sthr = xi * E ** nLIV / EPl ** (nLIV - 2)
		return np.amax([sthr, smin0])
	else:
		return smin0

def getEmin(sigma, field):
	""" 
	Return minimum required cosmic ray energy for interaction *sigma* with *field* 
	"""
	return getSmin(sigma) / 4. / field.getEmax()

def process(sigma, field, name, parametersLIV = LIVParameters(order = 0, eta = 0.), folder = '../data'):
	""" 
	Calculate the interaction rates for a given process on a given photon field .

	# Input
	. sigma   : crossection (function) of the EM-process
	. field   : photon field as defined in photonField.py
	. name    : name of the process which will be calculated. Necessary for the naming of the data folder
	. order   : order of the LIV (0 = symmetric)
	. energyQG: energy at which LIV sets in (defaults to Planck energy)
	. sign    : superluminal (+1), subluminal (-1)
	. parametersLIV: encapsulates LIV parameters
	"""
	eta = parametersLIV.getEta()
	sign = np.sign(eta)

	# output folder
	subfolder = 'Eqg_%2.1eeV-order_%i-%s' % (energyQG / eV, order, 'superluminal' if sign > 0 else 'subluminal')
	if not folder.endswith('/'):
		folder += '/'
	folder = folder + name + '/' + subfolder
	if not os.path.exists(folder):
		os.makedirs(folder)

	# tabulated energies, limit to energies where the interaction is possible
	Emin = getEmin(sigma, field)
	E = np.logspace(9, 23, 281) * eV
	E = E[E > Emin]
	
	# -------------------------------------------
	# calculate interaction rates
	# -------------------------------------------
	# tabulated values of s_kin = s - mc^2
	# Note: integration method (Romberg) requires 2^n + 1 log-spaced tabulation points
	s_kin = np.logspace(4, 23, 2 ** 20 + 1) * eV ** 2
	xs = getTabulatedXS(sigma, s_kin)

	sThrIn = np.array([getSthrInner(sigma, e, parametersLIV) for e in E])
	sThrOut = np.array([getSthrOuter(sigma, e, parametersLIV) for e in E])

	rate = calc_rate_s_liv(s_kin, xs, E, field, (sThrIn, sThrOut), parametersLIV = parametersLIV)

	# print('---------------- ', s_kin.shape, xs.shape, rate.shape)
	# for i in range(len(s_kin)):
	# 	if xs[i] <= 0. or xs[i] > 1e10:
	# 		print(s_kin[i], xs[i])
	# for j in range(len(rate)):
	# 	print(E[j] / eV, rate[j])

	# save
	fname = folder + '/rate_%s.txt' % field.name
	data = np.c_[np.log10(E / eV), rate]
	fmt = '%.2f\t%8.7e'
	try:
		git_hash = gh.get_git_revision_hash()
		headerStr = '%s interaction rates\nphoton field: %s\n' % (name, field.info)
		headerStr += 'Produced with crpropa-data version: ' + git_hash + '\n'
		headerStr += 'log10(E/eV), 1/lambda [1/Mpc]'
		header = (headerStr)
	except:
		headerStr = '%s interaction rates\n' % (name)
		headerStr += 'photon field: %s\n' % (field.info)
		headerStr += 'log10(E/eV), 1/lambda [1/Mpc]'
		header = (headerStr)
	np.savetxt(fname, data, fmt = fmt, header = header)

	# -------------------------------------------
	# calculate cumulative differential interaction rates for sampling s values
	# -------------------------------------------

	skin_min = getSmin(sigma, parametersLIV, Erange = (E[0], E[-1]))

	# tabulated values of s_kin = s - mc^2, limit to relevant range
	# Note: use higher resolution and then downsample
	skin = np.logspace(4, 23, 380000 + 1) * eV ** 2
	skin = skin[skin > skin_min]

	xs = getTabulatedXS(sigma, skin)
	rate = calc_rate_s_liv(skin, xs, E, field, (sThrIn, sThrOut), cdf = True, parametersLIV = parametersLIV)

	# downsample
	skin_save = np.logspace(4, 23, 190 + 1) * eV ** 2
	skin_save = skin_save[skin_save > skin_min]

	rate_save = np.array([np.interp(skin_save, skin, r) for r in rate])
	
	# save
	data = np.c_[np.log10(E / eV), rate_save]  # prepend log10(E/eV) as first column
	row0 = np.r_[0, np.log10(skin_save / eV ** 2)][np.newaxis]
	data = np.r_[row0, data]  # prepend log10(s_kin/eV^2) as first row

	fname = folder + '/cdf_%s.txt' % field.name
	fmt = '%.2f' + '\t%6.5e' * np.shape(rate_save)[1]
	try:
		git_hash = gh.get_git_revision_hash()
		headerStr = '%s cumulative differential rate\n' % name
		headerStr += 'photon field: %s\n' % field.info
		headerStr += 'Produced with crpropa-data version: ' + git_hash + '\n'
		headerStr += 'log10(E/eV), d(1/lambda)/ds_kin [1/Mpc/eV^2] for log10(s_kin/eV^2) as given in first row'
		header = (headerStr)
	except:
		headerStr = '%s cumulative differential rate' % name
		headerStr += 'photon field: %s\n' % field.info
		headerStr +='log10(E/eV), d(1/lambda)/ds_kin [1/Mpc/eV^2] for log10(s_kin/eV^2) as given in first row' 
		header = (header)
	np.savetxt(fname, data, fmt = fmt, header = header)

	del data, rate, skin, skin_save, rate_save


###########################################################################################
###########################################################################################
if __name__ == "__main__":

	interactions = ['PairProductionLIV', 'InverseComptonScatteringLIV']
	fields = [
		photonField.CMB(),
		# photonField.EBL_Gilmore12(),
		# photonField.URB_Protheroe96()
	]
	orders = [1, 2]
	energiesQG = np.logspace(-9., 1., 11, endpoint = True) * EPl
	orders = [1]

	# interactions = ['PairProductionLIV']
	interactions = ['InverseComptonScatteringLIV']

	for interaction in interactions:
		print('-------------------------------')
		print('=> interaction = %s' % interaction)
		for order in orders:
			for energyQG in energiesQG:
				for field in fields:

					livM = LIVParameters(eta = -1., energyQG = energyQG, order = order)
					print(livM, 'photon field =', field.name)
					process(sigmaPP, field, interaction, livM)

					livP = LIVParameters(eta = 1., energyQG = energyQG, order = order)
					print(livP, 'photon field =', field.name)
					process(sigmaPP, field, interaction, livP)