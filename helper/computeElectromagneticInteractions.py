import os
import sys
import numpy as np

# import from CRPropa3-data folder
sys.path.append('CRPropa3-data/')
import gitHelp as gh
import photonField

# import modified interaction rate
from interactionRateLIV import *

# import units and constants from CRPropa
from crpropa import eV, mass_electron, c_light, c_squared, sigma_thomson, alpha_finestructure



me2 = (mass_electron * c_squared) ** 2  # squared electron mass [J^2/c^4]
sigmaThomson = sigma_thomson  # Thomson cross section [m^2]
alpha = alpha_finestructure  # fine structure constant
MPl = 1.9561e9 # J



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

def getSmin(sigma):
	""" 
	Return minimum required s_kin = s - (mc^2)^2 for interaction 
	"""
	sigmas = {}
	sigmas[sigmaPP] = 4 * me2
	sigmas[sigmaDPP] = 16 * me2
	sigmas[sigmaTPP] = np.exp((218 / 27) / (28 / 9)) * me2 - me2
	sigmas[sigmaICS] = 1e-40 * me2 # arbitrary
	return sigmas[sigma]

def getEmin(sigma, field):
	""" 
	Return minimum required cosmic ray energy for interaction *sigma* with *field* 
	"""
	return getSmin(sigma) / 4. / field.getEmax()

def process(sigma, field, name, order = 1, sign = 1, energyQG = MPl, folder = '../data'):
	""" 
	Calculate the interaction rates for a given process on a given photon field .

	# Input
	. sigma   : crossection (function) of the EM-process
	. field   : photon field as defined in photonField.py
	. name    : name of the process which will be calculated. Necessary for the naming of the data folder
	. order   : order of the LIV (0 = symmetric)
	. energyQG: energy at which LIV sets in (defaults to Planck energy)
	. sign    : superluminal (+1), subluminal (-1)
	"""
	# adjust conventions
	xi = sign * MPl / energyQG
	nLIV = order + 2
	eta = 0.

	# output folder
	subfolder = 'Eqg_%2.1e-order_%i-%s' % (energyQG / eV, order, 'superluminal' if sign > 0 else 'subluminal')
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
	s_kin = np.logspace(4, 23, 2 ** 18 + 1) * eV ** 2
	xs = getTabulatedXS(sigma, s_kin)
	rate = calc_rate_s_liv(s_kin, xs, E, field, energyQG = energyQG, order = order, sign = sign)

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
	np.savetxt(fname, data, fmt=fmt, header=header)

	# -------------------------------------------
	# calculate cumulative differential interaction rates for sampling s values
	# -------------------------------------------
	# find minimum value of s_kin
	#skin1 = getSmin(sigma)  # s threshold for interaction
	#skin2 = np.maximum(np.maximum([xi * E[0]**2],[4*me2 + eta*E[0]**2]),[0])[0]
	#skin2 = 4 * field.getEmin() * E[0] + xi * E[0]**2  # minimum achievable s in collision with background photon (at any tabulated E)
	#skin_min = max(skin1, skin2)
	#skin_min = np.maximum(np.maximum([xi * E[0]**2],[4*me2 + eta*E[0]**2]),[0])[0]

	EminLIV = E[0] ** nLIV / MPl ** (nLIV - 2)
	EmaxLIV = E[-1] ** nLIV / MPl ** (nLIV - 2)
	skin_min = np.amin([xi * EminLIV, 4 * me2 + np.amin([eta * EminLIV, eta * EmaxLIV]), 4 * me2 + np.amin([xi * EminLIV, xi * EmaxLIV])])
	skin_min = np.amax([skin_min, 0.])

	# tabulated values of s_kin = s - mc^2, limit to relevant range
	# Note: use higher resolution and then downsample
	skin = np.logspace(4, 23, 380000 + 1) * eV ** 2
	skin = skin[skin > skin_min]

	xs = getTabulatedXS(sigma, skin)
	rate = calc_rate_s_liv(skin, xs, E, field, cdf = True, energyQG = energyQG, order = order, sign = sign)

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

	fields = [
		photonField.CMB(),
		photonField.EBL_Gilmore12(),
		photonField.URB_Protheroe96()
	]
	orders = [1, 2]
	energiesQG = [0.1 * MPl, MPl, 10 * MPl]


	for order in orders:
		for energyQG in energiesQG:
			for field in fields:
				print(order, energyQG / MPl, field.name)
				process(sigmaPP, field, 'PairProductionLIV', order = order, sign =  1, energyQG = energyQG)
				process(sigmaPP, field, 'PairProductionLIV', order = order, sign = -1, energyQG = energyQG)

				#        process(sigmaDPP, field, 'EMDoublePairProduction')
				#        process(sigmaTPP, field, 'EMTripletPairProduction')
				#        process(sigmaICS, field, 'EMInverseComptonScattering')
