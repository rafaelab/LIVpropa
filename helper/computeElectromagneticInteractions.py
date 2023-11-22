import os
import sys
import numpy as np

# import from CRPropa3-data folder
sys.path.append('CRPropa3-data/')
import gitHelp as gh
import photonField

from kinematics import *
from interactionsElectromagnetic import *
from meanFreePath import computeInteractionRate
from constantsUnits import *
from warnings import filterwarnings

filterwarnings('ignore')


###############################################################################
###############################################################################
def process(interaction, field, kinematics = SpecialRelativity(), folder = '../data'):
	""" 
	Calculate the interaction rates for a given process on a given photon field .

	# Input
	. interaction: EM process
	. field   : photon field as defined in photonField.py
	. name    : name of the process which will be calculated. Necessary for the naming of the data folder
	. order   : order of the LIV (0 = symmetric)
	. energyQG: energy at which LIV sets in (defaults to Planck energy)
	. sign    : superluminal (+1), subluminal (-1)
	. kinematics: encapsulates LIV parameters
	"""
	particle = interaction.incidentParticle
	name = interaction.label + kinematics.label

	# output folder and kinematics-specific values
	if kinematics.label == 'SR':
		subfolder = ''

	elif kinematics.label == 'LIV':
		chi = kinematics.getChi(particle = particle)
		chiEl = kinematics.getChi(particle = 11)
		chiPh = kinematics.getChi(particle = 22)
		order = kinematics.getOrder()
		subfolder = 'chiEl_%+2.1e-chiPh_%+2.1e-order_%i' % (chiEl, chiPh, order)

	else:
		raise TypeError('Unknown type of kinematics.')

	if not folder.endswith('/'):
		folder += '/'
	folder = folder + name + '/' + subfolder
	if not os.path.exists(folder):
		os.makedirs(folder)

	# tabulated energies, limit to energies where the interaction is possible
	E = np.logspace(9, 23, 281) * eV
	Erange = (E[0], E[-1])

	Emin = interaction.minimumEnergyLab(field, Erange, kinematics = kinematics)
	E = E[E > Emin]
	Erange = (E[0], E[-1])
	
	# -------------------------------------------
	# calculate interaction rates
	# -------------------------------------------
	# tabulated values of s_kin = s - mc^2
	# Note: integration method (Romberg) requires 2^n + 1 log-spaced tabulation points
	sKin = np.logspace(4, 23, 2 ** 21 + 1) * eV ** 2
	xs = interaction.computeCrossSections(sKin)

	sThrIn = interaction.thresholdEnergy2Inner(E, Erange, kinematics = kinematics) 
	sThrOut = interaction.thresholdEnergy2Outer(E, kinematics = kinematics)

	rate = computeInteractionRate(sKin, xs, E, field, sThrIn, sThrOut, kinematics = kinematics, particle = particle)

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

	sKinMin = interaction.thresholdEnergy2(Erange, kinematics = kinematics)

	# tabulated values of s_kin = s - mc^2, limit to relevant range
	# Note: use higher resolution and then downsample
	sKin = np.logspace(4, 23, 380000 + 1) * eV ** 2
	sKin = sKin[sKin > sKinMin]

	xs = interaction.computeCrossSections(sKin)
	rate = computeInteractionRate(sKin, xs, E, field, sThrIn, sThrOut, cdf = True, kinematics = kinematics, particle = particle)

	# downsample
	sKin_save = np.logspace(4, 23, 190 + 1) * eV ** 2
	sKin_save = sKin_save[sKin_save > sKinMin]
	rate_save = np.array([np.interp(sKin_save, sKin, r) for r in rate])
	
	# save
	data = np.c_[np.log10(E / eV), rate_save]  # prepend log10(E/eV) as first column
	row0 = np.r_[0, np.log10(sKin_save / eV ** 2)][np.newaxis]
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

	del data, rate, sKin, sKin_save, rate_save




###########################################################################################
###########################################################################################
if __name__ == "__main__":

	pp = PairProduction()
	ics = InverseComptonScattering()

	interactions = [pp]
	fields = [
		photonField.CMB(),
		photonField.EBL_Gilmore12(),
		# photonField.URB_Protheroe96() # breaking down (why?)
	]
	orders = [1, 2]
	chis = np.logspace(-10., 10., 21, endpoint = True)
	chis = [1e-10, 1e-5, 1., 1e5]
	
	chisMix = [1e-3, 1e-2, 1e-1, 1., 1e1]

	# define the range within which chi values for photons and electrons will be allowed to differ
	chiRange = (1e-3, 1e3)

	for interaction in interactions:
		print('##########################################')
		print('=> interaction = ', interaction.name)

		for field in fields:
			print('  ==> photon field = ', field.name)

			## LIV 
			print('   ===> LIV')
			for order in orders:

				if interaction.label == 'PairProduction':
					for chiEl in chis:
						chiPh = chiEl
						chiDictP = {-11: chiEl, 11: chiEl, 22: chiPh}
						chiDictM = {-11: chiEl, 11: -1. * chiEl, 22: -1. * chiPh}
						livP = MonochromaticLIV(chi = chiDictP, order = order)
						livM = MonochromaticLIV(chi = chiDictM, order = order)	
						print('       ', livP)
						process(interaction, field, livP)
						print('       ', livM)
						process(interaction, field, livM)

						for chiPh in chisMix: # 
							chiRatio = np.abs(chiEl / chiPh)
							chiDictPP = {-11: chiEl, 11: chiEl, 22: chiPh}
							chiDictMM = {-11: chiEl, 11: -1. * chiEl, 22: -1. * chiPh}
							chiDictPM = {-11: chiEl, 11: chiEl, 22: -1. * chiPh}
							chiDictMP = {-11: chiEl, 11: -1. * chiEl, 22: chiPh}
							livPM = MonochromaticLIV(chi = chiDictPM, order = order)
							livMP = MonochromaticLIV(chi = chiDictMP, order = order)	
							livPP = MonochromaticLIV(chi = chiDictPP, order = order)
							livMM = MonochromaticLIV(chi = chiDictMM, order = order)	

							# run for different value
							if chiPh == chiEl and chiEl > chisMix[0] and chiEl < chisMix[-1]:
								print('       ', livPM)
								process(interaction, field, livPM)
								print('       ', livMP)
								process(interaction, field, livMP)
								print('       ', livPP)
								process(interaction, field, livPP)
								print('       ', livMM)
								process(interaction, field, livMM)


				elif interaction.label == 'InverseComptonScattering':
					for chiEl in chis:
						chiDictP = {-11: chiEl, 11: chiEl, 22: chiEl}
						chiDictM = {-11: chiEl, 11: -1. * chiEl, 22: -1. * chiEl}
						livP = MonochromaticLIV(chi = chiDictP, order = order)
						livM = MonochromaticLIV(chi = chiDictM, order = order)		
						print('       ', livP)
						process(interaction, field, livM)
						print('       ', livM)
						process(interaction, field, livP)
							

			## SR (for cross checks with CRPropa)
			sr = SpecialRelativity()
			print('  ===> SR')
			process(interaction, field, sr)