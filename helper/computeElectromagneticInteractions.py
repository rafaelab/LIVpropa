import argparse
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
def process(interaction, field, kinematics, folder = '../data'):
	""" 
	Calculate the interaction rates for a given process on a given photon field .

	# Input
	. interaction: EM process
	. field: photon field as defined in photonField.py
	. kinematics: kinematics of the process
	"""
	particle = interaction.incidentParticle
	name = interaction.label
	subfolder = kinematics.getIdentifier()

	if not folder.endswith('/'):
		folder += '/'
	folder = folder + name + '/' + subfolder
	if not os.path.exists(folder):
		os.makedirs(folder)

	# tabulated energies, limit to energies where the interaction is possible
	E = np.logspace(9, 23, 281) * eV
	Erange = (E[0], E[-1])

	Emin = interaction.minimumEnergyLab(field, Erange, kinematics)
	E = E[E > Emin]
	Erange = (E[0], E[-1])
	
	# -------------------------------------------
	# calculate interaction rates
	# -------------------------------------------
	# tabulated values of s_kin = s - mc^2
	# Note: integration method (Romberg) requires 2^n + 1 log-spaced tabulation points
	sKin = np.logspace(4, 23, 2 ** 21 + 1) * eV ** 2
	xs = interaction.computeCrossSections(sKin)

	sThrIn = interaction.thresholdEnergy2Inner(E, Erange, kinematics) 
	sThrOut = interaction.thresholdEnergy2Outer(E, kinematics)

	rate = computeInteractionRate(sKin, xs, E, field, sThrIn, sThrOut, kinematics, particle = particle)

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

	sKinMin = interaction.thresholdEnergy2(Erange, kinematics)

	# tabulated values of s_kin = s - mc^2, limit to relevant range
	# Note: use higher resolution and then downsample
	sKin = np.logspace(4, 23, 380000 + 1) * eV ** 2
	sKin = sKin[sKin > sKinMin]

	xs = interaction.computeCrossSections(sKin)
	rate = computeInteractionRate(sKin, xs, E, field, sThrIn, sThrOut, kinematics, cdf = True, particle = particle)

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

	interactions = [pp, ics]
	fields = [
		# photonField.CMB(),
		photonField.EBL_Gilmore12(),
		# photonField.URB_Protheroe96() # breaking down (why?)
	]
	orders = [1, 2]

	# define the range within which chi values for all particles are the same
	chiRange = np.array([1e-3, 1e-2, 1e-1, 1., 10.])
	chiRange = np.hstack([-chiRange[::-1], chiRange])

	sr = SpecialRelativity()
	kinSR = {-11: sr, 11: sr, 22: sr}

	listKinematics = []
	listKinematics.append(KinematicsMap(kinSR))

	#  assume all particles have the same chi value
	for order in orders:
		for chi in chiRange:
			liv = MonochromaticLIV(chi = chi, order = order)
			kinDict = {-11: liv, 11: liv, 22: liv}
			kinLIV = KinematicsMap(kinematicsDict = kinDict)
			listKinematics.append(kinLIV)

	#  assume chis photons and electrons/positrons are different
	chiRange = (-1e-1, 1e-1)
	for order in orders:
		for chiPh in chiRange:
			for chiEl in chiRange:
				livEl = MonochromaticLIV(chi = chiEl, order = order)
				livPh = MonochromaticLIV(chi = chiPh, order = order)
				livPo = livEl
				kinLIV = KinematicsMap(kinematicsDict = {22: livPh, 11: livEl, -11: livPo})
				listKinematics.append(kinLIV)


	for interaction in interactions:
		print('##########################################')
		print('=> interaction = ', interaction.name)

		for field in fields:
			print('. photon field = ', field.name)

			for kin in listKinematics:
				print('. ', kin.getIdentifier())
				process(interaction, field, kin)
							