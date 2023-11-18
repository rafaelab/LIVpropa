import os
import sys
import numpy as np
import warnings

from abc import ABC, abstractmethod

warnings.filterwarnings('ignore')

# import from CRPropa3-data folder
sys.path.append('CRPropa3-data/')
import gitHelp as gh
import photonField

# import LIV parameters definitions
from scenarioLIV import *

# import modified interaction rate
from interactionRateLIV import *

# units and constants
from constantsUnits import *


###############################################################################
###############################################################################
class ElectromagneticInteraction(ABC):
	"""
	Abstract base class for electromagnetic processes of the type:
	  X + gamma -> ...,
	wherein X is a photon or charged lepton.
	"""
	@property
	def name(self):
		return self._name
	
	@name.setter
	def name(self, n):
		self._name = n
	
	@property
	def label(self):
		return self._label
	
	@label.setter
	def label(self, l):
		self._label = l

	@property
	def sMin0(self):
		return self._sMin0
	
	@sMin0.setter
	def sMin0(self, s):
		self._sMin0 = s

	@property
	def incidentParticle(self):
		return self._incidentParticle
	
	@incidentParticle.setter
	def incidentParticle(self, p):
		self._incidentParticle = p

	@property
	def massInitial(self):
		return self._massInitial
	
	@massInitial.setter
	def massInitial(self, m):
		self._massInitial = m

	def thresholdEnergy2(self, Erange, kinematics = SpecialRelativity()):
		if kinematics.label == 'SR':
			return self.sMin0
		
		elif kinematics.label == 'LIV':
			particle = self.incidentParticle
			
			sMin = kinematics.computeDispersionCorrection(Erange[0], particle)
			sMax = kinematics.computeDispersionCorrection(Erange[1], particle)
			if sMin > sMax:
				sMin, sMax = sMax, sMin

			return np.maximum(sMin, self.sMin0)
		
		else:
			raise TypeError('Unknown type of kinematics.')
	
	@abstractmethod
	def thresholdEnergy2Outer(self, E, kinematics = SpecialRelativity()):
		"""
		Calculate threshold s for the outer interaction rate integral.
		"""
		pass

	def thresholdEnergy2Inner(self, E, Erange, kinematics = SpecialRelativity()):
		"""
		Calculate threshold s for the inner interaction rate integral.
		"""
		ds = kinematics.computeDispersionCorrection(E, self.incidentParticle)
		sMin = self.thresholdEnergy2(Erange, kinematics = kinematics)
		
		return np.maximum(sMin, ds)

	@abstractmethod
	def crossSection(self, s):
		pass

	@abstractmethod
	def computeCrossSections(self, sKin):
		""" 
		Get cross section for tabulated s_kin 
		"""
		pass

	def minimumEnergyLab(self, field, Erange, kinematics = SpecialRelativity()):
		""" 
		Return minimum required cosmic ray energy for interaction *sigma* with *field* 
		"""
		return self.thresholdEnergy2(Erange, kinematics = kinematics) / 4. / field.getEmax()

	# @abstractmethod
	# def minimumEnergyBackground(self, field, kinematics = SpecialRelativity()):
	# 	pass


###############################################################################
###############################################################################
class PairProduction(ElectromagneticInteraction):
	"""
	Breit-Wheeler pair production:
	  gamma + gamma -> e+ + e-
	"""
	def __init__(self):
		self.name = 'pair production'
		self.label = 'PairProduction'
		self.sMin0 = 4. * me2
		self.incidentParticle = 22
		self.massInitial = 0.

	def crossSection(self, s):
		""" 
		Pair production cross section (Breit-Wheeler), see Lee 1996 
		"""
		smin = 4 * me2
		if (s < smin):
			return 0.

		b = np.sqrt(1. - smin / s)
		return sigmaThomson * 3 / 16 * (1 - b * b) * ((3 - b ** 4) * (np.log1p(b) - np.log1p(-b)) - 2. * b * (2 - b * b))

	def computeCrossSections(self, sKin):
		""" 
		Get cross section for tabulated s_kin 
		"""
		return np.array([self.crossSection(s) for s in sKin])

	def thresholdEnergy2Outer(self, E, kinematics = SpecialRelativity()):		
		if kinematics.label == 'SR':
			return self.sMin0
		
		elif kinematics.label == 'LIV':
			chiPh = kinematics.getChi(particle = 22)
			if chiPh == 0.:
				return self.sMin0
			
			dsPh = kinematics.computeDispersionCorrection(E, 22)
			dsEl = kinematics.computeDispersionCorrection(E, 11)

			sThrPh = np.maximum(dsPh, self.sMin0)
			sThrEl = 4 * me2 + dsEl / 2 ** (kinematics.getNLIV() - 2)

			return np.maximum(sThrPh, sThrEl)

		else:
			raise TypeError('Unknown type of kinematics.')



###############################################################################
###############################################################################
class InverseComptonScattering(ElectromagneticInteraction):
	"""
	Inverse Compton Scattering:
	  e + gamma -> e + gamma
	"""
	def __init__(self):
		self.name = 'inverse Compton scattering'
		self.label = 'InverseComptonScattering'
		self.sMin0 = 1e-40 * me2
		self.incidentParticle = 11
		self.massInitial = mass_electron

	def crossSection(self, s):
		"""
		Inverse Compton scattering cross sections, see Lee 1996
		"""
		smin = me2
		if (s < smin):  # numerically unstable close to smin
			return 0.

		# note: formula unstable for (s - smin) / smin < 1E-5
		b = (s - smin) / (s + smin)
		A = 2. / b / (1 + b) * (2 + 2 * b - b * b - 2 * b * b * b)
		B = (2 - 3 * b * b - b * b * b) / b ** 2 * (np.log1p(b) - np.log1p(-b))

		return sigmaThomson * 3. / 8. * smin / s / b * (A - B)

	def computeCrossSections(self, sKin):
		""" 
		Get cross section for tabulated s_kin 
		"""
		return np.array([self.crossSection(sK + me2) for sK in (sKin)])

	def thresholdEnergy2Outer(self, E, kinematics = SpecialRelativity()):
		if kinematics.label == 'SR':
			return me2

		elif kinematics.label == 'LIV':
			chi = kinematics.getChi(particle = 11)
			if chi == 0.:
				return me2
			
			ds = kinematics.computeDispersionCorrection(E, 11)
			sThr = np.maximum(ds, self.sMin0)
			sThr = np.maximum(me2 + ds, sThr)
			
			return sThr
		
		else:
			raise TypeError('Unknown type of kinematics.')


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
		order = kinematics.getOrder()
		subfolder = 'chi_%+2.1e-order_%i' % (chi, order)

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
	sKin = np.logspace(4, 23, 2 ** 20 + 1) * eV ** 2
	xs = interaction.computeCrossSections(sKin)

	sThrIn = interaction.thresholdEnergy2Inner(E, Erange, kinematics = kinematics) 
	sThrOut = interaction.thresholdEnergy2Outer(E, kinematics = kinematics)

	rate = calc_rate_s_liv(sKin, xs, E, field, sThrIn, sThrOut, kinematics = kinematics, particle = particle)

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
	rate = calc_rate_s_liv(sKin, xs, E, field, sThrIn, sThrOut, cdf = True, kinematics = kinematics, particle = particle)

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
		photonField.CMB(),
		photonField.EBL_Gilmore12(),
		# photonField.URB_Protheroe96()
	]
	orders = [2]
	chis = np.logspace(-10., 10., 21, endpoint = True)
	chis = [1e-10, 1e-5, 1., 1e5]
	# interactions = [ics]
	# chis = [1e-10]

	for interaction in interactions:
		print('-------------------------------')
		print('=> interaction = ', interaction.name)
		for order in orders:
			for chi in chis:
				for field in fields:
					chiDictP = {-11: chi, 11: chi, 22: chi}
					chiDictM = {-11: chi, 11: -1. * chi, 22: -1. * chi}
					livP = MonochromaticLIV(chi = chiDictP, order = order)
					livM = MonochromaticLIV(chi = chiDictM, order = order)
					print(livP, ', photon field =', field.name)
					process(interaction, field, livM)
					print(livM, ', photon field =', field.name)
					process(interaction, field, livP)