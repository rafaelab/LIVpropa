import numpy as np
from abc import ABC, abstractmethod
from kinematics import *
from constantsUnits import *
from warnings import filterwarnings

filterwarnings('ignore')

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