import numpy as np
from abc import ABC, abstractmethod
from constantsUnits import *



###############################################################################
###############################################################################
class Kinematics(ABC):
	"""
	General abstract class holding information about a given type of kinematics.
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

	@abstractmethod
	def computeDispersionCorrectionMomentum(self, p, m):
		"""
		Compute correction to special-relativistic dispersion relation given the momentum
		"""
		pass

	@abstractmethod
	def computeDispersionCorrectionEnergy(self, E, m):
		"""
		Compute correction to special-relativistic dispersion relation given the energy.
		"""
		pass

	def computeDispersionCorrection(self, E, m):
		"""
		"""
		return self.computeDispersionCorrectionEnergy(E, m)

	@abstractmethod
	def computeEnergy2FromMomentum(self, p, m):
		pass

	def computeEnergyFromMomentum(self, p, m):
		return np.sqrt(self.computeEnergy2FromMomentum(p, m))

	# @abstractmethod
	# def computeMomentumFromEnergy(self, m, E):
	# 	"""
	# 	Not all equations can be solved analytically.
	# 	This method will not be implemented for the timebeing.
	# 	"""
	# 	pass

	@abstractmethod
	def getIdentifier(self):
		pass

	@abstractmethod
	def __str__(self):
		pass


###############################################################################
###############################################################################
class SpecialRelativity(Kinematics):
	"""
	The usual special-relativistic case.
	"""
	def __init__(self):
		self.name = 'special relativity'
		self.label = 'SR'

	def computeEnergy2FromMomentum(self, p, m):
		"""
		Returns m^2 + p^2
		"""
		return (m * c_squared) ** 2 + (p * c_light) ** 2
		
	def computeDispersionCorrectionEnergy(self, E, m):
		"""
		"""
		return 0.
	
	def computeDispersionCorrectionMomentum(self, p, m):
		"""
		"""
		return 0.
	
	def getIdentifier(self):
		"""
		"""
		return 'SR'

	def __str__(self):
		"""
		"""
		s = self.name
		return s


###############################################################################
###############################################################################
class MonochromaticLIV(Kinematics):
	"""
	A simple phenomenological implementation of LIV.
	This preserves energy-momentum conservations and modifies the dispersion relations as:
		E^2 = m^2 + p^2 + f(p).
	The fact that only a single value of n is taken into account defines the naming choice `monochromatic'.
	"""
	def __init__(self, chi, order = 0):
		"""
		"""
		self.chi = chi
		self.order = order
		self.name = 'monochromatic LIV'
		self.label = 'LIVmono'

	def getOrder(self):
		"""
		"""
		return self.order

	def isSuperluminal(self):
		"""
		"""
		if chi == 0:
			return False
		return chi > 0.
	
	def isSubluminal(self):
		"""
		"""
		if chi == 0:
			return False
		return chi < 0.

	def getChi(self):
		"""
		"""
		return self.chi

	def getNLIV(self):
		"""
		"""
		return self.order + 2
	
	def getSign(self):
		"""
		"""
		return np.sign(self.getChi())

	def computeDispersionCorrectionEnergy(self, E, m):
		"""
		Calculation of the term that modifies the SR dispersion relation in the ultrarelativistic limit
		Calculate: dE2 =  chi * (p / Eliv) ^ n
		"""		
		return self.chi * (E) ** (self.getNLIV()) / EPl ** (self.getNLIV() - 2)
	
	def computeDispersionCorrectionMomentum(self, p, m):
		"""
		Calculation of the term that modifies the SR dispersion relation.
		Calculate: dE2 =  chi * (p / Eliv) ^ n
		"""
		return self.chi * (p * c_light) ** (self.getNLIV()) / EPl ** (self.getNLIV() - 2)

	def computeEnergy2FromMomentum(self, p, m):
		"""
		Calculate: E^2 = m^2 + p^2 + chi * (p / Eliv) ^ n
		"""
		sr = SpecialRelativity()
		E2 = sr.computeEnergy2FromMomentum(p, m)
		dE2 = self.computeDispersionCorrection(p, m)
		return E2 + dE2

	def getIdentifier(self):
		"""
		"""
		s = self.label
		s += '%i' % self.order
		s += '_chi_%+2.1e' % self.chi
		return s

	def __str__(self):
		"""
		"""
		s = 'order = %i' % self.order
		s += ', chi_%i = %+2.1e' % self.chi
		return s


###############################################################################
###############################################################################
class KinematicsMap():
	"""
	This class is a simple dictionary that maps particles to their respective kinematics.
	"""
	def __init__(self, kinematicsDict = {}):
		self.kinematicsDict = kinematicsDict
	
	def add(self, particle, kinematics):
		"""
		Add an entry to the dictionary (a particle and its corresponding kinematics).
		"""
		self.kinematicsDict[particle] = kinematics

	def getParticles(self):
		"""
		Returns the list of particles.
		"""
		return list(self.kinematicsDict.keys())

	def getIdentifier(self):
		"""
		Returns a string with the format:
			Id_+{ID}-{KIN}_{KIN_DETAILS}
		Example: 
		. Id_+11-SR
		. Id_+22-LIVmono0_chi_+1.0e-5
		"""
		entries = []
		for particle in self.getParticles():
			kinematics = self.kinematicsDict[particle]
			s = 'Id_%+i_' % particle
			s += kinematics.getIdentifier()
			entries.append(s)

		return '-'.join(entries)
	
	def isSpecialRelativity(self):
		"""
		Returns True if all particles obey special relativity.
		"""
		for kinematics in self.kinematicsDict.values():
			if kinematics.label != 'SR':
				return False
		return True

	def isMonochromaticLIV(self):
		"""
		Returns True if at least one particle obeys LIV.
		"""
		for kinematics in self.kinematicsDict.values():
			if kinematics.label == 'LIVmono':
				return True
		return False

	def getKinematicsForParticle(self, particle):
		"""
		Returns the chi value for a given particle.
		"""
		return self.kinematicsDict[particle]

	

###############################################################################
###############################################################################
