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
	def getParticles(self):
		"""
		Returns list of all particles affected by this type of kinematics.
		All other particles obey, by default, special relativity
		"""
		pass

	@abstractmethod
	def computeDispersionCorrectionMomentum(self, p, particle):
		"""
		Compute correction to special-relativistic dispersion relation given the momentum
		"""
		pass

	@abstractmethod
	def computeDispersionCorrectionEnergy(self, E, particle):
		"""
		Compute correction to special-relativistic dispersion relation given the energy.
		"""
		pass

	def computeDispersionCorrection(self, E, particle):
		return self.computeDispersionCorrectionEnergy(E, particle)

	@abstractmethod
	def computeEnergy2FromMomentum(self, p, particle):
		pass

	def computeEnergyFromMomentum(self, p, particle):
		return np.sqrt(self.computeEnergy2FromMomentum(p, mass = mass, particle = particle))

	# @abstractmethod
	# def computeMomentumFromEnergy(self, m, E):
	# 	"""
	# 	Not all equations can be solved analytically.
	# 	This method will not be implemented for the timebeing.
	# 	"""
	# 	pass

	@abstractmethod
	def __str__(self):
		pass


###############################################################################
###############################################################################
class SpecialRelativity(Kinematics):
	"""
	The usual special-relativistic case.
	"""
	def __init__(self, particles = None):
		self.name = 'special relativity'
		self.label = 'SR'
		self.particles = particles

	def getParticles(self):
		return self.particles

	def computeEnergy2FromMomentum(self, p, particle):
		m = particleMassesDictionary[particle]
		return (m * c_squared) ** 2 + (p * c_light) ** 2
		
	def computeDispersionCorrectionEnergy(self, E, particle):
		return 0.
	
	def computeDispersionCorrectionMomentum(self, p, particle):
		return 0.
	
	def __str__(self):
		s = self.name
		return s


###############################################################################
###############################################################################
class MonochromaticLIV(Kinematics):
	"""
	A simple phenomenological implementation of LIV.
	This preserves energy-momentum conservations and modifies the dispersion relations as:
		E^2 = m^2 + p^2 + f^n.
	The fact that only a single value of n is taken into account defines the naming choice `monochromatic'.
	"""
	def __init__(self, chi = {11: 0., 22: 0.}, order = 0):
		self.chi = chi
		self.order = order
		self.name = 'monochromatic LIV'
		self.label = 'LIV'

	def getOrder(self):
		return self.order

	def getParticles(self):
		return self.chi.keys()
	
	def isParticleLorentzInvariant(self, particle):
		if particle in getParticles():
			return False
		return True

	def isSuperluminal(self, particle = None):
		if self.isParticleLorentzInvariant(particle):
			return False
		
		if particle is None:
			return np.any([chi >= 0. for chi in self.chi.values()])
		else:
			return self.chi[particle] >= 0. 
	
	def isSubluminal(self):
		if self.isParticleLorentzInvariant(particle):
			return False
		
		if particle is None:
			return np.any([chi < 0. for chi in self.chi.values()])
		else:
			return self.chi[particle] < 0. 

	def getChi(self, particle = None):
		if particle is None:
			return self.chi
		else:
			if particle in self.getParticles():
				return self.chi[particle]
			else:
				return 0.

	def getNLIV(self):
		return self.order + 2
	
	def getSign(self, particle = None):
		return np.sign(self.getChi(particle = particle))

	def computeDispersionCorrectionEnergy(self, E, particle):
		"""
		Calculation of the term that modifies the SR dispersion relation.
		"""		
		m = particleMassesDictionary[particle]
		return self.chi[particle] * (E) ** (self.getNLIV()) / EPl ** (self.getNLIV() - 2)
	
	def computeDispersionCorrectionMomentum(self, p, particle):
		"""
		Calculation of the term that modifies the SR dispersion relation.
		"""
		return self.chi[particle] * (p * c_light) ** (self.getNLIV()) / EPl ** (self.getNLIV() - 2)

	def computeEnergy2FromMomentum(self, p, particle):
		"""
		Note that particle and `m' should refer to the same particle
		"""
		sr = SpecialRelativity()
		E2 = sr.computeEnergy2FromMomentum(p, particle)
		dE2 = self.computeDispersionCorrection(p, particle)
		return E2 + dE2

	def __str__(self):
		s = 'order=%i' % self.order
		for particle in self.getParticles():
			s += ', chi_%i=%+2.1e' % (particle, self.chi[particle])
		return s


###############################################################################
###############################################################################
