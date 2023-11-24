#include "livpropa/LorentzSymmetry.h"



namespace livpropa {


LorentzSymmetry::LorentzSymmetry() {
}

LorentzSymmetry::LorentzSymmetry(int particle, unsigned int orderLIV, double parameterLIV) {
	setOrder(orderLIV);
	setCoefficient(parameterLIV);
	setParticleId(particle);
}

void LorentzSymmetry::setOrder(unsigned int orderLIV) {
	order = orderLIV;
}

void LorentzSymmetry::setCoefficient(double parameterLIV) {
	coefficient = parameterLIV;
}

void LorentzSymmetry::setParticleId(int id) {
	particleId = id;
}

unsigned int LorentzSymmetry::getOrder() const {
	return order;
}

double LorentzSymmetry::getCoefficient() const {
	return coefficient;
}

int LorentzSymmetry::getParticleId() const {
	return particleId;
}

double LorentzSymmetry::getEnergyScale() const {
	return abs(coefficient) * energy_planck;
}

bool LorentzSymmetry::isSuperluminal() const {
	return coefficient > 0 ? true : false; 
}

bool LorentzSymmetry::isSubluminal() const {
	return coefficient < 0 ? true : false; 
}

double LorentzSymmetry::getSymmetryBreakingShift(const double& p) const {
	return coefficient * pow((p * c_light) / energy_planck, order);
}

double LorentzSymmetry::computeEnergyFromMomentum(const double& p) const {
	double m = 0;
	switch (particleId) {
		case 22:
			m = 0;
			break;
		case 11:
			m = mass_electron;
			break;
		default:
			throw std::runtime_error("LorentzSymmetry: unknown particle id.");
	}
	
	double dE = getSymmetryBreakingShift(p);

	return sqrt(pow_integer<2>(m * c_squared) + pow_integer<2>(p * c_light) + dE);
}


} // namespace livpropa