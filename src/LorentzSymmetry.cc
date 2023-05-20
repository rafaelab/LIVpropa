#include "livpropa/LorentzSymmetry.h"



namespace livpropa {


LorentzSymmetry::LorentzSymmetry(double energyLIV, unsigned int orderLIV, double parameterLIV) {
	setEnergyScale(energyLIV);
	setOrder(orderLIV);
	setParameter(parameterLIV);
}

void LorentzSymmetry::setOrder(unsigned int orderLIV) {
	order = orderLIV;
}

void LorentzSymmetry::setParameter(double parameterLIV) {
	parameter = parameterLIV;
}

void LorentzSymmetry::setEnergyScale(double energyQG) {
	energyScale = energyQG;
}

unsigned int LorentzSymmetry::getOrder() const {
	return order;
}

double LorentzSymmetry::getParameter() const {
	return parameter;
}

double LorentzSymmetry::getEnergyScale() const {
	return energyScale;
}

bool LorentzSymmetry::isSuperluminal() const {
	return parameter > 0 ? true : false; 
}

bool LorentzSymmetry::isSubluminal() const {
	return parameter < 0 ? true : false; 
}


} // namespace livpropa