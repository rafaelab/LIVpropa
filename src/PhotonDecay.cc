#include "livpropa/PhotonDecay.h"


namespace livpropa {



PhotonDecay::PhotonDecay(Kinematics kinematics, bool haveElectrons, double thinning, double limit) {
	setHaveElectrons(haveElectrons);
	setLimit(limit);
	setThinning(thinning);
	setKinematics(kinematics);
	setInteractionTag("PD");
}

void PhotonDecay::setKinematics(Kinematics kin) {
	kinematics = kin;
}

void PhotonDecay::setHaveElectrons(bool electrons) {
	haveElectrons = electrons;
}

void PhotonDecay::setLimit(double l) {
	limit = l;
}

void PhotonDecay::setThinning(double t) {
	thinning = t;
}

void PhotonDecay::setInteractionTag(string tag) {
	interactionTag = tag;
}

bool PhotonDecay::isImplemented() const {
	std::string typeEl = kinematics[ 11]->getNameTag();
	std::string typePo = kinematics[-11]->getNameTag();
	std::string typePh = kinematics[ 22]->getNameTag();

	if (typeEl != typePo)
		return false;

	if (typePh != "SR" and (typePh != typeEl or typePh != typePo))
		return false;

	if (kinematics[11]->isLorentzViolatingMonochromatic() && kinematics[-11]->isLorentzViolatingMonochromatic()) {
		auto kinEl = kinematics[ 11]->toLorentzViolatingKinematicsMonochromatic();
		auto kinPo = kinematics[-11]->toLorentzViolatingKinematicsMonochromatic();
		if (kinEl->getOrder() != kinPo->getOrder())
			return false;

		if (kinematics[22]->isLorentzViolatingMonochromatic()) {
			auto kinPh = kinematics[ 22]->toLorentzViolatingKinematicsMonochromatic();
			if (kinEl->getOrder() != kinPh->getOrder())
				return true;
		}
	} 
	
	return false;
}

string PhotonDecay::getInteractionTag() const {
	return interactionTag;
}

double PhotonDecay::computeThresholdMomentum() const {
	double pThr = std::numeric_limits<double>::max();
	if (! isImplemented())
		return pThr;

	double chiEl = 0;
	double chiPo = 0;
	double chiPh = 0;

	if (kinematics[22]->isLorentzViolatingMonochromatic())
		chiPh = kinematics[22]->toLorentzViolatingKinematicsMonochromatic()->getCoefficient();
	else if (kinematics[22]->isSpecialRelativistic())
		chiPh = 0;

	if (kinematics[11]->isLorentzViolatingMonochromatic())
		chiEl = kinematics[11]->toLorentzViolatingKinematicsMonochromatic()->getCoefficient();
	else if (kinematics[11]->isSpecialRelativistic())
		chiEl = 0;

	// current limitation
	chiPo = chiEl;

	double m = particleMasses.at(11);
	unsigned int order = kinematics[22]->toLorentzViolatingKinematicsMonochromatic()->getOrder();

	switch (order) {
		case 0:
			if (chiPh > chiEl) 
				pThr = 2. * mass_electron * c_light / sqrt(chiPh - chiEl);
			break;
		case 1:
			if (chiPh >= 0.) 
				pThr = cbrt(8 * pow_integer<2>(mass_electron * c_squared) * energy_planck / (2 * chiPh - chiEl)) / c_light;
			break;
		default:	
			throw std::runtime_error("PhotonDecay: only LIV of orders 0, 1 are implemented.");
	}

	return pThr;
}

double PhotonDecay::computeThresholdEnergy() const {
	int id = 22;
	double p = computeThresholdMomentum();
	return kinematics[id]->computeEnergyFromMomentum(p, id);
}

void PhotonDecay::process(Candidate* candidate) const {
	// check if photon
	int id = candidate->current.getId();
	if (id != 22)
		return;

	// get and scale the particle energy
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	// compute threshold energy
	double Ethr = computeThresholdEnergy() * (1 + z);

	// photon emission if above threshold
	if (E > Ethr) {
		candidate->setActive(false);

		if (haveElectrons) {
			double Epair = E - Ethr;
			Vector3d pos = candidate->current.getPosition();
			candidate->addSecondary( 11, 0.5 * Epair / (1 + z), pos);
			candidate->addSecondary(-11, 0.5 * Epair / (1 + z), pos);
		}
	}

	return;
}


} // namespace livpropa
