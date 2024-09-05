#include "livpropa/PhotonDecay.h"


namespace livpropa {



PhotonDecay::PhotonDecay(ref_ptr<Kinematics> kinematics, bool haveElectrons, double thinning, double limit) {
	setHaveElectrons(haveElectrons);
	setLimit(limit);
	setThinning(thinning);
	setKinematics(kinematics);
	setInteractionTag("PD");
}

void PhotonDecay::setKinematics(ref_ptr<Kinematics> kin) {
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

void PhotonDecay::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string PhotonDecay::getInteractionTag() const {
	return interactionTag;
}

double PhotonDecay::computeThresholdMomentum() const {
	// check type of kinematics
	std::string livType = kinematics->getShortIdentifier();

	double pThr = std::numeric_limits<float>::max();
	if (livType == "LIV") {
		// type casting
		MonochromaticLIV* kin = static_cast<MonochromaticLIV*>(kinematics.get()); 

		// LIV order
		unsigned int n = kin->getOrder();

		// get coefficients
		double chiEl = kin->getCoefficientForParticle(11);
		double chiPh = kin->getCoefficientForParticle(22);
		
		if (n == 0) {
			if (chiPh > chiEl) 
				pThr = 2. * mass_electron * c_light / sqrt(chiPh - chiEl);
		} else if (n == 1) {
			if (chiPh >= 0.) 
				pThr = cbrt(8 * pow_integer<2>(mass_electron * c_squared) * energy_planck / (2 * chiPh - chiEl)) / c_light;
		} else {
			throw std::runtime_error("PhotonDecay: only LIV of orders 0 and 1 are implemented.");
		}
	}

	return pThr;
}

double PhotonDecay::computeThresholdEnergy() const {
	double p = computeThresholdMomentum();
	return kinematics->computeEnergyFromMomentum(p, 22);
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
