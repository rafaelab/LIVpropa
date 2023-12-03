#include "livpropa/VacuumCherenkov.h"


namespace livpropa {



VacuumCherenkov::VacuumCherenkov(ref_ptr<Kinematics> kinematics, bool havePhotons, double thinning, double limit) {
	setHavePhotons(havePhotons);
	setLimit(limit);
	setThinning(thinning);
	setKinematics(kinematics);
}

void VacuumCherenkov::setKinematics(ref_ptr<Kinematics> kin) {
	kinematics = kin;
}

void VacuumCherenkov::setHavePhotons(bool photons) {
	havePhotons = photons;
}

void VacuumCherenkov::setLimit(double l) {
	limit = l;
}

void VacuumCherenkov::setThinning(double t) {
	thinning = t;
}

void VacuumCherenkov::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

std::string VacuumCherenkov::getInteractionTag() const {
	return interactionTag;
}

double VacuumCherenkov::computeThresholdMomentum(const int& id) const {
	double m = particleMasses.at(id);

	// check type of kinematics
	std::string livType = kinematics->getShortIdentifier();

	double pThr = 0;
	if (livType == "LIV") {
		// type casting
		MonochromaticLIV* kin = static_cast<MonochromaticLIV*>(kinematics.get()); 

		// LIV order
		unsigned int n = kin->getOrder();

		// get coefficients
		double chiOt = kin->getCoefficientForParticle(id);
		double chiPh = kin->getCoefficientForParticle(22);
		
		if (n == 0) {
			if (chiOt > chiPh) 
				pThr = m * c_light / sqrt(chiOt - chiPh);

		} else if (n == 1) {
			if (chiOt > 0. and chiPh >= -3. * chiOt) {
				pThr = cbrt(pow_integer<2>(mass_electron * c_squared) * energy_planck / (2 * chiOt)) / c_light;
			} else if (chiOt <= 0 and chiOt > chiPh and chiPh <= -3. * chiOt) {
				pThr = cbrt(- 4. * pow_integer<2>(mass_electron * c_squared) * energy_planck * (chiOt + chiPh) / pow_integer<2>(chiOt - chiPh)) / c_light;
			}

		} else {
			throw std::runtime_error("VacuumCherenkov: only LIV of orders 0 and 1 are implemented.");
		}
	}

	return pThr;
}

double VacuumCherenkov::computeThresholdEnergy(const int& id) const {
	double p = computeThresholdMomentum(id);
	return kinematics->computeEnergyFromMomentum(p, id);
}

void VacuumCherenkov::process(Candidate* candidate) const {
		// check if photon
	int id = candidate->current.getId();
	if (abs(id) != 11)
		return;

	// get and scale the particle energy
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	// compute threshold energy
	double Ethr = computeThresholdEnergy(id) * (1 + z);

	// photon emission if above threshold
	if (E > Ethr) {
		candidate->current.setEnergy(Ethr / (1 + z));
		
		if (havePhotons) {
			double Ephoton = E - Ethr;
			Vector3d pos = candidate->current.getPosition();
			candidate->addSecondary(22, Ephoton / (1 + z), pos);
		}
	}

	return;
}


} // namespace livpropa
