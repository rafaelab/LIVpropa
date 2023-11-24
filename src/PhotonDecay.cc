#include "livpropa/PhotonDecay.h"


namespace livpropa {



PhotonDecay::PhotonDecay(ref_ptr<Kinematics> kinematics, bool haveElectrons, double thinning, double limit) {
	setHaveElectrons(haveElectrons);
	setLimit(limit);
	setThinning(thinning);
	setKinematics(kinematics);
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

double PhotonDecay::computeThresholdEnergy() const {
	return 0;
	// // LIV order
	// unsigned int n0 = livPhoton.getOrder() - 1;
	// unsigned int n1 = livPrimary.getOrder() - 1;
	// if (n0 != n1) 
	// 	throw std::runtime_error("PhotonDecay: the code can only compute cases where both the particle involved and the photon have the same LIV order.");
	// unsigned int n = n0;

	// // LIV coefficients
	// double chi0 = livPhoton.getCoefficient();
	// double chi1 = livPrimary.getCoefficient();
	
	// double pThr = 0;
	// double m = 0;
	// if (abs(livPrimary.getParticleId()) == 11)
	// 	m = mass_electron;
	// else
	// 	throw std::runtime_error("PhotonDecay: cannot compute energy threshold for particle other than electrons/positrons.");

	// if (n == 0) {
	// 	pThr = 2. * m * c_light / sqrt(chi1 - chi0);
	// } else if (n == 1) {
	// 	double p = m * c_squared * energy_planck;
	// 	if (chi1 > 0. && chi0 >= -3. * chi1) {
	// 		pThr = cbrt(8. * p / (2 * chi0 - chi1)) / c_light;
	// 	} else if (chi0 < chi1 && chi1 <= 0. && chi0 < -3. * chi1) {
	// 		pThr = cbrt(- 8. * p * chi1 / pow_integer<2>(chi1 - chi0)) / c_light;
	// 	}
	// } else {
	// 	throw std::runtime_error("PhotonDecay: only LIV of orders 1 and 2 are implemented.");
	// }

	// double pc = pThr * c_light;

	// return sqrt(pow_integer<2>(pc) + pow_integer<2>(m * c_squared) + chi1 * pow_integer<2>(pc) * pow(pc / energy_planck, n));
}

void PhotonDecay::process(Candidate* candidate) const {
	// check if electron / positron
	int id = candidate->current.getId();

	// if (id != livPrimary.getParticleId() || livPhoton.getOrder() == 0)
	// 	return;

	// // scale the particle energy
	// double z = candidate->getRedshift();
	// double E = candidate->current.getEnergy() * (1 + z);

	
	// // compute threshold energy
	// double Ethr = computeThresholdEnergy() * (1 + z);

	// // photon emission if above threshold
	// if (E > Ethr) {
	// 	double Ephoton = E - Ethr;
	// 	Vector3d pos = candidate->current.getPosition();
	// 	candidate->addSecondary(22, Ephoton / (1 + z), pos);
	// } else {
	// 	return;
	// }
	
}


} // namespace livpropa
