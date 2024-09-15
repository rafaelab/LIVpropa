#include "livpropa/VacuumCherenkov.h"


namespace livpropa {



VacuumCherenkov::VacuumCherenkov(int id, ref_ptr<AbstractKinematics> kinOt, ref_ptr<AbstractKinematics> kinPh, bool havePhotons, VacuumCherenkovSpectrum spec, ref_ptr<SamplerEvents> samplerEvents, ref_ptr<SamplerDistribution> samplerDistribution, int maxSamples, double limit) {
	setInteractionTag("VC");
	setParticle(id);
	setHavePhotons(havePhotons);
	setLimit(limit);
	setKinematicsPhoton(kinPh);
	setKinematicsParticle(kinOt);
	setSpectrum(spec);
	setSamplerEvents(samplerEvents);
	setSamplerDistribution(samplerDistribution);;
	setMaximumSamples(maxSamples);
}

VacuumCherenkov::VacuumCherenkov(int id, ref_ptr<AbstractKinematics> kin, bool havePhotons, VacuumCherenkovSpectrum spec, ref_ptr<SamplerEvents> samplerEvents, ref_ptr<SamplerDistribution> samplerDistribution, int maxSamples, double limit) {
	setInteractionTag("VC");
	setParticle(id);
	setHavePhotons(havePhotons);
	setLimit(limit);
	setKinematicsPhoton(kin);
	setKinematicsParticle(kin);
	setSpectrum(spec);
	setSamplerEvents(samplerEvents);
	setSamplerDistribution(samplerDistribution);
	setMaximumSamples(maxSamples);
}

void VacuumCherenkov::setParticle(int id) {
	particleId = id;
}

void VacuumCherenkov::setKinematicsPhoton(ref_ptr<AbstractKinematics> kin) {
	kinematicsPhoton = kin;
}

void VacuumCherenkov::setKinematicsParticle(ref_ptr<AbstractKinematics> kin) {
	kinematicsParticle = kin;
}

void VacuumCherenkov::setHavePhotons(bool photons) {
	havePhotons = photons;
}

void VacuumCherenkov::setLimit(double l) {
	limit = l;
}

void VacuumCherenkov::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

void VacuumCherenkov::setSpectrum(VacuumCherenkovSpectrum spec) {
	if (spec == VacuumCherenkovSpectrum::Default)
		spectrum = getDefaultSpectrum(kinematicsParticle);
	else
		spectrum = spec;

	string errorMesssage = "VacuumCherenkov: treatment for this particular combination of kinematics of photon + other particle is not implemented.";

	string kinType = kinematicsParticle->getNameTag();
	if ((kinType == "LIVmono0" or kinType == "LIVmono1") and spec != VacuumCherenkovSpectrum::Step) {
		throw runtime_error(errorMesssage);
	}

	if (spectrum == VacuumCherenkovSpectrum::Full) {
		distribution = buildSpectrum(kinematicsParticle, kinematicsPhoton);
	}
}

void VacuumCherenkov::setSamplerEvents(ref_ptr<SamplerEvents> s) {
	if (s == nullptr)
		samplerEvents = new SamplerEventsNull();
	else
		samplerEvents = s;
}

void VacuumCherenkov::setSamplerDistribution(ref_ptr<SamplerDistribution> s) {
	samplerDistribution = s;
}

void VacuumCherenkov::setMaximumSamples(int nSamples) {
	maximumSamples = nSamples;
}

int VacuumCherenkov::getParticle() const {
	return particleId;
}

string VacuumCherenkov::getInteractionTag() const {
	return interactionTag;
}

ref_ptr<Histogram1D> VacuumCherenkov::getDistribution() const {
	return distribution;
}

double VacuumCherenkov::computeThresholdMomentum() const {
	return thresholdMomentum(particleId, kinematicsParticle, kinematicsPhoton);
}

double VacuumCherenkov::computeThresholdEnergy() const {
	double pThr = thresholdMomentum(particleId, kinematicsParticle, kinematicsPhoton);
	return kinematicsParticle->computeEnergyFromMomentum(pThr, particleId);
}

double VacuumCherenkov::computeInteractionRate(const double& p) const {
	return interactionRate(p, kinematicsParticle, kinematicsPhoton);
}

void VacuumCherenkov::process(Candidate* candidate) const {
	// check if electron/positron (only particles implemented so far)
	int id = candidate->current.getId();

	if (id == 22 && particleId != id)
		return;

	// do not perform any calculations for the special-relativistic case
	if (kinematicsParticle->isLorentzInvariant())
		return;

	if (spectrum == VacuumCherenkovSpectrum::Absent) 
		return;

	// get and scale the particle energy
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	// compute threshold energy; emit photon only if E > Ethr
	double Ethr = computeThresholdEnergy() * (1 + z);
	if (E < Ethr)
		return;

		
	// cout << "process: E = " << E / eV << ", Ethr = " << Ethr / eV << ", dE = " << (E - Ethr) / eV << endl;		
	if (havePhotons) {
		// position where photon will be produced
		Vector3d pos = candidate->current.getPosition();
		double step = candidate->getCurrentStep() / (1 + z); 

		switch (spectrum) {
			case VacuumCherenkovSpectrum::Step: {
				emissionSpectrumStep(candidate, Ethr);
				break;
			}
			case VacuumCherenkovSpectrum::Full: {
				emissionSpectrumFull(candidate, Ethr);
				break;
			}
			default:
				break;
		}
	}

	candidate->current.setEnergy((E - Ethr) / (1 + z));

}

void VacuumCherenkov::emissionSpectrumStep(Candidate* candidate, const double& Ethr) const {
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	// // exactly one photon is emitted
	// double Ephoton = E - Ethr;	

	// // position where photon will be produced
	// Random &random = Random::instance();
	// Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
	
	// candidate->addSecondary(22, Ephoton / (1 + z), pos);
}

void VacuumCherenkov::emissionSpectrumFull(Candidate* candidate, const double& Ethr) const {
	int id = candidate->current.getId();
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);
	double step = candidate->getCurrentStep() / (1 + z);

	// double rate = computeInteractionRate(id, E);
	// if (rate <= 0)
	// 	return;

// 	auto distIt = distributions.find(id);
// 	if (distIt == distributions.end())
// 		return;
// 	const ref_ptr<Histogram1D>& distribution = (*distIt).second;



// 	Random &random = Random::instance();

// 	// change in differential rate
// 	double dE0 = E - Ethr; 

// // cout << "process: E = " << E / eV << ", Ethr = " << Ethr / eV << ", dE = " << (dE0) / eV << endl;

// 	if (samplerDistribution == nullptr) {
// 		double dE = dE0;
// 		while (dE > 0) {
// 			Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
// 			std::pair<double, double> range = xRange(kinematics[id], kinematics[22]);
// 			cout << "sampling" << endl;
// 			double x = distribution->getSampleInRange(range);
// 			cout << "differentialProbability" << endl;
// 			double dPdx = differentialProbability(x, kinematics[id], kinematics[22]);
// 			double Ephoton = x * E;

// 			cout << "E = " << E / eV << ", Ephoton = " << Ephoton / eV << endl;

// 			// if the energy drops below the threshold, then there is no emission
// 			Ephoton = std::min(Ephoton, E - Ethr);

// 			candidate->addSecondary(22, Ephoton / (1 + z), pos);
// 			dE -= Ephoton;
// 		}

// 	} else {
// 		samplerDistribution->transformToPDF();
// 		samplerDistribution->transformToCDF();
// 		std::vector<double> sampled = samplerDistribution->getSample(maximumSamples);
// 		double dEs = std::accumulate(sampled.begin(), sampled.end(), decltype(sampled)::value_type(0));
// 		if (samplerDistribution->getSize() > 0) {
// 			double w0 = dE0 / dEs;
// 			for (size_t i = 0; i < sampled.size(); i++) {
// 				Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
// 				double Es = sampled[i];
// 				double w = w0 * samplerEvents->computeWeight(-11, Es, Es / E, i);
// 				candidate->addSecondary(22, Es, pos, w);
// 			}
// 		}
// 		samplerDistribution->clear();
// 	}
}

VacuumCherenkovSpectrum VacuumCherenkov::getDefaultSpectrum(const ref_ptr<AbstractKinematics>& kin) {
	string kinType = kin->getNameTag();
	
	if (kinType == "SR") {
		return VacuumCherenkovSpectrum::Absent;
	} else if (kinType == "LIVmono0") {
		return VacuumCherenkovSpectrum::Step;
	} else if (kinType == "LIVmono1") {
		return VacuumCherenkovSpectrum::Step;
	} else if (kinType == "LIVmono2") {
		return VacuumCherenkovSpectrum::Full;
	} 
	
	throw runtime_error("VacuumCherenkov: default spectrum can be retrieved only for monochromatic LIV with orders 0, 1, and 2.");	
}

template<class KO, class KP>
double VacuumCherenkov::thresholdMomentum(const int& id, const KO& kinOt, const KP& kinPh) {
	// throw runtime_error("VacuumCherenkov treatment for this particular combination of kinematics of photon + other particle is not implemented.");
	return _defaultThresholdMomentum;
}

template<>
double VacuumCherenkov::thresholdMomentum(const int& id, const ref_ptr<AbstractKinematics>& kinOt, const ref_ptr<AbstractKinematics>& kinPh) {
	string kOt = kinOt->getNameTag();
	string kPh = kinPh->getNameTag();

	if (kOt == "LIVmono0") {
		const auto& kinOtNew = kinOt->toMonochromaticLorentzViolatingKinematics0();
		if (kPh == "SR") {
			const SpecialRelativisticKinematics& kinPhNew = kinPh->toSpecialRelativisticKinematics();
			return thresholdMomentum(id, kinOtNew, kinPhNew);		
		} else if (kPh == "LIVmono0") {
			const MonochromaticLorentzViolatingKinematics<0>& kinPhNew = kinPh->toMonochromaticLorentzViolatingKinematics0();
			return thresholdMomentum(id, kinOtNew, kinPhNew);	
		}
	} else if (kOt == "LIVmono1") {
		const auto& kinOtNew = kinOt->toMonochromaticLorentzViolatingKinematics1();
		if (kPh == "SR") {
			const SpecialRelativisticKinematics& kinPhNew = kinPh->toSpecialRelativisticKinematics();
			return thresholdMomentum(id, kinOtNew, kinPhNew);		
		} else if (kPh == "LIVmono1") {
			const MonochromaticLorentzViolatingKinematics<1>& kinPhNew = kinPh->toMonochromaticLorentzViolatingKinematics1();
			return thresholdMomentum(id, kinOtNew, kinPhNew);	
		}
	} else if (kOt == "LIVmono2") {
		const auto& kinOtNew = kinOt->toMonochromaticLorentzViolatingKinematics2();
		if (kPh == "SR") {
			const SpecialRelativisticKinematics& kinPhNew = kinPh->toSpecialRelativisticKinematics();
			return thresholdMomentum(id, kinOtNew, kinPhNew);		
		} else if (kPh == "LIVmono2") {
			const MonochromaticLorentzViolatingKinematics<2>& kinPhNew = kinPh->toMonochromaticLorentzViolatingKinematics2();
			return thresholdMomentum(id, kinOtNew, kinPhNew);	
		}
	}

	return _defaultThresholdMomentum;
}

template<>
double VacuumCherenkov::thresholdMomentum(const int& id, const MonochromaticLorentzViolatingKinematics<0>& kinOt, const MonochromaticLorentzViolatingKinematics<0>& kinPh) {
	double pThr = _defaultThresholdMomentum;

	double mass = particleMasses.at(id);

	double chiOt = kinOt.getCoefficient();
	double chiPh = kinPh.getCoefficient();
	if (chiOt > chiPh) 
		pThr = mass * c_light / sqrt(chiOt - chiPh);

	return pThr;
}

template<>
double VacuumCherenkov::thresholdMomentum(const int& id, const MonochromaticLorentzViolatingKinematics<1>& kinOt, const MonochromaticLorentzViolatingKinematics<1>& kinPh) {
	double pThr = _defaultThresholdMomentum;

	double chiOt = kinOt.getCoefficient();
	double chiPh = kinPh.getCoefficient();

	double mass = particleMasses.at(id);
	double m2 = pow_integer<2>(mass * c_squared);

	if (chiOt > 0. and chiPh >= -3. * chiOt) {
		pThr = cbrt(m2 * energy_planck / (2 * chiOt)) / c_light;
	} else if (chiOt <= 0 and chiOt > chiPh and chiPh < -3. * chiOt) {
		pThr = cbrt(- 4. * m2 * energy_planck * (chiOt + chiPh) / pow_integer<2>(chiOt - chiPh)) / c_light;
	}

	return pThr;
}

template<>
double VacuumCherenkov::thresholdMomentum(const int& id, const MonochromaticLorentzViolatingKinematics<2>& kinOt,  const MonochromaticLorentzViolatingKinematics<2>& kinPh)  {
	double pThr = _defaultThresholdMomentum;

	double chiOt = kinOt.getCoefficient();
	double chiPh = kinPh.getCoefficient();

	double mass = particleMasses.at(id);
	double m2 = pow_integer<2>(mass * c_squared);
	double Epl2 = pow_integer<2>(energy_planck);
	double a = - 8 - 6 * sqrt(2);

	if (chiOt > 0. and chiPh >= a * chiOt) {
		pThr = pow(m2 * Epl2 / 3. / (chiOt * chiOt), 0.25);
	} else if ((a * chiOt <= 0 and chiPh < a * chiOt) or (chiOt <= 0 and chiPh < chiOt)) {
		double t = (chiOt + 2 * chiPh) / (chiOt - chiPh);
		double l = chiOt - chiPh;
		double F = 2. / 27. * l * (t * t * t + pow_integer<3>(sqrt(t * t - 3)) - 4.5 * t);
		pThr = pow(pow_integer<2>(m2 * energy_planck) / F, 0.25);
	}
	pThr /= c_light;

	return pThr;
}

template<class KO, class KP>
ref_ptr<Histogram1D> VacuumCherenkov::buildSpectrum(const KO& kinOt, const KP& kinPh) {
	throw runtime_error("Full treatment of vacuum Cherenkov spectrum for this particular combination of kinematics of photon + other particle is not implemented.");
}

template<>
ref_ptr<Histogram1D> VacuumCherenkov::buildSpectrum(const ref_ptr<AbstractKinematics>& kinOt, const ref_ptr<AbstractKinematics>& kinPh) {
	string kOt = kinOt->getNameTag();
	string kPh = kinPh->getNameTag();

	if (kOt == "LIVmono0") {
		const auto& kinOtNew = kinOt->toMonochromaticLorentzViolatingKinematics0();
		if (kPh == "SR") {
			const SpecialRelativisticKinematics& kinPhNew = kinPh->toSpecialRelativisticKinematics();
			return buildSpectrum(kinOtNew, kinPhNew);		
		} else if (kPh == "LIVmono0") {
			const MonochromaticLorentzViolatingKinematics<0>& kinPhNew = kinPh->toMonochromaticLorentzViolatingKinematics0();
			return buildSpectrum(kinOtNew, kinPhNew);	
		}
	} else if (kOt == "LIVmono1") {
		const auto& kinOtNew = kinOt->toMonochromaticLorentzViolatingKinematics1();
		if (kPh == "SR") {
			const SpecialRelativisticKinematics& kinPhNew = kinPh->toSpecialRelativisticKinematics();
			return buildSpectrum(kinOtNew, kinPhNew);		
		} else if (kPh == "LIVmono1") {
			const MonochromaticLorentzViolatingKinematics<1>& kinPhNew = kinPh->toMonochromaticLorentzViolatingKinematics1();
			return buildSpectrum(kinOtNew, kinPhNew);	
		}
	} else if (kOt == "LIVmono2") {
		const auto& kinOtNew = kinOt->toMonochromaticLorentzViolatingKinematics2();
		if (kPh == "SR") {
			const SpecialRelativisticKinematics& kinPhNew = kinPh->toSpecialRelativisticKinematics();
			return buildSpectrum(kinOtNew, kinPhNew);		
		} else if (kPh == "LIVmono2") {
			const MonochromaticLorentzViolatingKinematics<2>& kinPhNew = kinPh->toMonochromaticLorentzViolatingKinematics2();
			return buildSpectrum(kinOtNew, kinPhNew);	
		}
	}

	throw runtime_error("Full treatment of vacuum Cherenkov spectrum for this particular combination of kinematics of photon + other particle is not implemented.");
}

template<class KP>
ref_ptr<Histogram1D> VacuumCherenkov::buildSpectrum(const MonochromaticLorentzViolatingKinematics<2>& kinOt, const KP& kinPh) {
	ref_ptr<Histogram1D> dist = new Histogram1D(100, 0., 1.);

	double chiOt = kinOt.getCoefficient();
	double chiPh = kinPh.getCoefficient();

	double chiOt2 = chiOt * chiOt;
	double chiOt3 = chiOt2 * chiOt;
	double chiPh2 = chiPh * chiPh;
	double dChi = chiPh - chiOt;
	double S = sqrt(3. * chiOt * (4 * chiPh - chiOt));

	double Gp0 = chiOt * (S - 3 * chiOt) / 160. / pow_integer<4>(dChi);
	double Gp1 = 37. * (S - 6 * chiPh) * chiPh * chiOt - 64. * pow_integer<3>(chiOt);
	double Gp2 = - pow_integer<2>(chiOt) * (14 * S - 207. * chiPh) - 10 * (5 * S - 16 * chiPh) * pow_integer<2>(chiPh);
	double Gp = Gp0 * (Gp1 + Gp2);

	double Gm0 = chiOt * (- S - 3 * chiOt) / 160. / pow_integer<4>(dChi);
	double Gm1 = 37. * (- S - 6 * chiPh) * chiPh * chiOt - 64. * pow_integer<3>(chiOt);
	double Gm2 = - pow_integer<2>(chiOt) * (- 14 * S - 207. * chiPh) - 10 * (- 5 * S - 16 * chiPh) * pow_integer<2>(chiPh);
	double Gm = Gm0 * (Gm1 + Gm2);

	double P = 0;
	for (size_t i = 0; i < dist->getNumberOfBins(); i++) {
		double x = dist->getBinCentre(i);
		double x2 = x * x;
		double x3 = x2 * x;
		double y = x3 - 3. * x2 + 3. * x;
		double a =  (2. / x - 2. + x);
		double b1 = - chiPh * x3 / 2.;
		double b2 = chiOt * y / 2.;
		double b = b1 + b2;
		double c = (157. * chiOt - 22. * chiPh) / 120.;
		if (chiOt >= std::max(0., chiPh)) 
			P = a * b / c;
		else if (chiOt > 0 and  chiOt < chiPh)
			P = std::max(0., a * b / (c - Gm));
		else if (chiOt < chiPh and chiPh < 0.)
			P = std::max(0., a * b / Gp);
		else
			P = 0;
		dist->setBinContent(i, P);
	}

	dist->normalise(dist->integrate());
	dist->transformToPDF();
	dist->transformToCDF();

	return dist;
}

template<class KO, class KP>
double VacuumCherenkov::interactionRate(const double& p, const KO& kinOt, const KP& kinPh) {
	return _defaultInteractionRate;
}

template<>
double VacuumCherenkov::interactionRate(const double& p, const ref_ptr<AbstractKinematics>& kinOt, const ref_ptr<AbstractKinematics>& kinPh) {
	string kOt = kinOt->getNameTag();
	string kPh = kinPh->getNameTag();

	if (kOt == "LIVmono2") {
		const MonochromaticLorentzViolatingKinematics<2>& kinOtNew = kinOt->toMonochromaticLorentzViolatingKinematics2();
		const MonochromaticLorentzViolatingKinematics<2>& kinPhNew = kinPh->toMonochromaticLorentzViolatingKinematics2();
		return interactionRate(p, kinOtNew, kinPhNew);
	}

	return _defaultInteractionRate;
}

template<>
double VacuumCherenkov::interactionRate(const double& p, const MonochromaticLorentzViolatingKinematics<2>& kinOt, const MonochromaticLorentzViolatingKinematics<2>& kinPh) {
	double chiOt = kinOt.getCoefficient();
	double chiPh = kinPh.getCoefficient();

	double q = pow_integer<3>(p * c_light) / pow_integer<2>(energy_planck) / h_dirac;
	double Q = alpha_finestructure * q;
	double dChi = chiPh - chiOt;
	double S = sqrt(3. * chiOt * (4 * chiPh - chiOt));

	if (chiOt == 0) {
		if (chiOt > chiPh)
			return (157. * chiOt - 22. * chiPh) / 120. * Q;
	} else if (chiOt > 0) {
		if (chiOt < chiPh) {
			double Gp0 = chiOt * (S - 3 * chiOt) / 160. / pow_integer<4>(dChi);
			double Gp1 = 37. * (S - 6 * chiPh) * chiPh * chiOt - 64. * pow_integer<3>(chiOt);
			double Gp2 = - pow_integer<2>(chiOt) * (14 * S - 207. * chiPh) - 10 * (5 * S - 16 * chiPh) * pow_integer<2>(chiPh);
			double Gp = Gp0 * (Gp1 + Gp2);
			return Q * Gp;
		} else {
			return (157. * chiOt - 22. * chiPh) / 120. * Q;
		}
	} else {
		if (chiOt > chiPh) {
			double Gm0 = chiOt * (- S - 3 * chiOt) / 160. / pow_integer<4>(dChi);
			double Gm1 = 37. * (- S - 6 * chiPh) * chiPh * chiOt - 64. * pow_integer<3>(chiOt);
			double Gm2 = - pow_integer<2>(chiOt) * (- 14 * S - 207. * chiPh) - 10 * (- 5 * S - 16 * chiPh) * pow_integer<2>(chiPh);
			double Gm = Gm0 * (Gm1 + Gm2);
			return Q * ((157. * chiOt - 22. * chiPh) / 120. - Gm);
		}
	}

	return _defaultInteractionRate;;
}

template<>
double VacuumCherenkov::interactionRate(const double& p, const MonochromaticLorentzViolatingKinematics<2>& kinOt, const SpecialRelativisticKinematics& kinPh) {
	const MonochromaticLorentzViolatingKinematics<2> kinPhNew(0.);
	return interactionRate(p, kinOt, kinPhNew);
}

// template<class KO, class KP>
// double VacuumCherenkov::differentialInteractionRate(const double& p, const double& x, const KO& kinOt, const KP& kinPh) {
// 	return _defaultInteractionRate;
// }

// template<>
// double VacuumCherenkov::differentialInteractionRate(const double& p, const double& x, const MonochromaticLorentzViolatingKinematics<2>& kinOt, const MonochromaticLorentzViolatingKinematics<2>& kinPh) {
// 	double chiOt = kinOt.getCoefficient();
// 	double chiPh = kinPh.getCoefficient();

// 	double q = pow_integer<3>(p * c_light) / pow_integer<2>(energy_planck) / h_dirac;
// 	double Q = alpha_finestructure * q;
// 	double dChi = chiPh - chiOt;

// 	double y = (x * x * x - 3. * x * x + 3. * x);
// 	double omega = - chiPh * Q / 2. + chiOt * Q / 2. y; 

// 	return omega * (2. / x - 2. + x);
// }

// template<class KO, class KP>
// double VacuumCherenkov::differentialProbability(const double& x, const KO& kinOt, const KP& kinPh) {
// 	// throw runtime_error("Full treatment of vacuum Cherenkov spectrum for this particular combination of kinematics of photon + other particle is not implemented.");
// 	return 0;
// }

// template<>
// double VacuumCherenkov::differentialProbability(const double& x, const MonochromaticLorentzViolatingKinematics<2>& kinOt, const  MonochromaticLorentzViolatingKinematics<2>& kinPh) {
// 	double chiOt = kinOt.getCoefficient();
// 	double chiPh = kinPh.getCoefficient();

// 	double chiOt2 = chiOt * chiOt;
// 	double chiOt3 = chiOt2 * chiOt;
// 	double chiPh2 = chiPh * chiPh;
// 	double dChi = chiPh - chiOt;
// 	double S = sqrt(3. * chiOt * (4 * chiPh - chiOt));

// 	double Gp0 = chiOt * (S - 3 * chiOt) / 160. / pow_integer<4>(dChi);
// 	double Gp1 = 37. * (S - 6 * chiPh) * chiPh * chiOt - 64. * pow_integer<3>(chiOt);
// 	double Gp2 = - pow_integer<2>(chiOt) * (14 * S - 207. * chiPh) - 10 * (5 * S - 16 * chiPh) * pow_integer<2>(chiPh);
// 	double Gp = Gp0 * (Gp1 + Gp2);

// 	double Gm0 = chiOt * (- S - 3 * chiOt) / 160. / pow_integer<4>(dChi);
// 	double Gm1 = 37. * (- S - 6 * chiPh) * chiPh * chiOt - 64. * pow_integer<3>(chiOt);
// 	double Gm2 = - pow_integer<2>(chiOt) * (- 14 * S - 207. * chiPh) - 10 * (- 5 * S - 16 * chiPh) * pow_integer<2>(chiPh);
// 	double Gm = Gm0 * (Gm1 + Gm2);

// 	double x2 = x * x;
// 	double x3 = x2 * x;
// 	double y = x3 - 3. * x2 + 3. * x;
// 	double a =  (2. / x - 2. + x);
// 	double b1 = - chiPh * x3 / 2.;
// 	double b2 = chiOt * y / 2.;
// 	double b = b1 + b2;
// 	double c = (157. * chiOt - 22. * chiPh) / 120.;

// 	double P = 0;
// 	if (chiOt >= std::max(0., chiPh)) 
// 		P = a * b / c;
// 	else if (chiOt > 0 and  chiOt < chiPh)
// 		P = std::max(0., a * b / (c - Gm));
// 	else if (chiOt < chiPh and chiPh < 0.)
// 		P = std::max(0., a * b / Gp);
// 	else
// 		P = 0;
	
// 	return P;
// }

// template<>
// double VacuumCherenkov::differentialProbability(const double& x, const MonochromaticLorentzViolatingKinematics<2>& kinOt, const  SpecialRelativisticKinematics& kinPh) {
// 	const MonochromaticLorentzViolatingKinematics<2> kinPhNew(0.);
// 	return differentialProbability(x, kinOt, kinPhNew);
// }

// template<class KO, class KP>
// std::pair<double, double> VacuumCherenkov::xRange(const KO& kinOt, const KP& kinPh) {
// 	return std::pair<double, double>(1., 1.);
// }

// template<>
// std::pair<double, double> VacuumCherenkov::xRange(const MonochromaticLorentzViolatingKinematics<2>& kinOt, const MonochromaticLorentzViolatingKinematics<2>& kinPh)  {
// 	std::pair<double, double> x = std::make_pair(1, 1);

// 	double chiOt = kinOt.getCoefficient();
// 	double chiPh = kinPh.getCoefficient();

// 	double xMin = 1;
// 	double xMax = 1;

// 	if (chiOt >= 0 and chiPh <= chiOt) {
// 		xMin = 0;
// 		xMax = 1;
// 	} else if (chiOt > 0 and chiPh > chiOt) {
// 		xMin = 0;
// 		if (chiOt = chiPh and chiPh != 0)
// 			xMax = 1;
// 		else
// 			xMax = -1.5 * chiOt / (chiPh - chiOt) + 0.5 * sqrt(3. * chiOt * (4. * chiPh - chiOt)) / (chiPh - chiOt);
// 	} else if (chiOt < 0 and chiPh < chiOt) {
// 		if (chiOt = chiPh and chiPh != 0)
// 			xMin = 1;
// 		else
// 			xMin = -1.5 * chiOt / (chiPh - chiOt) + 0.5 * sqrt(3. * chiOt * (4. * chiPh - chiOt)) / (chiPh - chiOt);
// 		xMax = 1;
// 	} else {
// 		xMin = -1;
// 		xMax = -1;
// 	}

// 	return std::make_pair(xMin, xMax);
// }



} // namespace livpropa
