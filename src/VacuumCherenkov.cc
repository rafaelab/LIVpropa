#include "livpropa/VacuumCherenkov.h"


namespace livpropa {



VacuumCherenkov::VacuumCherenkov(int id, KinematicsMap kin, VacuumCherenkovSpectrum spec, bool havePhotons, bool angularCorrection, bool continuousEnergyLoss, ref_ptr<Sampler> sampler, ref_ptr<RuntimeWeighter> runtimeWeighter, ref_ptr<PosterioriWeighter> posterioriWeighter, double limit) {
	setInteractionTag("VC");
	setParticle(id);
	setHavePhotons(havePhotons);
	setAngularCorrection(angularCorrection);
	setContinuousEnergyLoss(continuousEnergyLoss);
	setLimit(limit);

	if (not kin.exists(id)) {
		throw runtime_error("VacuumCherenkov: kinematics for the desired particle is not specified in `KinematicsMap`.");
	}
	setKinematicsParticle(kin[id]);

	if (not kin.exists(22)) {
		KISS_LOG_WARNING << "`VacuumCherenkov` is being iniatialised from `KinematicsMap` without specifying the photon kinematics. It will be set to the same as the desired particle kinematics, by default." << endl;
		setKinematicsPhoton(kin[id]);
	} else {
		setKinematicsPhoton(kin[22]);
	}

	setSampler(sampler);
	setRuntimeWeighter(runtimeWeighter);
	setPosterioriWeighter(posterioriWeighter);

	setMinimumEnergyFractionSpectrum(minEnergyFractionSpectrum);
	setNumberOfBinsSpectrum(nBinsSpectrum);
	setSpectrum(spec, sampler);
}

VacuumCherenkov::VacuumCherenkov(int id, ref_ptr<Kinematics> kinOt, ref_ptr<Kinematics> kinPh, VacuumCherenkovSpectrum spec, bool havePhotons, bool angularCorrection, bool continuousEnergyLoss, ref_ptr<Sampler> sampler, ref_ptr<RuntimeWeighter> runtimeWeighter, ref_ptr<PosterioriWeighter> posterioriWeighter, double limit) {
	setInteractionTag("VC");
	setParticle(id);
	setHavePhotons(havePhotons);
	setAngularCorrection(angularCorrection);
	setLimit(limit);
	setContinuousEnergyLoss(continuousEnergyLoss);

	setKinematicsPhoton(kinPh);
	setKinematicsParticle(kinOt);

	setSampler(sampler);
	setRuntimeWeighter(runtimeWeighter);
	setPosterioriWeighter(posterioriWeighter);

	setMinimumEnergyFractionSpectrum(minEnergyFractionSpectrum);
	setNumberOfBinsSpectrum(nBinsSpectrum);
	setSpectrum(spec, sampler);
}

VacuumCherenkov::VacuumCherenkov(int id, ref_ptr<Kinematics> kin, VacuumCherenkovSpectrum spec, bool havePhotons, bool angularCorrection, bool continuousEnergyLoss, ref_ptr<Sampler> sampler, ref_ptr<RuntimeWeighter> runtimeWeighter, ref_ptr<PosterioriWeighter> posterioriWeighter, double limit) {
	setInteractionTag("VC");
	setParticle(id);
	setHavePhotons(havePhotons);
	setAngularCorrection(angularCorrection);
	setContinuousEnergyLoss(continuousEnergyLoss);
	setLimit(limit);

	setKinematicsPhoton(kin);
	setKinematicsParticle(kin);

	setSampler(sampler);
	setRuntimeWeighter(runtimeWeighter);
	setPosterioriWeighter(posterioriWeighter);

	setMinimumEnergyFractionSpectrum(minEnergyFractionSpectrum);
	setNumberOfBinsSpectrum(nBinsSpectrum);
	setSpectrum(spec, sampler);
}

void VacuumCherenkov::setParticle(int id) {
	if (id == 22) {
		throw runtime_error("VacuumCherenkov: cannot set photon as the particle.");
	}
	particleId = id;
}

void VacuumCherenkov::setKinematicsPhoton(ref_ptr<Kinematics> kin) {
	kinematicsPhoton = kin;
}

void VacuumCherenkov::setKinematicsParticle(ref_ptr<Kinematics> kin) {
	kinematicsParticle = kin;
}

void VacuumCherenkov::setRuntimeWeighter(ref_ptr<RuntimeWeighter> w) {
	runtimeWeighter = w;
	runtimeWeighter->reset();
}

void VacuumCherenkov::setPosterioriWeighter(ref_ptr<PosterioriWeighter> w) {
	posterioriWeighter = w;
	posterioriWeighter->reset();
}

void VacuumCherenkov::setSampler(ref_ptr<Sampler> s) {
	if (s == nullptr) {
		sampler = new InverseSampler();
	} else {
		sampler = s;
	}
}

void VacuumCherenkov::setHavePhotons(bool photons) {
	havePhotons = photons;
}

void VacuumCherenkov::setAngularCorrection(bool correction) {
	angularCorrection = correction;
}

void VacuumCherenkov::setContinuousEnergyLoss(bool loss) {
	continuousEnergyLoss = loss;
}

void VacuumCherenkov::setLimit(double l) {
	limit = l;
}

void VacuumCherenkov::setInteractionTag(std::string tag) {
	interactionTag = tag;
}

void VacuumCherenkov::setSpectrum(VacuumCherenkovSpectrum spec, ref_ptr<Sampler> s) {
	string kinType = kinematicsParticle->getNameTag();

	if (spec == VacuumCherenkovSpectrum::Default) {
		if (kinType == "SR") {
			spectrum = VacuumCherenkovSpectrum::Absent;
		} else if (kinType == "LIVmono0" or kinType == "LIVmono1") {
			spectrum = VacuumCherenkovSpectrum::Step;
		} else if (kinType == "LIVmono2") {
			spectrum = VacuumCherenkovSpectrum::Full;
		} else {
			throw runtime_error("VacuumCherenkov: treatment for this particular combination of kinematics of photon + other particle is not implemented.");
		}
	} else {
		spectrum = spec;
	}

	string errorMessage = "VacuumCherenkov: treatment for this particular combination of kinematics of photon + other particle is not implemented.";
	if ((kinType == "LIVmono0" or kinType == "LIVmono1") and spectrum != VacuumCherenkovSpectrum::Step) {
		throw runtime_error(errorMessage);
	}

	if (spectrum == VacuumCherenkovSpectrum::Full) {
		if (sampler == nullptr) 
			sampler = new NestedSampler(1000);
		buildSpectrum(kinematicsParticle, kinematicsPhoton);
	}
}

void VacuumCherenkov::setMinimumEnergyFractionSpectrum(double xmin) {
	xMin = xmin;
}

void VacuumCherenkov::setNumberOfBinsSpectrum(unsigned int n) {
	nBins = n;
}

int VacuumCherenkov::getParticle() const {
	return particleId;
}

string VacuumCherenkov::getInteractionTag() const {
	return interactionTag;
}

ref_ptr<Kinematics> VacuumCherenkov::getKinematicsParticle() const {
	return kinematicsParticle;
}

ref_ptr<Kinematics> VacuumCherenkov::getKinematicsPhoton() const {
	return kinematicsPhoton;
}

ref_ptr<Histogram1D> VacuumCherenkov::getDistribution() const {
	return distribution;
}

ref_ptr<Sampler> VacuumCherenkov::getSampler() const {
	return sampler;
}

ref_ptr<RuntimeWeighter> VacuumCherenkov::getRuntimeWeighter() const {
	return runtimeWeighter;
}

ref_ptr<PosterioriWeighter> VacuumCherenkov::getPosterioriWeighter() const {
	return posterioriWeighter;
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
	if (particleId != id)
		return;

	// do not perform any calculations for the special-relativistic case
	if (kinematicsParticle->isLorentzInvariant())
		return;

	// get and scale the particle energy
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	// compute threshold energy; emit photon only if E > Ethr
	double Ethr = computeThresholdEnergy() * (1 + z);

	if (E < Ethr or isnan(Ethr) or isinf(Ethr))
		return;

	switch (spectrum) {
		case VacuumCherenkovSpectrum::Step: 
			emissionSpectrumStep(candidate, Ethr);
			break;
		case VacuumCherenkovSpectrum::Full:
			emissionSpectrumFull(candidate, Ethr);
			break;
		default:
			break;
	}
}

void VacuumCherenkov::emissionSpectrumStep(Candidate* candidate, const double& Ethr) const {
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	// exactly one photon is emitted
	double Ephoton = E - Ethr;	

	// position where photon will be produced
	Random &random = Random::instance();
	Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
	
	candidate->addSecondary(22, Ephoton / (1 + z), pos);
	candidate->current.setEnergy(Ethr / (1 + z));
}

void VacuumCherenkov::emissionSpectrumFull(Candidate* candidate, const double& Ethr) const {
	Random& random = Random::instance();
	int id = candidate->current.getId();
	double z = candidate->getRedshift();
	double step = candidate->getCurrentStep() / (1 + z);
	double E = candidate->current.getEnergy() * (1 + z);

	double p = kinematicsParticle->computeMomentumFromEnergy(E, id);
	if (p <= 0) {
		KISS_LOG_WARNING << "VacuumCherenkov: particle momentum is negative." << endl;
		return;
	}

	std::pair<double, double> range = xRange(kinematicsParticle, kinematicsPhoton);
	if (range.first < 0. or range.second < 0. or range.second > 1) {
		KISS_LOG_WARNING << "VacuumCherenkov: particle momentum is negative." << endl;
		return;
	}
	
	// define "threshold" energy; it is different depending on whethr  CEL is assumed
	double Et = Ethr;

	// two implementations: continuous energy loss and instantaneous energy loss (the latter is faster)
	// the CEL approach has not been thoroughly tested
	if (continuousEnergyLoss) {
		double rate = computeInteractionRate(E);
		if (rate <= 0)
			return;

		double lossLength = c_light / rate;
		double loss = step / lossLength; // relative energy loss
		double dE = E * loss;

		if (E - dE < Ethr)
			return;

		// within one step, we can think as if the threshold changed
		Et = std::max(E - dE, Ethr);
	} 

	// change variable name
	double Eo = E;

	// current position; photons may be emitted at different positions
	Vector3d pos = candidate->current.getPosition();


	do {
		// perform sampling of the energy fraction with a specific sampler
		std::pair<double, double> range = xRange(kinematicsParticle, kinematicsPhoton);
		std::pair<double, double> sample = sampler->getSample(random);
		double x = sample.first;
		double w = sample.second;

		// compute photon energy
		double Ephoton = x * Eo;

		if (havePhotons and w > 0) {
			// decide whether to weight the secondaries at runtime or a posteriori
			if (posterioriWeighter != nullptr) {
				posterioriWeighter->push(Ephoton / (1 + z), w);
			} else {
				if (continuousEnergyLoss)
					pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
				candidate->addSecondary(22, Ephoton / (1 + z), pos, w);
			}
		} 
		
		Eo -= Ephoton;
	} while (Eo > Et);

	candidate->current.setEnergy(Eo / (1 + z));


	// add photon spectrum if using a posteriori weighter
	if (havePhotons and posterioriWeighter != nullptr) {
		vector<std::pair<double, double>> samples = posterioriWeighter->getEvents(random);

		for (auto sample : samples) {
			if (continuousEnergyLoss)
				Vector3d pos = random.randomInterpolatedPosition(candidate->previous.getPosition(), candidate->current.getPosition());
			if (! isnan(sample.second) and sample.second > 0)
				candidate->addSecondary(22, sample.first / (1 + z), pos, sample.second);
		}

		// reset the weighter
		posterioriWeighter->reset();
	}

	// reset the runtime weighters
	if (runtimeWeighter != nullptr)
		runtimeWeighter->reset();
}

VacuumCherenkovSpectrum VacuumCherenkov::getDefaultSpectrum(const ref_ptr<Kinematics>& kin) {
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
double VacuumCherenkov::thresholdMomentum(const int& id, const ref_ptr<Kinematics>& kinOt, const ref_ptr<Kinematics>& kinPh) {
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
	double eta = - 8 - 6 * sqrt(2);
	double kappa = 0.25;

	double p0 = pow(m2 * Epl2, kappa);
	if (chiOt > 0. and chiPh >= eta * chiOt) {
		pThr = pow(1. / 3. / chiOt, kappa) * p0;
	} else if ((eta * chiOt <= 0 and chiPh < eta * chiOt) or (chiOt <= 0 and chiPh < chiOt)) {
		double t = (chiOt + 2 * chiPh) / (chiOt - chiPh);
		double l = chiOt - chiPh;
		double F = 2. / 27. * l * (t * t * t + pow_integer<3>(sqrt(t * t - 3)) - 4.5 * t);
		pThr = pow(1. / F, kappa) * p0;
	}
	pThr /= c_light;

	return pThr;
}

template<class KO, class KP>
void VacuumCherenkov::buildSpectrum(const KO& kinOt, const KP& kinPh) {
	throw runtime_error("Full treatment of vacuum Cherenkov spectrum for this particular combination of kinematics of photon + other particle is not implemented.");
}

template<>
void VacuumCherenkov::buildSpectrum(const ref_ptr<Kinematics>& kinOt, const ref_ptr<Kinematics>& kinPh) {
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
void VacuumCherenkov::buildSpectrum(const MonochromaticLorentzViolatingKinematics<2>& kinOt, const KP& kinPh) {

	std::pair<double, double> range = xRange(kinOt, kinPh);
	double xMin = std::max(range.first, minEnergyFractionSpectrum);
	double xMax = std::min(range.second, 1.);

	ref_ptr<Histogram1D> dist = new Histogram1DLog10(xMin, xMax, nBinsSpectrum);

	double chiOt = kinOt.getCoefficient();
	double chiPh = kinPh.getCoefficient();

	double g0 = 0;
	if ((chiOt >= 0 and chiPh <= chiOt) or (chiOt == 0 and chiOt > chiPh)) {
		g0 = _G0(chiOt, chiPh);
	} else if (chiOt < 0 and chiPh < chiOt) {
		g0 = (_G0(chiOt, chiPh) - _Gm(chiOt, chiPh));
	} else if (chiOt > 0 and chiOt < chiPh) {
		g0 = _Gp(chiOt, chiPh);
	}

	double P;
	for (size_t i = 0; i < dist->getNumberOfBins(); i++) {
		double x = dist->getBinCentre(i);
		P = 0;
		if (x >= range.first and x < range.second) {
			double w0 = 0.5 * chiOt * (x * x * x - 3 * x * x + 3 * x) - 0.5 * chiPh * x * x * x;
			double P = (2. / x - 2. + x)  * w0 / g0;
			P = std::max(P, 0.);
			dist->setBinContent(i, P);
		}
	}

	dist->setIsPDF(true);
	sampler->setDistribution(dist);
	distribution = dist;

	switch (sampler->getType()) {
		case SamplerType::Inverse: {
			InverseSampler* s = static_cast<InverseSampler*>(sampler.get());
			sampler = new InverseSampler(dist);
			break;
		}
		case SamplerType::Importance: {
			ImportanceSampler* s = static_cast<ImportanceSampler*>(sampler.get());
			std::function<double(double)> wFunc = s->getWeightFunction();
			ref_ptr<Histogram1D> aux = new Histogram1DLog10(xMin, 1., nBins);
			for (size_t i = 0; i < aux->getNumberOfBins(); i++) {
				double x = dist->getBinCentre(i);
				if (x >= range.first and x <= range.second) {
					aux->setBinContent(i, wFunc(aux->getBinCentre(i)));
				}
			}
			sampler = new ImportanceSampler(dist, aux, wFunc);
			break;
		}
		case SamplerType::Nested: {
			NestedSampler* s = static_cast<NestedSampler*>(sampler.get());
			sampler = new NestedSampler(dist, s->getNumberOfLivePoints());
			break;
		}
		default:
			throw runtime_error("VacuumCherenkov: unknown sampler type.");
	}
}

template<class KO, class KP>
double VacuumCherenkov::interactionRate(const double& p, const KO& kinOt, const KP& kinPh) {
	return _defaultInteractionRate;
}

template<>
double VacuumCherenkov::interactionRate(const double& p, const ref_ptr<Kinematics>& kinOt, const ref_ptr<Kinematics>& kinPh) {
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

	if ((chiOt >= 0 and chiPh <= chiOt) or (chiOt == 0 and chiOt > chiPh)) {
		return _G0(chiOt, chiPh) * Q;
	} else if (chiOt < 0 and chiPh < chiOt) {
		return Q * (_G0(chiOt, chiPh) - _Gm(chiOt, chiPh));
	} else if (chiOt > 0 and chiOt < chiPh) {
		return Q * _Gp(chiOt, chiPh);
	}

	return _defaultInteractionRate;;
}

template<>
double VacuumCherenkov::interactionRate(const double& p, const MonochromaticLorentzViolatingKinematics<2>& kinOt, const SpecialRelativisticKinematics& kinPh) {
	const MonochromaticLorentzViolatingKinematics<2> kinPhNew(0.);
	return interactionRate(p, kinOt, kinPhNew);
}

template<class KO, class KP>
std::pair<double, double> VacuumCherenkov::xRange(const KO& kinOt, const KP& kinPh) {
	return std::pair<double, double>(0, 0);
}

template<>
std::pair<double, double> VacuumCherenkov::xRange(const ref_ptr<Kinematics>& kinOt, const ref_ptr<Kinematics>& kinPh) {
	string kOt = kinOt->getNameTag();
	string kPh = kinPh->getNameTag();

	if (kOt == "LIVmono2") {
		const MonochromaticLorentzViolatingKinematics<2>& kinOtNew = kinOt->toMonochromaticLorentzViolatingKinematics2();
		const MonochromaticLorentzViolatingKinematics<2>& kinPhNew = kinPh->toMonochromaticLorentzViolatingKinematics2();
		return xRange(kinOtNew, kinPhNew);
	}

	return std::make_pair(-1, -1.);
}

template<>
std::pair<double, double> VacuumCherenkov::xRange(const MonochromaticLorentzViolatingKinematics<2>& kinOt, const MonochromaticLorentzViolatingKinematics<2>& kinPh)  {
	double chiOt = kinOt.getCoefficient();
	double chiPh = kinPh.getCoefficient();

	double xMin = -1;
	double xMax = -1;

	if (chiOt >= 0 and chiPh <= chiOt) {
		xMin = 0;
		xMax = 1;
	} else if (chiOt > 0 and chiPh > chiOt) {
		xMin = 0;
		xMax = -1.5 * chiOt / (chiPh - chiOt) + 0.5 * sqrt((3. * chiOt * (4. * chiPh - chiOt))) / abs(chiPh - chiOt);
	} else if (chiOt < 0 and chiPh < chiOt) {
		xMin = -1.5 * chiOt / (chiPh - chiOt) + 0.5 * sqrt((3. * chiOt * (4. * chiPh - chiOt))) / abs(chiPh - chiOt);
		xMax = 1;
	}

	return std::make_pair(xMin, xMax);
}


double VacuumCherenkov::_Gp(const double& chiOt, const double& chiPh) {
	double S = sqrt(3. * chiOt * (4 * chiPh - chiOt));
	double Gp0 = chiOt * (S - 3 * chiOt) / 160. / pow_integer<4>(chiPh - chiOt);
	double Gp1 = 37. * (S - 6 * chiPh) * chiPh * chiOt - 64. * pow_integer<3>(chiOt);
	double Gp2 = - pow_integer<2>(chiOt) * (14 * S - 207. * chiPh) - 10 * (5 * S - 16 * chiPh) * pow_integer<2>(chiPh);
	return Gp0 * (Gp1 + Gp2);
}

double VacuumCherenkov::_Gm(const double& chiOt, const double& chiPh) {
	double S = sqrt(3. * chiOt * (4 * chiPh - chiOt));
	double Gm0 = chiOt * (- S - 3 * chiOt) / 160. / pow_integer<4>(chiPh - chiOt);
	double Gm1 = 37. * (- S - 6 * chiPh) * chiPh * chiOt - 64. * pow_integer<3>(chiOt);
	double Gm2 = - pow_integer<2>(chiOt) * (- 14 * S - 207. * chiPh) - 10 * (- 5 * S - 16 * chiPh) * pow_integer<2>(chiPh);
	return Gm0 * (Gm1 + Gm2);
}

double VacuumCherenkov::_G0(const double& chiOt, const double& chiPh) {
	return (157. * chiOt - 22. * chiPh) / 120.;
}


} // namespace livpropa
