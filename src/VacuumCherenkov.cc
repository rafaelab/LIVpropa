#include "livpropa/VacuumCherenkov.h"


namespace livpropa {



VacuumCherenkov::VacuumCherenkov(Kinematics kinematics, bool havePhotons, ref_ptr<SamplerEvents> samplerEvents, ref_ptr<SamplerDistribution> samplerDistribution, int maxSamples, double limit) {
	setInteractionTag("VC");
	setHavePhotons(havePhotons);
	setLimit(limit);
	setKinematics(kinematics);
	setSpectra(spectra);
	setSamplerEvents(samplerEvents);
	setSamplerDistribution(samplerDistribution);
	setMaximumSamples(maxSamples);
}

void VacuumCherenkov::setKinematics(Kinematics kin) {
	kinematics = kin;
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

void VacuumCherenkov::setSpectrumTypeForParticle(int id, VacuumCherenkovSpectrum spec) {
	if (! isTreatmentImplemented(kinematics[id].get(), spec))
		throw runtime_error("VacuumCherenkov treatment for this particular combination of kinematics of photon + other particle is not implemented.");

	if (spec == VacuumCherenkovSpectrum::Default)
		spectra[id] = getDefaultSpectrum(kinematics[id]);
	else
		spectra[id] = spec;
}

void VacuumCherenkov::setSpectra(unordered_map<int, VacuumCherenkovSpectrum> emissionSpectra, const unsigned int nPoints) {
	if (emissionSpectra.empty()) {
		for (const auto& s : emissionSpectra) {
			int id = s.first;
			if (id != 22)
				emissionSpectra[id] = VacuumCherenkovSpectrum::Default;
		}
	} 

	for (const auto& s : emissionSpectra) {
		int id = s.first;
		VacuumCherenkovSpectrum spec = s.second;
		if (id != 22)
			setSpectrumTypeForParticle(id, spec);		
	}

	return;
}

void VacuumCherenkov::setSamplerEvents(ref_ptr<SamplerEvents> s) {
	if (s == NULL)
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

std::string VacuumCherenkov::getInteractionTag() const {
	return interactionTag;
}

double VacuumCherenkov::computeThresholdMomentum(const int& id) const {
	auto kinOt = kinematics[id];
	auto kinPh = kinematics[22];
	return thresholdMomentum(id, kinOt, kinPh);
}

double VacuumCherenkov::computeThresholdEnergy(const int& id) const {
	auto kinOt = kinematics[id];
	auto kinPh = kinematics[22];
	double p = thresholdMomentum(id, kinOt, kinPh);
	return kinematics[id]->computeEnergyFromMomentum(p, id);
}

double VacuumCherenkov::computeInteractionRate(const int& id, const double& p) const {
	auto kinOt = kinematics[id];
	auto kinPh = kinematics[22];
	return interactionRate(p, kinOt, kinPh);
}

void VacuumCherenkov::process(Candidate* candidate) const {
	// check if electron/positron (only particles implemented so far)
	int id = candidate->current.getId();
	if (! kinematics.exists(id))
		return;

	// do not perform any calculations for the special-relativistic case
	std::string k = kinematics[id]->getNameTag();
	if (k == "SR")
		return;
		
	else if (! kinematics[id]->isLorentzViolatingMonochromatic())
		throw runtime_error("VacuumCherenkov treatment for this particular combination of kinematics of photon + other particle is not implemented.");
	
	VacuumCherenkovSpectrum spectrum = VacuumCherenkovSpectrum::Default;
	if (spectra.find(id) == spectra.end()) {
		spectrum = getDefaultSpectrum(kinematics[id]);
	} else {
		spectrum = spectra.at(id);
	}

	// get and scale the particle energy
	double z = candidate->getRedshift();
	double E = candidate->current.getEnergy() * (1 + z);

	// compute threshold energy
	double Ethr = computeThresholdEnergy(id) * (1 + z);


	// photon emission if above threshold
	if (E > Ethr) {
		candidate->current.setEnergy(Ethr / (1 + z));
		
		if (havePhotons) {
			// position where photon will be produced
			Vector3d pos = candidate->current.getPosition();

			if (spectrum == VacuumCherenkovSpectrum::Step) {
				double Ephoton = E - Ethr;	
				candidate->addSecondary(22, Ephoton / (1 + z), pos);

			} else if (spectrum == VacuumCherenkovSpectrum::Full) {
				double dE0 = E - Ethr;	// energy available for photons
				if (samplerDistribution == NULL) {
					double dE = dE0;
					while (dE > 0) {
						double Ephoton = distribution.getSample() * E;
						candidate->addSecondary(22, Ephoton / (1 + z), pos);
						dE -= Ephoton;
					}
				} else {
					samplerDistribution->transformToPDF();
					samplerDistribution->transformToCDF();
					std::vector<double> sampled = samplerDistribution->getSample(maximumSamples);
					double dEs = std::accumulate(sampled.begin(), sampled.end(), decltype(sampled)::value_type(0));
					if (samplerDistribution->getSize() > 0) {
						double w0 = dE0 / dEs;
						for (size_t i = 0; i < sampled.size(); i++) {
							double Es = sampled[i];
							double w = w0 * samplerEvents->computeWeight(-11, Es, Es / E, i);
							candidate->addSecondary(22, Es, pos, w);
						}
					}
					samplerDistribution->clear();
				}
			}

		}
	}

	return;
}

bool VacuumCherenkov::isTreatmentImplemented(const ref_ptr<AbstractKinematics>& kin, VacuumCherenkovSpectrum spec) {
	string k = kin.get()->getNameTag();
	if (k == "SR") {
		return true;
	} else if (k == "LIVmono0") {
		return (spec == VacuumCherenkovSpectrum::Step) ? true : false;
	} else if (k == "LIVmono1") {
		return (spec == VacuumCherenkovSpectrum::Step) ? true : false;
	} else if (k == "LIVmono2") {
		return (spec == VacuumCherenkovSpectrum::Step or spec == VacuumCherenkovSpectrum::Full) ? true : false;
	} 
	throw runtime_error("VacuumCherenkov: treatment for this particular combination of kinematics of photon + other particle is not implemented.");	
}

VacuumCherenkovSpectrum VacuumCherenkov::getDefaultSpectrum(const ref_ptr<AbstractKinematics>& kin) {
	string kinType = kin.get()->getNameTag();
	
	if (kinType == "SR") {
		return VacuumCherenkovSpectrum::Absent;
	} else if (kinType == "LIVmono0") {
		return VacuumCherenkovSpectrum::Step;
	} else if (kinType == "LIVmono1") {
		return VacuumCherenkovSpectrum::Step;
	} else if (kinType == "LIVmono2") {
		return VacuumCherenkovSpectrum::Full;
	} 
	
	throw runtime_error("VacuumCherenkov: default spectrum can be retrieved only for monochromatic LIV with n=0,1, and 2.");	
}

template<class KO, class KP>
double VacuumCherenkov::thresholdMomentum(const int& id, const KO& kinOt, const KP& kinPh) {
	// throw runtime_error("VacuumCherenkov treatment for this particular combination of kinematics of photon + other particle is not implemented.");
	return std::numeric_limits<double>::max();
}

template<>
double VacuumCherenkov::thresholdMomentum(const int& id, const ref_ptr<AbstractKinematics>& kinOt, const ref_ptr<AbstractKinematics>& kinPh) {
	double pThr = std::numeric_limits<double>::max();

	string kOt = kinOt->getNameTag();
	string kPh = kinPh->getNameTag();

	if (kOt == "LIVmono0") {
		const auto& kinOtNew = kinOt->toMonochromaticLorentzViolatingKinematics0();
		if (kPh == "SR") {
			const auto& kinPhNew = kinPh->toSpecialRelativisticKinematics();
			return thresholdMomentum(id, kinOtNew, kinPhNew);		
		} else if (kPh == "LIVmono0") {
			const auto& kinPhNew = kinPh->toMonochromaticLorentzViolatingKinematics0();
			return thresholdMomentum(id, kinOtNew, kinPhNew);	
		}
	} else if (kOt == "LIVmono1") {
		const auto& kinOtNew = kinOt->toMonochromaticLorentzViolatingKinematics1();
		if (kPh == "SR") {
			const auto& kinPhNew = kinPh->toSpecialRelativisticKinematics();
			return thresholdMomentum(id, kinOtNew, kinPhNew);		
		} else if (kPh == "LIVmono1") {
			const auto& kinPhNew = kinPh->toMonochromaticLorentzViolatingKinematics1();
			return thresholdMomentum(id, kinOtNew, kinPhNew);	
		}
	} else if (kOt == "LIVmono2") {
		const auto& kinOtNew = kinOt->toMonochromaticLorentzViolatingKinematics2();
		if (kPh == "SR") {
			const auto& kinPhNew = kinPh->toSpecialRelativisticKinematics();
			return thresholdMomentum(id, kinOtNew, kinPhNew);		
		} else if (kPh == "LIVmono2") {
			const auto& kinPhNew = kinPh->toMonochromaticLorentzViolatingKinematics2();
			return thresholdMomentum(id, kinOtNew, kinPhNew);	
		}
	}

	return pThr;
}

template<>
double VacuumCherenkov::thresholdMomentum(const int& id, const MonochromaticLorentzViolatingKinematics<0>& kinOt, const MonochromaticLorentzViolatingKinematics<0>& kinPh) {
	double pThr = std::numeric_limits<double>::max();

	double mass = particleMasses.at(id);

	double chiOt = kinOt.getCoefficient();
	double chiPh = kinPh.getCoefficient();

	if (chiOt > chiPh) 
		pThr = mass * c_light / sqrt(chiOt - chiPh);

	return pThr;
}

template<>
double VacuumCherenkov::thresholdMomentum(const int& id, const MonochromaticLorentzViolatingKinematics<1>& kinOt, const MonochromaticLorentzViolatingKinematics<1>& kinPh) {
	double pThr = std::numeric_limits<double>::max();

	double chiOt = kinOt.getCoefficient();
	double chiPh = kinPh.getCoefficient();

	double mass = particleMasses.at(id);
	double m2 = pow_integer<2>(mass * c_squared);

	if (chiOt > 0. and chiPh >= -3. * chiOt) {
		pThr = cbrt(m2 * energy_planck / (2 * chiOt)) / c_light;
	} else if (chiOt <= 0 and chiOt > chiPh and chiPh <= -3. * chiOt) {
		pThr = cbrt(- 4. * m2 * energy_planck * (chiOt + chiPh) / pow_integer<2>(chiOt - chiPh)) / c_light;
	}

	return pThr;
}

template<>
double VacuumCherenkov::thresholdMomentum(const int& id, const MonochromaticLorentzViolatingKinematics<2>& kinOt,  const MonochromaticLorentzViolatingKinematics<2>& kinPh) {
	double pThr = std::numeric_limits<double>::max();

	double chiOt = kinOt.getCoefficient();
	double chiPh = kinPh.getCoefficient();

	double mass = particleMasses.at(id);
	double m2 = pow_integer<2>(mass * c_squared);
	double Epl2 = pow_integer<2>(energy_planck);
	double a = - 8 - 6 * sqrt(2);

	if (chiOt > 0. && chiPh >= a * chiOt) {
		double chiOt2 = chiOt * chiOt;
		pThr = pow(m2 * Epl2 / 3. / chiOt2, 0.25);
	} else if ((a * chiOt <= 0 && chiPh < a * chiOt) || (chiOt <= 0 && chiPh < chiOt)) {
		double t = (chiOt + 2 * chiPh) / (chiPh - chiOt);
		double l = chiOt - chiPh;
		double F = 2. / 27. * l * (t * t * t + pow_integer<3>(sqrt(t * t - 3)) - 4.5 * t);
		pThr = pow(pow_integer<2>(m2 * energy_planck) / F, 0.25);
	}

	return pThr;
}

template<class KO, class KP>
Histogram1D VacuumCherenkov::buildSpectrum(const int& id, const KO& kinOt, const KP& kinPh) {
	throw runtime_error("Full treatment of vacuum Cherenkov spectrum for this particular combination of kinematics of photon + other particle is not implemented.");
}

template<class KP>
Histogram1D VacuumCherenkov::buildSpectrum(const int& id, const MonochromaticLorentzViolatingKinematics<2>& kinOt, const KP& kinPh) {
	Histogram1D distribution(100, 0., 1.);

	double chiOt = kinOt.getCoefficient();
	double chiPh = kinPh.getCoefficient();

	double chiOt2 = chiOt * chiOt;
	double chiOt3 = chiOt2 * chiOt;
	double chiPh2 = chiPh * chiPh;
	double dChi = chiPh - chiOt;
	double S = sqrt(3. * chiOt * (4 * chiPh - chiOt));

	double Gp0 = chiOt * (S - 3 * chiOt) / 160. / pow_integer<4>(dChi);
	double Gp1 = 37. * (S - 6 * chiPh) * chiPh * chiOt - 64. * chiOt3;
	double Gp2 = - chiOt2 * (14 * S - 207. * chiPh) - 10 * (5 * S - 16 * chiPh) * chiPh2;
	double Gp = Gp0 * (Gp1 + Gp2);

	double Gm0 = chiOt * (- S - 3 * chiOt) / 160. / pow_integer<4>(dChi);
	double Gm1 = 37. * (- S - 6 * chiPh) * chiPh * chiOt - 64. * chiOt3;
	double Gm2 = - chiOt2 * (- 14 * S - 207. * chiPh) - 10 * (- 5 * S - 16 * chiPh) * chiPh2;
	double Gm = Gm0 * (Gm1 + Gm2);

	double P = 0;
	for (size_t i = 0; i < distribution.getNumberOfBins(); i++) {
		double x = distribution.getBinCentre(i);
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
		else if (chiOt > 0 &&  chiOt < chiPh)
			P = std::max(0., a * b / (c - Gm));
		else if (chiOt < chiPh && chiPh < 0.)
			P = std::max(0., a * b / Gp);
		distribution.setBinContent(i, P);
	}

	distribution.normalise(distribution.integrate());
	distribution.transformToPDF();
	distribution.transformToCDF();

	return distribution;
}

template<class KO, class KP>
double VacuumCherenkov::interactionRate(const double& p, const KO& kinOt, const KP& kinPh) {
	return 0;
}

template<>
double VacuumCherenkov::interactionRate(const double& p, const ref_ptr<AbstractKinematics>& kinOt, const ref_ptr<AbstractKinematics>& kinPh) {
	string kOt = kinOt.get()->getNameTag();
	string kPh = kinPh.get()->getNameTag();

	if (kOt == "SR")
		return 0;

	if (kOt == "LIVmono2") {
		auto kinOtNew = kinOt->toMonochromaticLorentzViolatingKinematics2();
		return interactionRate(p, kinOtNew, kinPh);
	}

	return 0;
}

template<>
double VacuumCherenkov::interactionRate(const double& p, const MonochromaticLorentzViolatingKinematics<2>& kinOt, const MonochromaticLorentzViolatingKinematics<2>& kinPh) {
	double chiOt = kinOt.getCoefficient();
	double chiPh = kinPh.getCoefficient();

	double q = pow_integer<3>(p * c_light) / pow_integer<2>(energy_planck);
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
			double Gp2 = - pow_integer<2>(chiPh) * (14 * S - 207. * chiPh) - 10 * (5 * S - 16 * chiPh) * pow_integer<2>(chiPh);
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

	return 0;
}

template<>
double VacuumCherenkov::interactionRate(const double& p, const MonochromaticLorentzViolatingKinematics<2>& kinOt, const SpecialRelativisticKinematics& kinPh) {
	const MonochromaticLorentzViolatingKinematics<2> kinPhNew(0.);
	return interactionRate(p, kinOt, kinPhNew);
}


} // namespace livpropa
