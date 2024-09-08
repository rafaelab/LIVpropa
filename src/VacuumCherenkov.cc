#include "livpropa/VacuumCherenkov.h"


namespace livpropa {



VacuumCherenkov::VacuumCherenkov(Kinematics kinematics, EmissionSpectraTable spec,  bool havePhotons, ref_ptr<SamplerEvents> samplerEvents, ref_ptr<SamplerDistribution> samplerDistribution, int maxSamples, double limit) {
	setInteractionTag("VC");
	setHavePhotons(havePhotons);
	setLimit(limit);
	setKinematics(kinematics);
	setSpectra(spec);
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

void VacuumCherenkov::setSpectra(EmissionSpectraTable spec, const unsigned int nPoints) {
	if (spec.empty()) {
		spectra[ 11] = EmissionSpectrum::Default;
		spectra[-11] = EmissionSpectrum::Default;
	}


	// vector<int> particles = kinematics.getParticles();
	// if (std::find(particles.begin(), particles.end(), x) != v.end())


	// auto kinPh = kinematics[ 22];


	// std::string kin = kinematics.getNameTag();

	// if (kin == "SR") {
	// 	if (spec != EmissionSpectrum::Default)
	// 		throw std::runtime_error("No vacuum Cherenkov implemented for special-relativistic kinematics.");

	// } else if (kin == "LIV") {
	// 	MonochromaticLIV* kin = static_cast<MonochromaticLIV*>(kinematics.get()); 
	// 	unsigned int order = kin->getOrder();

	// 	if (order == 0 || order == 1) {
	// 		if (spec == EmissionSpectrum::Default)
	// 			spec = EmissionSpectrum::Step;

	// 		if (spec == EmissionSpectrum::Step) {
	// 			return;
	// 		} else {
	// 			throw std::runtime_error("Only step-like spectrum is implemented for LIV of order 0 and 1.");
	// 		}
	// 	}

	// 	else if (order == 2) {
	// 		if (spec == EmissionSpectrum::Default)
	// 			spec = EmissionSpectrum::Full;

	// 		if (spec == EmissionSpectrum::Full) {
	// 			vc::monoLIV2::loadSpectrumDistribution(distribution, kin->getCoefficientForParticle(11), kin->getCoefficientForParticle(22));
	// 			return;
	// 		} else {
	// 			throw std::runtime_error("Only full spectrum is implemented for LIV of order 2.");
	// 		}
	// 	}

	// 	throw std::runtime_error("VacuumCherenkov: only LIV of orders 0 and 1 (step-like spectrum) and 2 (full spectrum) are implemented.");
	// }	

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

bool VacuumCherenkov::isImplemented(const int& id) const {
	std::string kinOt = kinematics[id]->getNameTag();
	std::string kinPh = kinematics[22]->getNameTag();

	if (kinOt == "LIVmono" && kinPh == "LIVmono") {
		auto kinOt = kinematics[id]->toLorentzViolatingKinematicsMonochromatic();
		auto kinPh = kinematics[22]->toLorentzViolatingKinematicsMonochromatic();
		if (kinOt->getOrder() == kinPh->getOrder())
			return true;

	} else if (kinOt == "LIVmono" && kinPh == "SR") {
		return true;
	} 
	
	return false;
}

double VacuumCherenkov::computeThresholdMomentum(const int& id) const {
	double pThr = std::numeric_limits<double>::max();
	if (! isImplemented(id))
		return pThr;

	if (kinematics[id]->isLorentzViolatingMonochromatic()) {
		auto kinOt = kinematics[id]->toLorentzViolatingKinematicsMonochromatic();
		double chiOt = kinOt->getCoefficient();

		// note that the SR case can be seen as LIV with coefficient equal to 0
		double chiPh = 0.;
		if (kinematics[22]->isLorentzViolatingMonochromatic()) {
			auto kinPh = kinematics[22]->toLorentzViolatingKinematicsMonochromatic();
			chiPh = kinPh->getCoefficient();
		}

		double m = particleMasses.at(id);
		pThr = vc::monoLIV0::computeThresholdMomentum(chiOt, chiPh, m);

		switch (kinOt->getOrder()) {
			case 0:
				pThr = vc::monoLIV0::computeThresholdMomentum(chiOt, chiPh, m);
				break;
			case 1:
				pThr = vc::monoLIV1::computeThresholdMomentum(chiOt, chiPh, m);
				break;
			case 2:
				pThr = vc::monoLIV2::computeThresholdMomentum(chiOt, chiPh, m);
				break;
			default:	
				throw std::runtime_error("VacuumCherenkov: only LIV of orders 0, 1, and 2 are implemented.");
		}
	}

	return pThr;
}

double VacuumCherenkov::computeThresholdEnergy(const int& id) const {
	double p = computeThresholdMomentum(id);
	return kinematics[id]->computeEnergyFromMomentum(p, id);
}

void VacuumCherenkov::process(Candidate* candidate) const {
	// check if electron/positron (only particles implemented so far)
	int id = candidate->current.getId();
	if (abs(id) != 11)
		return;

	// do not perform any calculations for the special-relativistic case
	std::string kin = kinematics[id]->getNameTag();
	if (kin == "SR")
		return;

	// get the emission spectrum
	EmissionSpectrum spectrum;
	try {
		spectrum = spectra.at(id);
	} catch (const int& id) {
		throw runtime_error("No emission spectrum defined for particle with id " + std::to_string(id) + ". Using default behaviour.\n");
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

			if (spectrum == EmissionSpectrum::Step) {
				double Ephoton = E - Ethr;	
				candidate->addSecondary(22, Ephoton / (1 + z), pos);

			} else if (spectrum == EmissionSpectrum::Full) {
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


///////////////////////////////////////////////////////////////////////////////////////////////////

namespace vc::monoLIV0 {

	double computeThresholdMomentum(const double& chiOt, const double& chiPh, const double& mass) {
		double pThr = std::numeric_limits<double>::max();

		if (chiOt > chiPh) 
			pThr = mass * c_light / sqrt(chiOt - chiPh);

		return pThr;
	}

}

namespace vc::monoLIV1 {

	double computeThresholdMomentum(const double& chiOt, const double& chiPh, const double& mass) {
		double pThr = std::numeric_limits<double>::max();
		double m2 = pow_integer<2>(mass * c_squared);

		if (chiOt > 0. and chiPh >= -3. * chiOt) {
			pThr = cbrt(m2 * energy_planck / (2 * chiOt)) / c_light;
		} else if (chiOt <= 0 and chiOt > chiPh and chiPh <= -3. * chiOt) {
			pThr = cbrt(- 4. * m2 * energy_planck * (chiOt + chiPh) / pow_integer<2>(chiOt - chiPh)) / c_light;
		}

		return pThr;
	}

}

namespace vc::monoLIV2 {

	double _S(const double& chiOt, const double& chiPh) {
		return sqrt(3. * chiOt * (4 * chiPh - chiOt));
	}

	double _tau(const double& chiOt, const double& chiPh) {
		return (chiOt + 2 * chiPh) / (chiPh - chiOt);
	}

	double _F(const double& chiOt, const double& chiPh) {
		double t = _tau(chiOt, chiPh);
		double l = chiOt - chiPh;
		return 2. / 27. * l * (t * t * t + pow_integer<3>(sqrt(t * t - 3)) - 4.5 * t);
	}

	double _G(const double& chiOt, const double& chiPh, const double& sign) {
		double chiOt2 = chiOt * chiOt;
		double chiOt3 = chiOt2 * chiOt;
		double chiPh2 = chiPh * chiPh;
		double dChi = chiPh - chiOt;
		double S = sign * _S(chiOt, chiPh);

		double f0 = chiOt * (S - 3 * chiOt) / 160. / pow_integer<4>(dChi);
		double f1 = 37. * (S - 6 * chiPh) * chiPh * chiOt - 64. * chiOt3;
		double f2 = - chiOt2 * (14 * S - 207. * chiPh) - 10 * (5 * S - 16 * chiPh) * chiPh2;

		return f0 * (f1 + f2);
	}

	double _x(const double& chiOt, const double& chiPh, const double& sign) {
		if (chiOt == chiPh && chiOt != 0) {
			return 1;
		} else {
			double dChi = chiPh - chiOt;
			double f1 = - 1.5 * chiOt / dChi;
			double f2 = sqrt(3 * chiOt * (4 * chiPh - chiOt)) / abs(dChi);
			return f1 + sign * f2;
		}
	}

	double _rateQ(const double& p) {
		double q = pow_integer<3>(p * c_light) / pow_integer<2>(energy_planck);
		return alpha_finestructure * q;
	}

	double _rate1(const double& p, const double& chiOt, const double& chiPh) {
		return (157. * chiOt - 22. * chiPh) / 120. * _rateQ(p);
	}

	double _rate2(const double& p, const double& chiOt, const double& chiPh) {
		return _rateQ(p) * _G(chiOt, chiPh, +1.);
	}

	double _rate3(const double& p, const double& chiOt, const double& chiPh) {
		return _rateQ(p) * ((157. * chiOt - 22. * chiPh) / 120. - _G(chiOt, chiPh, -1.));
	}

	double computeSpectrum(const double& x, const double& chiOt, const double& chiPh) {
		double S = sqrt(3. * chiOt * (4 * chiPh - chiOt));
		double Gp = _G(chiOt, chiPh, +1.);
		double Gm = _G(chiOt, chiPh, -1.);

		double x2 = x * x;
		double x3 = x2 * x;

		double y = x3 - 3. * x2 + 3. * x;
		double a =  (2. / x - 2. + x);
		double b1 = - chiPh * x3 / 2.;
		double b2 = chiOt * y / 2.;
		double b = b1 + b2;
		double c = (157. * chiOt - 22. * chiPh) / 120.;

		if (chiOt >= std::max(0., chiPh)) 
			return a * b / c;
		else if (chiOt > 0 &&  chiOt < chiPh)
			return std::max(0., a * b / (c - Gm));
		else if (chiOt < chiPh && chiPh < 0.)
			return std::max(0., a * b / Gp);
		
		return 0.;
	}

	void loadSpectrumDistribution(Histogram1D& distribution, const double& chiOt, const double& chiPh) {
		for (size_t i = 0; i < distribution.getNumberOfBins(); i++) {
			double x = distribution.getBinCentre(i);
			double p = computeSpectrum(x, chiOt, chiPh);
			distribution.setBinContent(i, p);
		}

		distribution.normalise(distribution.integrate());
		distribution.transformToPDF();
		distribution.transformToCDF();
	}

	double computeThresholdMomentum(const double& chiOt, const double& chiPh, const double& mass) {
		double pThr = std::numeric_limits<double>::max();
		double a = - 8 - 6 * sqrt(2);
		double m2 = pow_integer<2>(mass * c_squared);
		double Epl2 = pow_integer<2>(energy_planck);

		if (chiOt > 0. && chiPh >= a * chiOt) {
			double chiOt2 = chiOt * chiOt;
			pThr = pow(m2 * Epl2 / 3. / chiOt2, 0.25);
		} else if ((a * chiOt <= 0 && chiPh < a * chiOt) || (chiOt <= 0 && chiPh < chiOt)) {
			double F = _F(chiOt, chiPh);
			pThr = pow(pow_integer<2>(m2 * energy_planck) / F, 0.25);
		}

		return pThr;
	}

	double computeInteractionRate(const double& p, const double& chiOt, const double& chiPh) {
		if (chiOt == 0) {
			if (chiOt > chiPh)
				return _rate1(p, chiOt, chiPh);
		} else if (chiOt > 0) {
			return (chiOt < chiPh) ? _rate2(p, chiOt, chiPh) : _rate1(p, chiOt, chiPh);
		} else {
			if (chiOt > chiPh)
				return _rate3(p, chiOt, chiPh);
		}
		

		return 0;
	}

} // namespace vc::monoLIV2



} // namespace livpropa
