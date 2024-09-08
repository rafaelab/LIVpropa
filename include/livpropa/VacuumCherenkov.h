#ifndef LIVPROPA_VACUUMCHERENKOV_H
#define LIVPROPA_VACUUMCHERENKOV_H

#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

#include "livpropa/Common.h"
#include "livpropa/Data.h"
#include "livpropa/Histogram.h"
#include "livpropa/Kinematics.h"
#include "livpropa/Sampler.h"
#include "livpropa/UnitsAndConstants.h"



namespace livpropa {
	
/**
 @class VacuumCherenkov
 @brief Simulate the emission of vacuum Cherenkov radiation due to LIV effects.

 This implementation follows the description from:
 	arXiv:2312.10803 (for n=0 and n=1)
	arXiv:2410.XXXXX (for n=2)
 Exactly one Cherenkov photon is emitted at a time, with an energy equal to the difference between the primary's energy and the threshold energy.
 The emission is considered to be instantaneous (this might change in the future).
 While the structure here is meant to be generic, it was tested only for electrons.
 Note that the treatment of the spectrum is fixed according to the kinematics.
 Since there are no other options, it is treated as step-like for n=0 and n=1, and as full for n=2.

 Three emission spectra are implemented:Default: use default behaviour (step for n=0,1, full for n=2)
*/
class VacuumCherenkov: public Module {
	public:
		enum EmissionSpectrum {
			Default,
			Step,
			Full
		};
		bool isSpectrumDefault(EmissionSpectrum spec) const {
			return spec == EmissionSpectrum::Default	;
		}
		bool isSpectrumStep(EmissionSpectrum spec) const {
			return spec == EmissionSpectrum::Step;
		}
		bool isSpectrumFull(EmissionSpectrum spec) const {
			return spec == EmissionSpectrum::Full;
		}


	protected:
		string interactionTag;
		ref_ptr<AbstractKinematics> kinematics;
		EmissionSpectrum spectrum;
		bool havePhotons;
		double limit;
		double thinning;
		ref_ptr<SamplerEvents> samplerEvents;
		ref_ptr<SamplerDistribution> samplerDistribution;
		int maximumSamples;
		Histogram1D distribution;

	public:
		VacuumCherenkov(ref_ptr<AbstractKinematics> kinematics, EmissionSpectrum spectrum = EmissionSpectrum::Default, bool havePhotons = false, ref_ptr<SamplerEvents> samplerEvents = NULL, ref_ptr<SamplerDistribution> samplerDistribution = NULL, int maximumSamples = 0, double limit = 0.1);
		void setKinematics(ref_ptr<AbstractKinematics> kinematics);
		void setHavePhotons(bool photons);
		void setLimit(double limit);
		void setInteractionTag(string tag);
		void setSpectrum(EmissionSpectrum spectrum, unsigned int nPoints = 1000);
		void setSamplerEvents(ref_ptr<SamplerEvents> sampler);
		void setSamplerDistribution(ref_ptr<SamplerDistribution> sampler);
		void setMaximumSamples(int nSamples);
		string getInteractionTag() const;
		double computeThresholdMomentum(const int& id) const;
		double computeThresholdEnergy(const int& id) const;
		void process(Candidate* candidate) const;

};
/** @}*/



namespace vc::monoLIV0 {
	double computeThresholdMomentum(const double& chiOt, const double& chiPh, const double& mass);
} // namespace vc::monoLIV0

namespace vc::monoLIV1 {
	double computeThresholdMomentum(const double& chiOt, const double& chiPh, const double& mass);
} // namespace vc::monoLIV1

namespace vc::monoLIV2 {

	double _S(const double& chiOt, const double& chiPh);

	double _tau(const double& chiOt, const double& chiPh); 

	double _F(const double& chiOt, const double& chiPh);

	double _G(const double& chiOt, const double& chiPh, const double& sign);

	double _x(const double& chiOt, const double& chiPh, const double& sign);

	double _rateQ(const double& p);

	double _rate1(const double& p, const double& chiOt, const double& chiPh);

	double _rate2(const double& p, const double& chiOt, const double& chiPh);

	double _rate3(const double& p, const double& chiOt, const double& chiPh);

	double computeSpectrum(const double& x, const double& chiOt, const double& chiPh);

	void loadSpectrumDistribution(Histogram1D& distribution, const double& chiOt, const double& chiPh); 

	double computeThresholdMomentum(const double& chiOt, const double& chiPh, const double& mass);

	double computeInteractionRate(const double& p, const double& chiOt, const double& chiPh);

} // namespace vc::monoLIV2


} // namespace livpropa

#endif // LIVPROPA_VACUUMCHERENKOV_H
