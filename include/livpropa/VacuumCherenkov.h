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
 @class EmissionSpectrum
 @brief Enumerate the possible emission spectra for the Cherenkov radiation.
*/
enum class VacuumCherenkovSpectrum {
	Default,
	Step,
	Full,
	Absent
};


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
*/
class VacuumCherenkov: public Module {
	private:
		static constexpr double _defaultInteractionRate = 0;//std::numeric_limits<double>::infinity();
		static constexpr double _defaultThresholdMomentum = std::numeric_limits<double>::infinity();

	protected:
		int particleId;
		string interactionTag;
		bool havePhotons;
		double limit;
		double thinning;
		VacuumCherenkovSpectrum spectrum;
		ref_ptr<AbstractKinematics> kinematicsPhoton;
		ref_ptr<AbstractKinematics> kinematicsParticle;
		ref_ptr<SamplerEvents> samplerEvents;
		ref_ptr<SamplerDistribution> samplerDistribution;
		ref_ptr<Histogram1D> distribution;
		int maximumSamples;
		

	public:
		VacuumCherenkov(int id, Kinematics kin, VacuumCherenkovSpectrum spec = VacuumCherenkovSpectrum::Default, bool havePhotons = true, ref_ptr<SamplerEvents> samplerEvents = nullptr, ref_ptr<SamplerDistribution> samplerDistribution = nullptr, int maximumSamples = 10, double limit = 0.1);
		VacuumCherenkov(int id, ref_ptr<AbstractKinematics> kinOt, ref_ptr<AbstractKinematics> kinPh, VacuumCherenkovSpectrum spec = VacuumCherenkovSpectrum::Default, bool havePhotons = true, ref_ptr<SamplerEvents> samplerEvents = nullptr, ref_ptr<SamplerDistribution> samplerDistribution = nullptr, int maximumSamples = 10, double limit = 0.1);
		VacuumCherenkov(int id, ref_ptr<AbstractKinematics> kin, VacuumCherenkovSpectrum spec = VacuumCherenkovSpectrum::Default, bool havePhotons = true, ref_ptr<SamplerEvents> samplerEvents = nullptr, ref_ptr<SamplerDistribution> samplerDistribution = nullptr, int maximumSamples = 10, double limit = 0.1);
		void setParticle(int id);
		void setKinematicsParticle(ref_ptr<AbstractKinematics> kin);
		void setKinematicsPhoton(ref_ptr<AbstractKinematics> kin);
		void setHavePhotons(bool photons);
		void setLimit(double limit);
		void setInteractionTag(string tag);
		void setSpectrum(VacuumCherenkovSpectrum spec);
		void setSamplerEvents(ref_ptr<SamplerEvents> sampler);
		void setSamplerDistribution(ref_ptr<SamplerDistribution> sampler);
		void setMaximumSamples(int nSamples);
		int getParticle() const;
		string getInteractionTag() const;
		ref_ptr<Histogram1D> getDistribution() const;
		double computeThresholdMomentum() const;
		double computeThresholdEnergy() const;
		double computeInteractionRate(const double& p) const;
		void process(Candidate* candidate) const;
		void emissionSpectrumStep(Candidate* candidate, const double& Ethr) const;
		void emissionSpectrumFull(Candidate* candidate, const double& Ethr) const;
		static VacuumCherenkovSpectrum getDefaultSpectrum(const ref_ptr<AbstractKinematics>& kin);
		template<class KO, class KP> static double thresholdMomentum(const int& id, const KO& kinOt, const KP& kinPh);
		template<> static double thresholdMomentum(const int& id, const ref_ptr<AbstractKinematics>& kinOt, const ref_ptr<AbstractKinematics>& kinPh);
		template<> static double thresholdMomentum(const int& id, const MonochromaticLorentzViolatingKinematics<0>& kinOt, const MonochromaticLorentzViolatingKinematics<0>& kinPh);
		template<> static double thresholdMomentum(const int& id, const MonochromaticLorentzViolatingKinematics<1>& kinOt, const MonochromaticLorentzViolatingKinematics<1>& kinPh);
		template<> static double thresholdMomentum(const int& id, const MonochromaticLorentzViolatingKinematics<2>& kinOt,  const MonochromaticLorentzViolatingKinematics<2>& kinPh);
		template<class KO, class KP> ref_ptr<Histogram1D> buildSpectrum(const KO& kinOt, const KP& kinPh);
		template<> ref_ptr<Histogram1D> buildSpectrum(const ref_ptr<AbstractKinematics>& kinOt, const ref_ptr<AbstractKinematics>& kinPh);
		template<class KP> ref_ptr<Histogram1D> buildSpectrum(const MonochromaticLorentzViolatingKinematics<2>& kinOt, const KP& kinPh);
		template<class KO, class KP> static double interactionRate(const double& p, const KO& kinOt, const KP& kinPh);
		template<> static double interactionRate(const double& p, const ref_ptr<AbstractKinematics>& kinOt, const ref_ptr<AbstractKinematics>& kinPh);
		template<> static double interactionRate(const double& p, const MonochromaticLorentzViolatingKinematics<2>& kinOt, const MonochromaticLorentzViolatingKinematics<2>& kinPh);
		template<> static double interactionRate(const double& p, const MonochromaticLorentzViolatingKinematics<2>& kinOt, const SpecialRelativisticKinematics& kinPh);
		// // template<class KO, class KP> double differentialInteractionRate(const double& p, const KO& kinOt, const KP& kinPh);
		// // template<> double differentialInteractionRate(const double& p, const MonochromaticLorentzViolatingKinematics<2>& kinOt, const SpecialRelativisticKinematics& kinPh);
		template<class KO, class KP> static std::pair<double, double> xRange(const KO& kinOt, const KP& kinPh); 
		template<> std::pair<double, double> xRange(const ref_ptr<AbstractKinematics>& kinOt, const ref_ptr<AbstractKinematics>& kinPh);
		template<> static std::pair<double, double> xRange(const MonochromaticLorentzViolatingKinematics<2>& kinOt, const MonochromaticLorentzViolatingKinematics<2>& kinPh);
		// template<class KO, class KP> static double differentialProbability(const double& x, const KO& kinOt, const KP& kinPh);
		// template<> static double differentialProbability(const double& x, const MonochromaticLorentzViolatingKinematics<2>& kinOt, const MonochromaticLorentzViolatingKinematics<2>& kinPh);
		// template<> static double differentialProbability(const double& x, const MonochromaticLorentzViolatingKinematics<2>& kinOt, const SpecialRelativisticKinematics& kinPh);
};
/** @}*/

// SWIG sometimes uses weird names; this exposes the right names to the Python interface
constexpr VacuumCherenkovSpectrum VacuumCherenkovSpectrumDefault = VacuumCherenkovSpectrum::Default;
constexpr VacuumCherenkovSpectrum VacuumCherenkovSpectrumStep = VacuumCherenkovSpectrum::Step;
constexpr VacuumCherenkovSpectrum VacuumCherenkovSpectrumFull = VacuumCherenkovSpectrum::Full;




// namespace vc::monoLIV0 {
// 	double computeThresholdMomentum(const double& chiOt, const double& chiPh, const double& mass);
// } // namespace vc::monoLIV0

// namespace vc::monoLIV1 {
// 	double computeThresholdMomentum(const double& chiOt, const double& chiPh, const double& mass);
// } // namespace vc::monoLIV1

// namespace vc::monoLIV2 {

// 	double _S(const double& chiOt, const double& chiPh);

// 	double _tau(const double& chiOt, const double& chiPh); 

// 	double _F(const double& chiOt, const double& chiPh);

// 	double _G(const double& chiOt, const double& chiPh, const double& sign);

// 	double _x(const double& chiOt, const double& chiPh, const double& sign);

// 	double _rateQ(const double& p);

// 	double _rate1(const double& p, const double& chiOt, const double& chiPh);

// 	double _rate2(const double& p, const double& chiOt, const double& chiPh);

// 	double _rate3(const double& p, const double& chiOt, const double& chiPh);

// 	double computeSpectrum(const double& x, const double& chiOt, const double& chiPh);

// 	void loadSpectrumDistribution(Histogram1D& distribution, const double& chiOt, const double& chiPh); 

// 	double computeThresholdMomentum(const double& chiOt, const double& chiPh, const double& mass);

// 	double computeInteractionRate(const double& p, const double& chiOt, const double& chiPh);

// } // namespace vc::monoLIV2


} // namespace livpropa

#endif // LIVPROPA_VACUUMCHERENKOV_H
