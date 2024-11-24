#ifndef LIVPROPA_VACUUMCHERENKOV_H
#define LIVPROPA_VACUUMCHERENKOV_H

#include <cmath>
#include <fstream>
#include <functional>
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
 The default behaviour is:
 . n=0, n=1: step
 . n=2: full spectrum.
 Therefore, whenever a full spectrum treatment is available, it is preferred.
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
 . arXiv:2312.10803 (for n=0 and n=1)
 . arXiv:2410.XXXXX (for n=2)

 The emission is considered to be instantaneous, although in the future a continuous energy loss implementation will be available.
 This can be controlled through the constructor flag `continuousEnergyLoss` which, by default, is false.

 The emission can be either step-like or the full spectrum can be used.
 For step-like spectra, exactly one Cherenkov photon is emitted at a time, with an energy equal to the difference between the primary's energy and the threshold energy.
 Note that several photons can be produced at a time for the full spectrum case, such that this treatment can be time-consuming.
 In this case, the user can control the type of sampling using custom functions. However, this does not work very well outside C++. The default implementation works reasonably.
*/
class VacuumCherenkov: public Module {
	protected:
		string interactionTag;
		int particleId;
		bool havePhotons;
		bool angularCorrection;
		bool continuousEnergyLoss;
		double limit;
		double thinning;
		VacuumCherenkovSpectrum spectrum;
		ref_ptr<Kinematics> kinematicsPhoton;
		ref_ptr<Kinematics> kinematicsParticle;
		ref_ptr<Histogram1D> distribution;
		ref_ptr<DistributionSampler> sampler;
		std::function<double(double)> weightFunction;

	public:
		VacuumCherenkov(int id, KinematicsMap kin, VacuumCherenkovSpectrum spec = VacuumCherenkovSpectrum::Default, bool havePhotons = true, bool angularCorrection = false, bool continuousEnergyLoss = false, ref_ptr<DistributionSampler> sampler = nullptr, double limit = 0.1);
		VacuumCherenkov(int id, ref_ptr<Kinematics> kinOt, ref_ptr<Kinematics> kinPh, VacuumCherenkovSpectrum spec = VacuumCherenkovSpectrum::Default, bool havePhotons = true, bool angularCorrection = false, bool continuousEnergyLoss = false, ref_ptr<DistributionSampler> sampler = nullptr, double limit = 0.1);
		VacuumCherenkov(int id, ref_ptr<Kinematics> kin, VacuumCherenkovSpectrum spec = VacuumCherenkovSpectrum::Default, bool havePhotons = true, bool angularCorrection = false, bool continuousEnergyLoss = false, ref_ptr<DistributionSampler> sampler = nullptr, double limit = 0.1);
		void setParticle(int id);
		void setKinematicsParticle(ref_ptr<Kinematics> kin);
		void setKinematicsPhoton(ref_ptr<Kinematics> kin);
		void setAngularCorrection(bool correction);
		void setContinuousEnergyLoss(bool loss);
		void setHavePhotons(bool photons);
		void setLimit(double limit);
		void setInteractionTag(string tag);
		void setWeightFunction(std::function<double(double)> func);
		void setSpectrum(VacuumCherenkovSpectrum spec, ref_ptr<DistributionSampler> sampler);
		int getParticle() const;
		string getInteractionTag() const;
		std::function<double(double)> getWeightFunction() const;
		ref_ptr<Kinematics> getKinematicsParticle() const;
		ref_ptr<Kinematics> getKinematicsPhoton() const;
		ref_ptr<Histogram1D> getDistribution() const;
		ref_ptr<DistributionSampler> getSampler() const;
		double computeThresholdMomentum() const;
		double computeThresholdEnergy() const;
		double computeInteractionRate(const double& p) const;
		void process(Candidate* candidate) const;
		void emissionSpectrumStep(Candidate* candidate, const double& Ethr) const;
		void emissionSpectrumFull(Candidate* candidate, const double& Ethr) const;
		static VacuumCherenkovSpectrum getDefaultSpectrum(const ref_ptr<Kinematics>& kin);
		template<class KO, class KP> static double thresholdMomentum(const int& id, const KO& kinOt, const KP& kinPh);
		template<> static double thresholdMomentum(const int& id, const ref_ptr<Kinematics>& kinOt, const ref_ptr<Kinematics>& kinPh);
		template<> static double thresholdMomentum(const int& id, const MonochromaticLorentzViolatingKinematics<0>& kinOt, const MonochromaticLorentzViolatingKinematics<0>& kinPh);
		template<> static double thresholdMomentum(const int& id, const MonochromaticLorentzViolatingKinematics<1>& kinOt, const MonochromaticLorentzViolatingKinematics<1>& kinPh);
		template<> static double thresholdMomentum(const int& id, const MonochromaticLorentzViolatingKinematics<2>& kinOt,  const MonochromaticLorentzViolatingKinematics<2>& kinPh);
		template<class KO, class KP> void buildSpectrum(const KO& kinOt, const KP& kinPh);
		template<class KP> void buildSpectrum(const MonochromaticLorentzViolatingKinematics<2>& kinOt, const KP& kinPh);
		template<> void buildSpectrum(const ref_ptr<Kinematics>& kinOt, const ref_ptr<Kinematics>& kinPh);
		template<class KO, class KP> static double interactionRate(const double& p, const KO& kinOt, const KP& kinPh);
		template<> static double interactionRate(const double& p, const ref_ptr<Kinematics>& kinOt, const ref_ptr<Kinematics>& kinPh);
		template<> static double interactionRate(const double& p, const MonochromaticLorentzViolatingKinematics<2>& kinOt, const MonochromaticLorentzViolatingKinematics<2>& kinPh);
		template<> static double interactionRate(const double& p, const MonochromaticLorentzViolatingKinematics<2>& kinOt, const SpecialRelativisticKinematics& kinPh);
		template<class KO, class KP> static std::pair<double, double> xRange(const KO& kinOt, const KP& kinPh); 
		template<> std::pair<double, double> xRange(const ref_ptr<Kinematics>& kinOt, const ref_ptr<Kinematics>& kinPh);
		template<> static std::pair<double, double> xRange(const MonochromaticLorentzViolatingKinematics<2>& kinOt, const MonochromaticLorentzViolatingKinematics<2>& kinPh);

	private:
		static double _Gp(const double& chiOt, const double& chiPh);
		static double _Gm(const double& chiOt, const double& chiPh);
		static double _G0(const double& chiOt, const double& chiPh);
		static ref_ptr<DistributionSampler> _getDefaultSampler();
		template<class KO, class KP> static std::function<double(double)> _getDefaultWeightFunction(const KO& kinOt, const KP& kinPh);
		template<class KP> static std::function<double(double)> _getDefaultWeightFunction(const MonochromaticLorentzViolatingKinematics<2>& kinOt, const KP& kinPh);
		template<> static std::function<double(double)> _getDefaultWeightFunction(const ref_ptr<Kinematics>& kinOt, const ref_ptr<Kinematics>& kinPh);

		static constexpr double _defaultInteractionRate = 0;
		static constexpr double _defaultThresholdMomentum = std::numeric_limits<double>::infinity();

};
/** @}*/

// SWIG sometimes uses weird names; this exposes the right names to the Python interface
constexpr VacuumCherenkovSpectrum VacuumCherenkovSpectrumDefault = VacuumCherenkovSpectrum::Default;
constexpr VacuumCherenkovSpectrum VacuumCherenkovSpectrumStep = VacuumCherenkovSpectrum::Step;
constexpr VacuumCherenkovSpectrum VacuumCherenkovSpectrumFull = VacuumCherenkovSpectrum::Full;



} // namespace livpropa

#endif // LIVPROPA_VACUUMCHERENKOV_H
