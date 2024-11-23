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
		static constexpr double _defaultInteractionRate = 0;
		static constexpr double _defaultThresholdMomentum = std::numeric_limits<double>::infinity();
		// ref_ptr<Sampler> _defaultSampler;

	protected:
		string interactionTag;
		int particleId;
		bool havePhotons;
		bool angularCorrection;
		bool continuousEnergyLoss;
		double limit;
		double thinning;
		VacuumCherenkovSpectrum spectrum;
		KinematicsMap kinematics;
		ref_ptr<Histogram1D> distribution;
		ref_ptr<Sampler> sampler;
		

	public:
		VacuumCherenkov(int id, KinematicsMap kin, VacuumCherenkovSpectrum spec = VacuumCherenkovSpectrum::Default, bool havePhotons = true, bool angularCorrection = false, bool continuousEnergyLoss = false, double limit = 0.1);
		VacuumCherenkov(int id, ref_ptr<Kinematics> kinOt, ref_ptr<Kinematics> kinPh, VacuumCherenkovSpectrum spec = VacuumCherenkovSpectrum::Default, bool havePhotons = true, bool angularCorrection = false, bool continuousEnergyLoss = false, double limit = 0.1);
		VacuumCherenkov(int id, ref_ptr<Kinematics> kin, VacuumCherenkovSpectrum spec = VacuumCherenkovSpectrum::Default, bool havePhotons = true, bool angularCorrection = false, bool continuousEnergyLoss = false, double limit = 0.1);
		void setParticle(int id);
		void setAngularCorrection(bool correction);
		void setContinuousEnergyLoss(bool loss);
		void setHavePhotons(bool photons);
		void setLimit(double limit);
		void setInteractionTag(string tag);
		void setSpectrum(VacuumCherenkovSpectrum spec);
		void setSampler(ref_ptr<Sampler> sampler);
		int getParticle() const;
		string getInteractionTag() const;
		ref_ptr<Kinematics> getKinematicsParticle() const;
		ref_ptr<Kinematics> getKinematicsPhoton() const;
		ref_ptr<Histogram1D> getDistribution() const;
		ref_ptr<Sampler> getSampler() const;
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
		template<class KO, class KP> ref_ptr<Histogram1D> buildSpectrum(const KO& kinOt, const KP& kinPh);
		template<> ref_ptr<Histogram1D> buildSpectrum(const ref_ptr<Kinematics>& kinOt, const ref_ptr<Kinematics>& kinPh);
		template<class KP> ref_ptr<Histogram1D> buildSpectrum(const MonochromaticLorentzViolatingKinematics<2>& kinOt, const KP& kinPh);
		template<class KO, class KP> static double interactionRate(const double& p, const KO& kinOt, const KP& kinPh);
		template<> static double interactionRate(const double& p, const ref_ptr<Kinematics>& kinOt, const ref_ptr<Kinematics>& kinPh);
		template<> static double interactionRate(const double& p, const MonochromaticLorentzViolatingKinematics<2>& kinOt, const MonochromaticLorentzViolatingKinematics<2>& kinPh);
		template<> static double interactionRate(const double& p, const MonochromaticLorentzViolatingKinematics<2>& kinOt, const SpecialRelativisticKinematics& kinPh);
		template<class KO, class KP> static std::pair<double, double> xRange(const KO& kinOt, const KP& kinPh); 
		template<> std::pair<double, double> xRange(const ref_ptr<Kinematics>& kinOt, const ref_ptr<Kinematics>& kinPh);
		template<> static std::pair<double, double> xRange(const MonochromaticLorentzViolatingKinematics<2>& kinOt, const MonochromaticLorentzViolatingKinematics<2>& kinPh);
		static double _Gp(const double& chiOt, const double& chiPh);
		static double _Gm(const double& chiOt, const double& chiPh);
		static double _G0(const double& chiOt, const double& chiPh);
};
/** @}*/

// SWIG sometimes uses weird names; this exposes the right names to the Python interface
constexpr VacuumCherenkovSpectrum VacuumCherenkovSpectrumDefault = VacuumCherenkovSpectrum::Default;
constexpr VacuumCherenkovSpectrum VacuumCherenkovSpectrumStep = VacuumCherenkovSpectrum::Step;
constexpr VacuumCherenkovSpectrum VacuumCherenkovSpectrumFull = VacuumCherenkovSpectrum::Full;



} // namespace livpropa

#endif // LIVPROPA_VACUUMCHERENKOV_H
