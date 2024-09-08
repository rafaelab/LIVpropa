#ifndef LIVPROPA_KINEMATICS_H
#define LIVPROPA_KINEMATICS_H

#include <algorithm>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <unordered_map>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/Polynomials>

#include "livpropa/Common.h"
#include "livpropa/UnitsAndConstants.h"



namespace livpropa {


/**
 @class Kinematics
 @brief Class holding information about the kinematics at play.
  The class is abstract.
  Note that the virtual methods are explicitly implemented with dummy values to workaround a SWIG bug.
 */
class AbstractKinematics: public crpropa::Referenced {
	public:
		virtual ~AbstractKinematics() = default;
		virtual string getShortIdentifier() const {
			return "KinematicsBase";
		};
		virtual string getLocationData(const vector<int>& particles) const {
			return "";
		};
		virtual double getSymmetryBreakingShift(const double& p) const {
			return 0;
		};
		virtual double computeMomentumFromEnergy(const double& E, const int& id) const {
			return 0;
		};
		virtual double computeEnergy2FromMomentum(const double& p, const int& id) const {
			return 0;
		};
		double computeEnergyFromMomentum(const double& p, const int& id) const;
		
};


/**
 @class SpecialRelativity
 @brief Class holding information about the usual Lorentz invariant case.
  This class corresponds to a general form of the dispersion relation:
    E^2 = (mc)^2 + (pc)^2
 */
class SpecialRelativity : public AbstractKinematics {
	public:
		SpecialRelativity();
		~SpecialRelativity();
		string getShortIdentifier() const;
		string getLocationData(const vector<int>& particles) const;
		double getSymmetryBreakingShift(const double& p, const int& id) const;
		double computeEnergy2FromMomentum(const double& p, const int& id) const;
		double computeMomentumFromEnergy(const double& E, const int& id) const;
};



/**
 @class LorentzViolating
 @brief Class holding information about the a Lorentz-breaking kinematics.
  The class is abstract.

  Sometimes we need to solve the equation for p and there are multiple real positive solutions.
  In this case, one of the following strategies can be adopted:
    (1) randomly pick one of the possible values;
	(2) compute the average;
	(3) the closest to the usual special-relativistic result;
	(4) the lowest among the values;
	(5) the largest among the values.
  Note that 1 and 2 are the most "sensible" ones.
  The third is phenomenologically motivated.
 */
class LorentzViolating : public AbstractKinematics {
	public:
		enum class SymmetryBreaking {
			Random,
			Average,
			Smallest,
			Largest,
			Closest
		};

	protected:
		SymmetryBreaking symmetryBreaking;
		static double selectFinalMomentumRandom(const vector<double>& ps);
		static double selectFinalMomentumSmallest(const vector<double>& ps);
		static double selectFinalMomentumLargest(const vector<double>& ps);
		static double selectFinalMomentumAverage(const vector<double>& ps);
		static double selectFinalMomentumClosest(const vector<double>& ps, const double& E, const int& id);

	public:
		virtual ~LorentzViolating() = default;
		void setSymmetryBreaking(SymmetryBreaking treatment);
		SymmetryBreaking getSymmetryBreaking() const;
		double selectFinalMomentum(const vector<double>& ps, const double& E, const int& id) const;
};


/**
 @class LorentzViolatingMonochromatic
 @brief Class holding information about a scenario with monochromatic LIV.
  A simple phenomenological implementation of LIV.
  This preserves energy-momentum conservations and modifies the dispersion relations as:
	E^2 = m^2 + p^2 + chi (pc / E_pl)^n.
  The fact that only a single value of n is taken into account defines the naming choice "monochromatic".
 */
class LorentzViolatingMonochromatic : public LorentzViolating {
	protected:
		unsigned int order; 
		double coefficient;
		static vector<double> computeMomentumFromEnergy0(const double& E, const double& m, const double& chi);
		static vector<double> computeMomentumFromEnergy1(const double& E, const double& m, const double& chi);
		static vector<double> computeMomentumFromEnergy2(const double& E, const double& m, const double& chi);
		static vector<double> computeMomentumFromEnergyN(const double& E, const double& m, const double& chi, const unsigned int& n);

	public:
		LorentzViolatingMonochromatic(SymmetryBreaking symmetryBreaking = SymmetryBreaking::Random);
		LorentzViolatingMonochromatic(unsigned int n, SymmetryBreaking symmetryBreaking = SymmetryBreaking::Random);
		LorentzViolatingMonochromatic(unsigned int n, double chi, SymmetryBreaking symmetryBreaking = SymmetryBreaking::Random);
		~LorentzViolatingMonochromatic();
		void setOrder(unsigned int n);
		void setCoefficient(double coefficient);
		unsigned int getOrder() const;
		double getCoefficient() const;
		string getShortIdentifier() const;
		double getSymmetryBreakingShift(const double& p) const;
		double computeEnergy2FromMomentum(const double& p, const int& id) const;
		double computeMomentumFromEnergy(const double& E, const int& id) const;

};



/**
 @class MonochromaticLIV
 @brief Class holding information about a scenario with monochromatic LIV.
  A simple phenomenological implementation of LIV.
  This preserves energy-momentum conservations and modifies the dispersion relations as:
	E^2 = m^2 + p^2 + chi (pc / E_pl)^n.
  The fact that only a single value of n is taken into account defines the naming choice "monochromatic".

  Sometimes we need to solve the equation for p and there are multiple real positive solutions.
  In this case, one of the following strategies can be adopted:
    (1) randomly pick one of the possible values;
	(2) compute the average;
	(3) the closest to the usual special-relativistic result;
	(4) the lowest among the values;
	(5) the largest among the values.
  Note that 1 and 2 are the most "sensible" ones.
  The third is phenomenologically motivated, but not physically justifiable.
 */
class MonochromaticLIV : public LorentzViolating {
	protected:
		unsigned int order; 
		unordered_map<int, double> coefficients = {};
		SymmetryBreaking symmetryBreaking;
		using CoefficientsIterator = typename std::unordered_map<int, double>::const_iterator;

	public:
		MonochromaticLIV(SymmetryBreaking symmetryBreaking = SymmetryBreaking::Random);
		MonochromaticLIV(unsigned int n, SymmetryBreaking symmetryBreaking = SymmetryBreaking::Random);
		MonochromaticLIV(unsigned int n, double chi, SymmetryBreaking symmetryBreaking = SymmetryBreaking::Random);
		MonochromaticLIV(unsigned int order, unordered_map<int, double> coeff, SymmetryBreaking symmetryBreaking = SymmetryBreaking::Random);
		MonochromaticLIV(unsigned int order, vector<int> particles, vector<double> chi, SymmetryBreaking symmetryBreaking = SymmetryBreaking::Random);
		~MonochromaticLIV();
		void setOrder(unsigned int n);
		// void setSymmetryBreaking(SymmetryBreaking treatment);
		void setCoefficients(unordered_map<int, double> coeffs);
		void addCoefficient(int particle, double coeff);
		unsigned int getOrder() const;
		std::unordered_map<int, double> getCoefficients() const;
		std::string getShortIdentifier() const;
		std::string getLocationData(const std::vector<int>& particles) const;
		std::vector<int> getParticles()  const;
		double getCoefficientForParticle(const int& particle) const;
		double getSymmetryBreakingShift(const double& p, const int& id) const;
		double computeEnergy2FromMomentum(const double& p, const int& id) const;
		double computeMomentumFromEnergy(const double& E, const int& id) const;

};



// class Kinematics {
// 	protected:
// 		// std::ref_ptr<AbstractKinematics> kinematics;

// 	public:
// 		Kinematics();

// };





} // namespace livpropa

#endif // LIVPROPA_KINEMATICS_H
