#ifndef LIVPROPA_KINEMATICS_H
#define LIVPROPA_KINEMATICS_H

#include <algorithm>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <unordered_map>

#include <crpropa/Common.h>
#include <crpropa/Random.h>
#include <crpropa/Referenced.h>
#include <crpropa/Units.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/Polynomials>

#include "livpropa/UnitsAndConstants.h"


using crpropa::Referenced;
using crpropa::pow_integer;
using crpropa::ref_ptr;



namespace livpropa {

/**
 @class Kinematics
 @brief Class holding information about the kinematics at play.
  The class is abstract.
  Note that the virtual methods are explicitly implemented with dummy values to workaround a SWIG bug.
 */
class Kinematics : public Referenced {
	public:
		virtual ~Kinematics() = default;
		virtual std::string getShortIdentifier() const {
			return "KinematicsBase";
		};
		virtual std::string getLocationData(const std::vector<int>& particles) const {
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
class SpecialRelativity : public Kinematics {
	public:
		SpecialRelativity();
		~SpecialRelativity();
		std::string getShortIdentifier() const;
		std::string getLocationData(const std::vector<int>& particles) const;
		double getSymmetryBreakingShift(const double& p, const int& id) const;
		double computeEnergy2FromMomentum(const double& p, const int& id) const;
		double computeMomentumFromEnergy(const double& E, const int& id) const;
};

/**
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
enum class SymmetryBreaking {
	Random,
	Average,
	Smallest,
	Largest,
	Closest
};


/**
 @class MonochromaticLIV
 @brief Class holding information about a scenario with monochromatic LIV.
  A simple phenomenological implementation of LIV.
  This preserves energy-momentum conservations and modifies the dispersion relations as:
	E^2 = m^2 + p^2 + chi (pc / E_pl)^n.
  The fact that only a single value of n is taken into account defines the naming choice "monochromatic".
 */
class MonochromaticLIV : public Kinematics {
	protected:
		unsigned int order; 
		std::unordered_map<int, double> coefficients = {};
		SymmetryBreaking symmetryBreaking;
		using CoefficientsIterator = typename std::unordered_map<int, double>::const_iterator;

	public:
		MonochromaticLIV(SymmetryBreaking symmetryBreaking = SymmetryBreaking::Random);
		MonochromaticLIV(unsigned int n, SymmetryBreaking symmetryBreaking = SymmetryBreaking::Random);
		MonochromaticLIV(unsigned int n, double chi, SymmetryBreaking symmetryBreaking = SymmetryBreaking::Random);
		MonochromaticLIV(unsigned int order, std::unordered_map<int, double> coeff, SymmetryBreaking symmetryBreaking = SymmetryBreaking::Random);
		MonochromaticLIV(unsigned int order, std::vector<int> particles, std::vector<double> chi, SymmetryBreaking symmetryBreaking = SymmetryBreaking::Random);
		~MonochromaticLIV();
		void setOrder(unsigned int n);
		void setSymmetryBreaking(SymmetryBreaking treatment);
		void setCoefficients(std::unordered_map<int, double> coeffs);
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



} // namespace livpropa

#endif // LIVPROPA_KINEMATICS_H
