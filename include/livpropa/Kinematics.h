#ifndef LIVPROPA_KINEMATICS_H
#define LIVPROPA_KINEMATICS_H

#include <algorithm>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <unordered_map>

#include <crpropa/Common.h>
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
 */
class Kinematics : public Referenced {
	public:
		virtual ~Kinematics() = default;
		virtual std::string getShortIdentifier() const = 0;
		virtual std::string getLocationData(std::vector<int> particles = {}) const = 0;
		virtual double getSymmetryBreakingShift(const double& p) const = 0;
		virtual double computeMomentumFromEnergy(const double& E, const int& id) const = 0;
		double computeEnergy2FromMomentum(const double& p, const int& id) const;
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
		std::string getShortIdentifier() const;
		std::string getLocationData(std::vector<int> particles = {}) const;
		double getSymmetryBreakingShift(const double& p, const int& id) const;
		// double computeEnergy2FromMomentum(const double& p, const int& id) const;
		// double computeEnergyFromMomentum(const double& p, const int& id) const;
		double computeMomentumFromEnergy(const double& E, const int& id) const;
};


/**
 @class MonochromaticLIV
 @brief Class holding information about a scenario with monochromatic LIV.
  A simple phenomenological implementation of LIV.
  This preserves energy-momentum conservations and modifies the dispersion relations as:
	E^2 = m^2 + p^2 + f^n.
  The fact that only a single value of n is taken into account defines the naming choice "monochromatic".
 */
class MonochromaticLIV : public Kinematics {
	protected:
		unsigned int order; 
		std::unordered_map<int, double> coefficients;
		using CoefficientsIterator = typename std::unordered_map<int, double>::const_iterator;

	public:
		MonochromaticLIV();
		MonochromaticLIV(unsigned int n);
		MonochromaticLIV(unsigned int order, std::unordered_map<int, double> coeff);
		MonochromaticLIV(unsigned int order, std::vector<int> particles, std::vector<double> chi);
		void setOrder(unsigned int n);
		void setCoefficients(std::unordered_map<int, double> coeffs);
		void addCoefficient(int particle, double coeff);
		unsigned int getOrder() const;
		std::unordered_map<int, double> getCoefficients() const;
		std::string getShortIdentifier() const;
		std::string getLocationData(std::vector<int> particles = {}) const;
		std::vector<int> getParticles()  const;
		double getCoefficientForParticle(const int& particle) const;
		double getSymmetryBreakingShift(const double& p, const int& id) const;
		double computeMomentumFromEnergy(const double& E, const int& id) const;
};



} // namespace livpropa

#endif // LIVPROPA_KINEMATICS_H
