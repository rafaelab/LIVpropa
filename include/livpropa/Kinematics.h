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

#include <kiss/logger.h>

#include "livpropa/Common.h"
#include "livpropa/UnitsAndConstants.h"



namespace livpropa {


// forward declaration
class SpecialRelativisticKinematics;
class LorentzViolatingKinematicsMonochromatic;


/**
 @class Kinematics
 @brief Class holding information about the kinematics at play.
  The class is abstract.
  Note that the virtual methods are explicitly implemented with dummy values to workaround a SWIG bug.
 */
class AbstractKinematics: public crpropa::Referenced {
	public:
		virtual ~AbstractKinematics() = default;
		virtual string getNameTag() const {
			return "AbstractKinematics";
		};
		virtual string getFilenamePart() const {
			return "";
		};
		virtual double getSymmetryBreakingShift(const double& p, const int& id) const {
			return 0;
		};
		virtual double computeMomentumFromEnergy(const double& E, const int& id) const {
			return 0;
		};
		virtual double computeEnergy2FromMomentum(const double& p, const int& id) const {
			return 0;
		};
		virtual string info() const {
			return "AbstractKinematics";
		};
		double computeEnergyFromMomentum(const double& p, const int& id) const;
		const SpecialRelativisticKinematics* toSpecialRelativisticKinematics() const;
		const LorentzViolatingKinematicsMonochromatic* toLorentzViolatingKinematicsMonochromatic() const;
		bool isSpecialRelativistic() const;
		bool isLorentzViolatingMonochromatic() const;
};



/**
 @class SpecialRelativisticKinematics
 @brief Class holding information about the usual Lorentz invariant case.
  This class corresponds to a general form of the dispersion relation:
    E^2 = (mc)^2 + (pc)^2
 */
class SpecialRelativisticKinematics : public AbstractKinematics {
	public:
		SpecialRelativisticKinematics();
		~SpecialRelativisticKinematics();
		string getNameTag() const;
		string getFilenamePart() const;
		// double getCoefficient() const;
		double getSymmetryBreakingShift(const double& p, const int& id) const;
		double computeEnergy2FromMomentum(const double& p, const int& id) const;
		double computeMomentumFromEnergy(const double& E, const int& id) const;
		string info() const;
};


/**
 @class LorentzViolatingKinematics
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
class LorentzViolatingKinematics : public AbstractKinematics {
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
		virtual ~LorentzViolatingKinematics() = default;
		void setSymmetryBreaking(SymmetryBreaking treatment);
		SymmetryBreaking getSymmetryBreaking() const;
		double selectFinalMomentum(const vector<double>& ps, const double& E, const int& id) const;
		string info() const;
};



/**
 @class LorentzViolatingKinematicsMonochromatic
 @brief Class holding information about a scenario with monochromatic LIV.
  A simple phenomenological implementation of LIV.
  This preserves energy-momentum conservations and modifies the dispersion relations as:
	E^2 = m^2 + p^2 + chi (pc / E_pl)^n.
  The fact that only a single value of n is taken into account defines the naming choice "monochromatic".
 */
class LorentzViolatingKinematicsMonochromatic : public LorentzViolatingKinematics {
	protected:
		unsigned int order; 
		double coefficient;
		static vector<double> computeMomentumFromEnergy0(const double& E, const double& m, const double& chi);
		static vector<double> computeMomentumFromEnergy1(const double& E, const double& m, const double& chi);
		static vector<double> computeMomentumFromEnergy2(const double& E, const double& m, const double& chi);
		static vector<double> computeMomentumFromEnergyN(const double& E, const double& m, const double& chi, const unsigned int& n);

	public:
		LorentzViolatingKinematicsMonochromatic(SymmetryBreaking symmetryBreaking = SymmetryBreaking::Closest);
		LorentzViolatingKinematicsMonochromatic(unsigned int n, SymmetryBreaking symmetryBreaking = SymmetryBreaking::Closest);
		LorentzViolatingKinematicsMonochromatic(unsigned int n, double chi, SymmetryBreaking symmetryBreaking = SymmetryBreaking::Closest);
		~LorentzViolatingKinematicsMonochromatic();
		void setOrder(unsigned int n);
		void setCoefficient(double coefficient);
		unsigned int getOrder() const;
		double getCoefficient() const;
		string getNameTag() const;
		string getFilenamePart() const;
		double getSymmetryBreakingShift(const double& p) const;
		double computeEnergy2FromMomentum(const double& p, const int& id) const;
		double computeMomentumFromEnergy(const double& E, const int& id) const;
		string info() const;
};


/**
 @class Kinematics
 @brief Class holding information about the kinematics at play.
  This class is a container for the kinematics of different particles.
 */
class Kinematics {
	public:
		typedef unordered_map<int, ref_ptr<AbstractKinematics>> ParticleKinematicsMap;

	private:
		using ParticleKinematicsIterator = typename ParticleKinematicsMap::const_iterator;
		ref_ptr<AbstractKinematics> specialRelativity;

	protected:
		ParticleKinematicsMap kinematics;

	public:
		Kinematics();
		Kinematics(vector<int> p, vector<ref_ptr<AbstractKinematics>> kin);
		Kinematics(vector<pair<int, ref_ptr<AbstractKinematics>>>);
		Kinematics(vector<int> p, ref_ptr<AbstractKinematics> kin);
		void add(const int& particle, const ref_ptr<AbstractKinematics>& kin);
		void remove(const int& id);
		bool isLorentzInvariant() const;
		bool isLorentzViolatingKinematics() const;
		bool exists(const int& pId) const;
		vector<int> getParticles() const;
		string getIdentifierForParticle(const int& pId, bool showParticleId = true) const; 
		string getIdentifier(const std::vector<int>& particles, bool simplify = false) const;
		ParticleKinematicsMap getParticleKinematicsMap() const;
		const ref_ptr<AbstractKinematics>& find(const int& id, bool showWarningInexistent = true) const;
		const ref_ptr<AbstractKinematics>& operator[](const int& pId);
		ref_ptr<AbstractKinematics> operator[](const int& pId) const;
};


/** 
This function prints Kinematics-type objects in a nice way
*/
std::ostream& operator<<(std::ostream& os, const AbstractKinematics& kin);
std::ostream& operator<<(std::ostream& os, const Kinematics& kin);




} // namespace livpropa

#endif // LIVPROPA_KINEMATICS_H
