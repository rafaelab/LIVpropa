#ifndef LIVPROPA_KINEMATICS_H
#define LIVPROPA_KINEMATICS_H

#include <algorithm>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/Polynomials>

#include <kiss/logger.h>

#include "livpropa/Common.h"
#include "livpropa/UnitsAndConstants.h"



namespace livpropa {


// forward declarations
class SpecialRelativisticKinematics;
class LorentzViolatingKinematics;
class AbstractMonochromaticLorentzViolatingKinematics;
template<int N> class MonochromaticLorentzViolatingKinematics;



/**
 @class Kinematics
 @brief Class holding information about the kinematics at play.
  The class is abstract.
  Note that the virtual methods are explicitly implemented with dummy values to workaround a SWIG bug.
 */
class Kinematics: public crpropa::Referenced {
	public:
		virtual ~Kinematics() = default;
		virtual string getNameTag() const {
			return "Kinematics";
		};
		virtual string getIdentifier() const {
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
			return "Kinematics";
		};
		double computeEnergyFromMomentum(const double& p, const int& id) const;
		bool isLorentzInvariant() const;
		bool isLorentzViolating() const;
		bool isSpecialRelativistic() const;
		bool isLorentzViolatingMonochromatic() const;
		const SpecialRelativisticKinematics& toSpecialRelativisticKinematics() const;
		const MonochromaticLorentzViolatingKinematics<0>& toMonochromaticLorentzViolatingKinematics0() const;
		const MonochromaticLorentzViolatingKinematics<1>& toMonochromaticLorentzViolatingKinematics1() const;
		const MonochromaticLorentzViolatingKinematics<2>& toMonochromaticLorentzViolatingKinematics2() const;
		friend std::ostream& operator<<(std::ostream& os, const Kinematics& kin);
};

using ParticleKinematicsDict = unordered_map<int, ref_ptr<Kinematics>>;
using ParticleKinematicsDictIterator = typename ParticleKinematicsDict::const_iterator;


/**
 @class SpecialRelativisticKinematics
 @brief Class holding information about the usual Lorentz invariant case.
  This class corresponds to a general form of the dispersion relation:
    E^2 = (mc)^2 + (pc)^2
 */
class SpecialRelativisticKinematics : public Kinematics {
	public:
		SpecialRelativisticKinematics();
		~SpecialRelativisticKinematics();
		string getNameTag() const;
		string getIdentifier() const;
		double getSymmetryBreakingShift(const double& p, const int& id) const;
		double getCoefficient() const;
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
class LorentzViolatingKinematics : public Kinematics {
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
		virtual string getNameTag() const {
			return "LorentzViolatingKinematics";
		};
		virtual string getIdentifier() const {
			return "";
		};
		virtual string info() const {
			return "LorentzViolatingKinematics";
		};
		void setSymmetryBreaking(SymmetryBreaking treatment);
		SymmetryBreaking getSymmetryBreaking() const;
		double selectFinalMomentum(const vector<double>& ps, const double& E, const int& id) const;
};



/**
 @class AbstractLorentzViolatingKinematicsMonochromatic
 @brief Class holding information about a scenario with monochromatic LIV.
  A simple phenomenological implementation of LIV.
  This preserves energy-momentum conservations and modifies the dispersion relations as:
	E^2 = m^2 + p^2 + chi (pc / E_pl)^n.
  The fact that only a single value of n is taken into account defines the naming choice "monochromatic".
 */
class AbstractMonochromaticLorentzViolatingKinematics : public LorentzViolatingKinematics {
	public:
		using LorentzViolatingKinematics::SymmetryBreaking;

	protected:
		int order;
		double coefficient;

	public:
		virtual ~AbstractMonochromaticLorentzViolatingKinematics() = default;
		void setOrder(int order);
		void setCoefficient(double coefficient);
		int getOrder() const;
		double getCoefficient() const;
		string getNameTag() const;
		string getIdentifier() const;
		double getSymmetryBreakingShift(const double& p) const;
		double computeEnergy2FromMomentum(const double& p, const int& id) const;
		string info() const;
		virtual double computeMomentumFromEnergy(const double& E, const int& id) const {
			return 0;
		};
};


template<int N>
class MonochromaticLorentzViolatingKinematics : public AbstractMonochromaticLorentzViolatingKinematics {
	public:
		MonochromaticLorentzViolatingKinematics(SymmetryBreaking symmetryBreaking = SymmetryBreaking::Closest);
		MonochromaticLorentzViolatingKinematics(double coefficient, SymmetryBreaking symmetryBreaking = SymmetryBreaking::Closest);
		~MonochromaticLorentzViolatingKinematics();
		double computeMomentumFromEnergy(const double& E, const int& id) const;
};



/**
 @class KinematicsMap
 @brief Class holding information about the kinematics at play for each type of particle.
  This class is a container for the kinematics of different particles.
  Particles whose kinematics are not specified will be treated as special relativistic.
 */
class KinematicsMap {
	protected:
		ParticleKinematicsDict kinematics;

	public:
		KinematicsMap();
		KinematicsMap(vector<int> p, vector<ref_ptr<Kinematics>> kin);
		KinematicsMap(vector<pair<int, ref_ptr<Kinematics>>>);
		KinematicsMap(vector<int> p, ref_ptr<Kinematics> kin);
		void add(const int& particle, ref_ptr<Kinematics> kin);
		void remove(const int& id);
		bool isLorentzInvariant() const;
		bool isLorentzViolating() const;
		bool exists(const int& pId) const;
		vector<int> getParticles() const;
		string getIdentifierForParticle(const int& pId, bool showParticleId = true) const; 
		string getIdentifier(const std::vector<int>& particles, bool simplify = false) const;
		string info() const;
		ParticleKinematicsDict getParticleKinematicsMap() const;
		const ref_ptr<Kinematics>& find(const int& id, bool showWarningInexistent = true) const;
		const ref_ptr<Kinematics>& operator[](const int& pId);
		ref_ptr<Kinematics> operator[](const int& pId) const;
		friend std::ostream& operator<<(std::ostream& os, const KinematicsMap& kin);
};




} // namespace livpropa

#endif // LIVPROPA_KINEMATICS_H
