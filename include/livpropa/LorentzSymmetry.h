#ifndef LIVPROPA_LORENTZSYMMETRY_H
#define LIVPROPA_LORENTZSYMMETRY_H

#include <cstdlib>
#include <cstring>
#include <fstream>


#include <crpropa/Common.h>
#include <crpropa/Units.h>

#include "livpropa/UnitsAndConstants.h"


using crpropa::pow_integer;


namespace livpropa {


/**
 @class LorentzSymmetry
 @brief Class holding information about a Lorentz symmetry scenario.
 This is done for individual particle types.
 This class corresponds to a general form of the dispersion relation:
   E^2 = (mc)^2 + (pc)^2 + chi * (pc / E_Pl)^n
 */
class LorentzSymmetry {
	protected:
		unsigned int order;
		double coefficient;
		int particleId;

	public:
		LorentzSymmetry();
		LorentzSymmetry(int particle, unsigned int order, double parameter);
		void setOrder(unsigned int order);
		void setCoefficient(double parameter);
		void setParticleId(int id);
		unsigned int getOrder() const;
		double getCoefficient() const;
		int getParticleId() const;
		double getEnergyScale() const;
		bool isSuperluminal() const;
		bool isSubluminal() const;\
		double getSymmetryBreakingShift(const double& p) const;
		double computeEnergyFromMomentum(const double& p) const;

};


} // namespace livpropa

#endif // LIVPROPA_LORENTZSYMMETRY_H
