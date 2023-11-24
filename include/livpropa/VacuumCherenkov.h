#ifndef LIVPROPA_VACUUMCHERENKOV_H
#define LIVPROPA_VACUUMCHERENKOV_H

#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <unordered_map>

#include <crpropa/Candidate.h>
#include <crpropa/Common.h>
#include <crpropa/Module.h>
#include <crpropa/Random.h>
#include <crpropa/Referenced.h>
#include <crpropa/Units.h>
#include <crpropa/Vector3.h>

#include "livpropa/Data.h"
#include "livpropa/Kinematics.h"
#include "livpropa/UnitsAndConstants.h"

using crpropa::Candidate;
using crpropa::Module;
using crpropa::Random;
using crpropa::Vector3d;
using crpropa::pow_integer;
using crpropa::ref_ptr;


namespace livpropa {

/**
 @class VacuumCherenkov
 @brief Simulate the emission of vacuum Cherenkov radiation due to LIV effects.

 This implementation follows the description from:
 	arXiv:2311.XXXXX
 Exactly one Cherenkov photon is emitted at a time, with an energy equal to the difference between the primary's energy and the threshold energy.
 The emission is consider to be instantaneous (this might change in the future).
 While the structure here is meant to be generic, it is likely applicable only to electrons.
*/
class VacuumCherenkov: public Module {
	private:
		ref_ptr<Kinematics> kinematics;
		bool havePhotons;
		double limit;
		double thinning;
		std::string interactionTag = "VC";

	public:
		VacuumCherenkov(ref_ptr<Kinematics> kinematics, bool havePhotons = false, double thinning = 0, double limit = 0.1);
		void setKinematics(ref_ptr<Kinematics> kinematics);
		void setHavePhotons(bool photons);
		void setLimit(double limit);
		void setThinning(double thinning);
		void setInteractionTag(std::string tag);
		std::string getInteractionTag() const;
		double computeThresholdEnergy() const;
		void process(Candidate* candidate) const;
};
/** @}*/

} // namespace livpropa

#endif // LIVPROPA_VACUUMCHERENKOV_H
