#ifndef LIVPROPA_PHOTONDECAY_H
#define LIVPROPA_PHOTONDECAY_H

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
using crpropa::Referenced;
using crpropa::Vector3d;
using crpropa::pow_integer;
using crpropa::ref_ptr;


namespace livpropa {

/**
 @class PhotonDecay
 @brief Simulate the decay of a photon into a pair, due to LIV effects.

 This implementation follows the description from:
 	arXiv:2312.XXXXX
 The emission is considered to be instantaneous (this might change in the future).
*/
class PhotonDecay: public Module {
	private:
		ref_ptr<Kinematics> kinematics;
		bool haveElectrons;
		double limit;
		double thinning;
		std::string interactionTag;

	public:
		PhotonDecay(ref_ptr<Kinematics> kinematics, bool haveElectrons = false, double thinning = 0, double limit = 0.1);
		void setKinematics(ref_ptr<Kinematics> kinematics);
		void setHaveElectrons(bool electrons);
		void setLimit(double limit);
		void setThinning(double thinning);
		void setInteractionTag(std::string tag);
		std::string getInteractionTag() const;
		double computeThresholdMomentum() const;
		double computeThresholdEnergy() const;
		void process(Candidate* candidate) const;
};
/** @}*/

} // namespace livpropa

#endif // LIVPROPA_PHOTONDECAY_H
