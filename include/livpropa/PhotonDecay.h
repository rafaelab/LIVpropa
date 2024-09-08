#ifndef LIVPROPA_PHOTONDECAY_H
#define LIVPROPA_PHOTONDECAY_H

#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <unordered_map>

#include "livpropa/Common.h"
#include "livpropa/Data.h"
#include "livpropa/Kinematics.h"
#include "livpropa/UnitsAndConstants.h"



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
		ref_ptr<AbstractKinematics> kinematics;
		bool haveElectrons;
		double limit;
		double thinning;
		string interactionTag;

	public:
		PhotonDecay(ref_ptr<AbstractKinematics> kinematics, bool haveElectrons = false, double thinning = 0, double limit = 0.1);
		void setKinematics(ref_ptr<AbstractKinematics> kinematics);
		void setHaveElectrons(bool electrons);
		void setLimit(double limit);
		void setThinning(double thinning);
		void setInteractionTag(string tag);
		string getInteractionTag() const;
		double computeThresholdMomentum() const;
		double computeThresholdEnergy() const;
		void process(Candidate* candidate) const;
};
/** @}*/

} // namespace livpropa

#endif // LIVPROPA_PHOTONDECAY_H
