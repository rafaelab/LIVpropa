#ifndef LIVPROPA_PAIRPRODUCTION_H
#define LIVPROPA_PAIRPRODUCTION_H

#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>

#include <crpropa/Candidate.h>
#include <crpropa/Common.h>
#include <crpropa/Module.h>
#include <crpropa/PhotonBackground.h>
#include <crpropa/Random.h>
#include <crpropa/Referenced.h>
#include <crpropa/Units.h>
#include <crpropa/Vector3.h>

#include "livpropa/Data.h"
#include "livpropa/Kinematics.h"
#include "livpropa/UnitsAndConstants.h"

using crpropa::Candidate;
using crpropa::Module;
using crpropa::PhotonField;
using crpropa::Random;
using crpropa::Referenced;
using crpropa::Vector3d;
using crpropa::ref_ptr;
using crpropa::pow_integer;
using crpropa::interpolate;
using crpropa::closestIndex;


namespace livpropa {


/**
 @class PairProduction
 @brief Breit-Wheeler pair production of photons with background photons considering Lorentz invariance violation (LIV).

 This module simulates electron-pair production of photons with background photons:
 	gamma + gamma_b --> e+ + e- 
 The resulting electron-positron pair is optionally created (default = false).
 The module limits the propagation step size to a fraction of the mean free path (default = 0.1).
 Thinning is available. A thinning of 0 means that all particles are tracked. 
 For the maximum thinning of 1, only a few representative particles are added to the list of secondaries.
 Note that for thinning>0 the output must contain the column "weights", which should be included in the post-processing.
 */
class PairProduction: public Module {
	protected:
		ref_ptr<PhotonField> photonField; 	// target photon field
		ref_ptr<Kinematics> kinematics; //!< the type of kinematics (SR, LIV)
		bool haveElectrons; // add secondary electrons to simulation
		double limit; // limit the step to a fraction of the mean free path
		double thinning; // factor of the thinning (0: no thinning, 1: maximum thinning)
		std::string interactionTag;  //!< tag corresponding to this interaction
		std::vector<double> tabEnergy;  //!< electron energy in [J]
		std::vector<double> tabRate;  //!< interaction rate in [1/m]
		std::vector<double> tabE;  //!< electron energy in [J]
		std::vector<double> tabs;  //!< s_kin = s - m^2 in [J**2]
		std::vector<std::vector<double>> tabCDF;  //!< cumulative interaction rate

	public:
		PairProduction(ref_ptr<PhotonField> photonField, ref_ptr<Kinematics> kinematics, bool haveElectrons = false, double thinning = 0, double limit = 0.1);
		void setPhotonField(ref_ptr<PhotonField> photonField);
		void setHaveElectrons(bool haveElectrons);
		void setLimit(double limit);
		void setThinning(double thinning);
		void setKinematics(ref_ptr<Kinematics> kin);
		void setInteractionTag(std::string tag);
		std::string getInteractionTag() const;
		void initRate(std::string filename);
		void initCumulativeRate(std::string filename);
		void performInteraction(Candidate* candidate) const;
		void process(Candidate* candidate) const;
};


} // namespace livpropa

#endif // LIVPROPA_PAIRPRODUCTION_H
